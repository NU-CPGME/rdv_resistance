#Analysis of Intra-Host Variation 


library(readxl)
library(lubridate)
library(ggplot2)
library(tidyr)
library(ggsci)
library(ggpubr)
library(reshape2)
library('tidyverse')
library(dplyr)
library(gtsummary)
library(MASS)
library(car)
library(gt)
library(OddsPlotty)
library(effects)
library(lme4)
library(emmeans)
library(lmerTest)



#Load and format Metadata

#RDV dataset

Metadata_RDV<-read.csv("hosp.csv", stringsAsFactors = T)
colnames(Metadata_RDV)
Metadata_RDV<-Metadata_RDV%>%dplyr::select(-X)
Metadata_RDV$date<-ymd(Metadata_RDV$sample_date)
Metadata_RDV$rdvdate<-ymd(Metadata_RDV$rdvdate)
Metadata_RDV$ca_hosp_admit_dt<-ymd(Metadata_RDV$ca_hosp_admit_dt)
Metadata_RDV$RDV<-plyr::mapvalues(Metadata_RDV$RDV,c("pre","post"),c("Pre","Post"))
Metadata_RDV$RDV<-relevel(Metadata_RDV$RDV,"Pre")
Metadata_RDV$TimeToRDV<-Metadata_RDV$daysdiff_sample_rdv
Metadata_RDV$TimeToHosp<-Metadata_RDV$daysdiff_sample_hosp
Metadata_RDV$mrn<-as.factor(Metadata_RDV$mrn)


#Control dataset
Metadata_Control<-read.csv("control_hosp.csv", stringsAsFactors = T)
colnames(Metadata_Control)
Metadata_Control<-Metadata_Control%>%dplyr::select(-X)
Metadata_Control$date<-ymd(Metadata_Control$sample_date)
Metadata_Control$rdvdate<-ymd(Metadata_Control$rdvdate)
Metadata_Control$ca_hosp_admit_dt<-ymd(Metadata_Control$ca_hosp_admit_dt)
Metadata_Control$Hosp<-plyr::mapvalues(Metadata_Control$Hosp,c("pre","post"),c("Pre","Post"))
Metadata_Control$Hosp<-relevel(Metadata_Control$Hosp,"Pre")
Metadata_Control$TimeToRDV<-Metadata_Control$daysdiff_sample_rdv
Metadata_Control$TimeToHosp<-Metadata_Control$daysdiff_sample_hosp
Metadata_Control$mrn<-as.factor(Metadata_Control$mrn)


#Get Variants from iVAR output

#Variants

myfiles = list.files(path="PathToDirectoryWith_iVAROutputs", pattern="*.tsv", full.names=TRUE)
dat_txt = plyr::ldply(myfiles, read.table, sep = "\t", fill=TRUE, header = TRUE)

myfilesFrame<-as.data.frame(myfiles)
myfilesFrame$name<-myfilesFrame$myfiles
myfilesFrame<-separate(myfilesFrame,col = name,into = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"name"),sep = "/")
myfilesFrame<-separate(myfilesFrame,col = name,into = c("name"),sep = "\\.")


files <- list.files(path = "PathToDirectoryWith_iVAROutputs", pattern = '.tsv$', full.names = T) %>%
  map(read.table, sep = "\t", fill=TRUE, header = TRUE)

fileNames<-myfilesFrame$name

samples=c()
for (i in seq_along(fileNames)) {
  value=rep(fileNames[[i]], length(files[[i]]$REGION))
  samples=append(samples,value)
}

dat_txt$Samples<-samples

dat_txt$study_id<-dat_txt$Sample

dat_txt<-merge(dat_txt,Metadata_RDV,by="study_id")

dat_txt_all<-dat_txt

dat_txt<-dat_txt%>%filter(PASS=="TRUE")%>%droplevels()

dat_txt_Positions<-dat_txt%>%filter(ALT%in%c("A","C","T","G"))

dat_txt_Positions_Counts<-as.data.frame(table(dat_txt_Positions$POS,dat_txt_Positions$RDV))

names(dat_txt_Positions_Counts)<-c("POS","RDV","Count")

dat_txt_Positions_Counts_Spread<-pivot_wider(dat_txt_Positions_Counts,names_from = RDV,values_from = Count)

#If Paired Data Add Time Difference Between Samples
dat_txt_Positions_Counts_Spread$Diff_Post_Pre<-dat_txt_Positions_Counts_Spread$Post-dat_txt_Positions_Counts_Spread$Pre


#Also If Paired, select only 1 of each pre and post samples 

All_Paired_Pre<-dat_txt_all%>%filter(RDV=="Pre"&mrn%in%unique(demo_paired$mrn))
All_Paired_Post<-dat_txt_all%>%filter(RDV=="Post"&mrn%in%unique(demo_paired$mrn))

All_Paired_Pre_Last<-All_Paired_Pre %>%
  group_by(mrn) %>%
  arrange(desc(sample_date)) %>%
  filter(row_number()==1)

All_Paired_Post_First<-All_Paired_Post %>%
  group_by(mrn) %>%
  arrange(sample_date) %>%
  filter(row_number()==1)


All_Paired_Pre_Post_First<-rbind(All_Paired_Pre_Last,All_Paired_Post_First)


All_Paired<-dat_txt_Positions%>%filter(study_id%in%All_Paired_Pre_Post_First$study_id)


All_Paired_Samples_Per_RDV<-All_Paired%>%dplyr::select(mrn,RDV)%>%distinct()

table(All_Paired_Samples_Per_RDV$mrn)

#Filter only positions that pass iVAR call and do not have ambigous positions or indels

All_Paired_Filtered<-All_Paired

All_Paired_Filtered<-All_Paired_Filtered%>%filter(PASS=="TRUE")%>%droplevels()

All_Paired_Filtered_Positions<-All_Paired_Filtered%>%filter(ALT%in%c("A","C","T","G"))%>%droplevels()


All_Paired_Filtered_Positions%>%
  ggplot(aes(x=POS,y=ALT_FREQ,color=ALT))+
  geom_point()+
  theme_pubr()+scale_color_lancet()+facet_grid(.~RDV)+ylab("Mutation Frequency Compared with Wu-1")+ guides(color=guide_legend(title="Change"))


All_Paired_Filtered_Positions_Appear<-All_Paired_Filtered_Positions%>%dplyr::select(study_id,POS,RDV,mrn)%>%distinct()

All_Paired_Filtered_Positions_Counts<-as.data.frame(table(All_Paired_Filtered_Positions_Appear$POS,All_Paired_Filtered_Positions_Appear$RDV))

names(All_Paired_Filtered_Positions_Counts)<-c("POS","RDV","Count")

All_Paired_Filtered_Positions_Counts_Spread<-pivot_wider(All_Paired_Filtered_Positions_Counts,names_from = RDV,values_from = Count)%>%distinct()

All_Paired_Filtered_Positions_Counts_Spread$Diff_Post_Pre<-All_Paired_Filtered_Positions_Counts_Spread$Post-All_Paired_Filtered_Positions_Counts_Spread$Pre


All_Paired_Filtered_Positions_Counts_Spread%>%filter(abs(Diff_Post_Pre)>=4)%>%
  ggplot(aes(x = POS,y= Diff_Post_Pre))+
  geom_col()+
  theme_pubr()+rotate_x_text(45)



#Calculate Entropy for all paired positions


GenomicPostions<-as.data.frame(read.table("../ORFs_POS.txt",header = T))

genelist=data.frame()

for (i in (1:length(All_Paired_Filtered_Positions$study_id))) {
  Mutation=All_Paired_Filtered_Positions[c(i),]$POS
  gene=(filter(GenomicPostions,S<=Mutation&E>=Mutation))$Label
  gene=ifelse(identical(gene, character(0)),"None",as.character(gene))
  MutationList=data.frame(Mutation,gene)
  genelist=rbind(genelist,MutationList)
  
}


All_Paired_Filtered_Positions$Gene<-genelist$gene

#Entropy per Position

All_Paired_Filtered_Positions_sh<-All_Paired_Filtered_Positions %>%
  dplyr::select(study_id,POS,ALT_FREQ,mrn,date,TimeToRDV,TimeToHosp,RDV,Gene) %>%
  group_by(study_id,POS) %>%
  mutate(Sh = (-(ALT_FREQ)*log2(ALT_FREQ))) 


All_Paired_Filtered_Positions_sh_persite<-aggregate(Sh~study_id+POS+mrn+date+TimeToRDV+TimeToHosp+RDV+Gene,All_Paired_Filtered_Positions_sh,sum)

#Entropy Overall

All_Paired_Filtered_Positions_sh_sum<-aggregate(Sh~study_id+mrn+date+TimeToRDV+TimeToHosp+RDV,All_Paired_Filtered_Positions_sh,sum)

#Entropy PerGene

All_Paired_Filtered_Positions_sh_Gene<-aggregate(Sh~study_id+ALT_FREQ+mrn+date+TimeToRDV+TimeToHosp+RDV+Gene,All_Paired_Filtered_Positions_sh,sum)


All_Paired_Pre_Date<-All_Paired_Filtered_Positions_sh_sum%>%filter(RDV=="Pre")%>%dplyr::select(mrn,date)
All_Paired_Post_Date<-All_Paired_Filtered_Positions_sh_sum%>%filter(RDV=="Post")%>%dplyr::select(mrn,date)


All_Paired_Time<-merge(All_Paired_Pre_Date,All_Paired_Post_Date,by="mrn")

All_Paired_Time$TimeBetweenSamples<-ymd(All_Paired_Time$date.y)-ymd(All_Paired_Time$date.x)


All_Paired_Filtered_Positions_sh_sum<-merge(All_Paired_Filtered_Positions_sh_sum,All_Paired_Time[,c(1,4)],by="mrn")

All_Paired_Filtered_Positions_sh_Gene<-merge(All_Paired_Filtered_Positions_sh_Gene,All_Paired_Time[,c(1,4)],by="mrn")


All_Paired_Filtered_Positions_sh_persite%>%
  ggplot(aes(x=POS,y=Sh,color=RDV))+
  geom_point()+
  theme_pubr()+scale_color_lancet()+rotate_x_text(45)+facet_grid(.~RDV)


All_Paired_Filtered_Positions_sh_sum%>%
  ggplot(aes(x=RDV,y=Sh,color=RDV))+geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_lancet()+rotate_x_text(45)#+stat_compare_means()


All_Paired_Filtered_Positions_sh_Gene%>%
  ggplot(aes(x=RDV,y=Sh,color=RDV))+geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_lancet()+rotate_x_text(45)+facet_grid(.~Gene)+stat_compare_means(aes(label = after_stat(p.signif)))

#Model for Entropy Changes Pre vs Post


M<-lmer(Sh ~ TimeBetweenSamples+RDV+
          (1 | mrn), data = All_Paired_Filtered_Positions_sh_sum)

pairs(lsmeans(M, ~ RDV))
tbl_regression(M, exponentiate = TRUE)


#Model for Entropy Changes Per Gene

M1<-lmer(Sh ~ TimeBetweenSamples+RDV*Gene+
           (1 | mrn), data = All_Paired_Filtered_Positions_sh_Gene)


tbl_regression(M1, exponentiate = TRUE)

tbl1<-as.data.frame(rbind(pairs(lsmeans(M1, ~ RDV | Gene)),adjust="FDR"))
gt(tbl1)
tbl1


#Test for significant changes in any position between pre and post with FDR

#Add Demo Data to Control

demo<-read.csv("demographics.csv",stringsAsFactors = T)
summary(demo)
colnames(demo)
demo$sample_date<-ymd(demo$sample_date)
demo$ca_dob<-ymd(demo$ca_dob)

demo$age.correct<-round((demo$sample_date-demo$ca_dob)/365,0)

demo$Vaccination.Status<-ifelse(demo$ca_vaccine1=="", "Unvaccinated","FullyVaccinated")
demo$Vaccination.Status<-ifelse(demo$ca_vaccine1!=""&demo$ca_vaccine2=="", "PartiallyVaccinated",as.character(demo$Vaccination.Status))


demo_RDV_paired<-demo%>%filter(Sample.Type=="RDV paired")

demo_Control_paired<-demo%>%filter(Sample.Type=="Control paired")


#Select Variables

demo_Simple<-demo%>%dplyr::select(sample_date,mrn,rdvdate,ca_sex,ca_race,ca_ethnicity,ca_bmi,
                                  ca_comorbid_sum,age.correct,ca_icu,ca_death,Vaccination.Status,Sample.Type)


demo_Simple$Group<-ifelse(demo_Simple$Sample.Type%in%c("Control paired","Control post"),"Control","RDV")

demo_Simple$Group<-as.factor(demo_Simple$Group)

demo_Simple<-demo_Simple%>%dplyr::select(-Sample.Type)%>%droplevels()

tbl_summary(
  demo_Simple[,-c(1:3)],
  by = Group, # split table by group
  missing = "no" # don't list missing data separately
) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p() %>% # test for a difference between groups
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() 



model_glm = glm(Group ~ . , family="binomial",data=na.omit(demo_Simple[,-c(1:3)]))
summary(model_glm)

tbl_regression(model_glm, exponentiate = TRUE)%>%bold_p(t = 0.05)%>%
  as_gt() %>%
  gt::tab_header(
    title = gt::md("**Control vs RDV**")
  )

odds_plot(model_glm)


#All Entropy per site with metadata

All_RDV_Paired_Filtered_Positions_Pre<-All_RDV_Paired_Filtered_Positions%>%filter(RDV=="Pre")%>%dplyr::select(mrn,POS,ALT_FREQ)
All_RDV_Paired_Filtered_Positions_Post<-All_RDV_Paired_Filtered_Positions%>%filter(RDV=="Post")%>%dplyr::select(mrn,POS,ALT_FREQ)

All_RDV_Paired_Filtered_Positions_Full<-merge(All_RDV_Paired_Filtered_Positions_Pre,All_RDV_Paired_Filtered_Positions_Post,by=c("mrn","POS"),all = T,suffixes = c("_Pre","_Post"))

All_RDV_Paired_Filtered_Positions_Full[is.na(All_RDV_Paired_Filtered_Positions_Full)] <- 0

All_RDV_Paired_Filtered_Positions_Full<-All_RDV_Paired_Filtered_Positions_Full%>%pivot_longer(ALT_FREQ_Pre:ALT_FREQ_Post,names_to = "RDV_Hosp",values_to = "ALT_FREQ")

All_RDV_Paired_Filtered_Positions_Full$RDV_Hosp<-plyr::mapvalues(All_RDV_Paired_Filtered_Positions_Full$RDV_Hosp,unique(All_RDV_Paired_Filtered_Positions_Full$RDV_Hosp),c("Pre","Post"))

All_Control_Paired_Filtered_Positions_Pre<-All_Control_Paired_Filtered_Positions%>%filter(Hosp=="Pre")%>%dplyr::select(mrn,POS,ALT_FREQ)
All_Control_Paired_Filtered_Positions_Post<-All_Control_Paired_Filtered_Positions%>%filter(Hosp=="Post")%>%dplyr::select(mrn,POS,ALT_FREQ)

All_Control_Paired_Filtered_Positions_Full<-merge(All_Control_Paired_Filtered_Positions_Pre,All_Control_Paired_Filtered_Positions_Post,by=c("mrn","POS"),all = T,suffixes = c("_Pre","_Post"))

All_Control_Paired_Filtered_Positions_Full[is.na(All_Control_Paired_Filtered_Positions_Full)] <- 0

All_Control_Paired_Filtered_Positions_Full<-All_Control_Paired_Filtered_Positions_Full%>%pivot_longer(ALT_FREQ_Pre:ALT_FREQ_Post,names_to = "RDV_Hosp",values_to = "ALT_FREQ")

All_Control_Paired_Filtered_Positions_Full$RDV_Hosp<-plyr::mapvalues(All_Control_Paired_Filtered_Positions_Full$RDV_Hosp,unique(All_Control_Paired_Filtered_Positions_Full$RDV_Hosp),c("Pre","Post"))

table(unique(All_RDV_Paired_Filtered_Positions_Full$POS)%in%unique(All_Control_Paired_Filtered_Positions_Full$POS))

table(unique(All_Control_Paired_Filtered_Positions_Full$POS)%in%unique(All_RDV_Paired_Filtered_Positions_Full$POS))


#RDV with all positions

Meta_All_RDV_Paired<-All_RDV_Paired_Pre_Post_First[,-c(2:20)]

names(Meta_All_RDV_Paired)[35]<-"RDV_Hosp"

Meta_All_Control_Paired<-All_Control_Paired_Pre_Post_Last[,-c(2:20)]
names(Meta_All_Control_Paired)[35]<-"RDV_Hosp"

Meta_All_Paired<-rbind(Meta_All_RDV_Paired,Meta_All_Control_Paired)

All_RDV_Paired_Filtered_Positions_Full<-merge(All_RDV_Paired_Filtered_Positions_Full,Meta_All_Paired,by=c("mrn","RDV_Hosp"),sort = F)

All_RDV_Paired_Filtered_Positions_Full_sh<-All_RDV_Paired_Filtered_Positions_Full %>%
  group_by(study_id,POS) %>%
  mutate(Sh = (-(ALT_FREQ)*log2(ALT_FREQ))) 

All_RDV_Paired_Filtered_Positions_Full_sh[is.na(All_RDV_Paired_Filtered_Positions_Full_sh)] <- 0



genelist2=data.frame()

for (i in (1:length(All_RDV_Paired_Filtered_Positions_Full_sh$study_id))) {
  Mutation=All_RDV_Paired_Filtered_Positions_Full_sh[c(i),]$POS
  gene=(filter(GenomicPostions,S<=Mutation&E>=Mutation))$Label
  gene=ifelse(identical(gene, character(0)),"None",as.character(gene))
  MutationList=data.frame(Mutation,gene)
  genelist2=rbind(genelist2,MutationList)
  
}


All_RDV_Paired_Filtered_Positions_Full_sh$Gene<-genelist2$gene

All_RDV_Paired_Filtered_Positions_Full_sh_persite<-aggregate(Sh~study_id+POS+mrn+date+Sample.Type..when.received.+TimeToHosp+RDV_Hosp+Gene,All_RDV_Paired_Filtered_Positions_Full_sh,sum)



All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta<-merge(All_RDV_Paired_Filtered_Positions_Full_sh_persite,(All_RDV_Paired_Filtered_Positions_Full[,c(5,19,28)]%>%distinct()),by="study_id")

PairedMetadataRDV_Control_Pre<-Meta_All_Paired%>%filter(RDV_Hosp=="Pre")
PairedMetadataRDV_Control_Post<-Meta_All_Paired%>%filter(RDV_Hosp=="Post")

PairedMetadataRDV_Control_Pre_Post<-merge(PairedMetadataRDV_Control_Pre[,c(6,7)],PairedMetadataRDV_Control_Post[,c(6,7)],by="mrn")

PairedMetadataRDV_Control_Pre_Post$TimeBetweenSamples<-ymd(PairedMetadataRDV_Control_Pre_Post$sample_date.y)-ymd(PairedMetadataRDV_Control_Pre_Post$sample_date.x)

All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta<-merge(All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta,PairedMetadataRDV_Control_Pre_Post[,c(1,4)],by="mrn")

All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta<-All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta%>%arrange(POS)

All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta$RDV_Hosp<-relevel(as.factor(All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta$RDV_Hosp),ref = "Pre")

All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta$TimeBetweenSamples<-as.numeric(All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta$TimeBetweenSamples)


#Test all positions with data for more than 3 patients

Pos_pValue_df=data.frame()

for (i in unique(All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta$POS)) {
  PositionData=(filter(All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta,POS==i))
  if (length(PositionData$study_id)>6&sum(PositionData$Sh)>0) {
    model_lmer<-lmer(Sh ~ TimeBetweenSamples+ct_n1+RDV_Hosp++
                       (1 | mrn), data = PositionData)
    pValue=summary(pairs(lsmeans(model_lmer, ~ RDV_Hosp)))$p.value
    POS_pValue=data.frame(POS=i,pValue)
    Pos_pValue_df=rbind(Pos_pValue_df,POS_pValue)
  }
  
}

library(qvalue)

Pos_pValue_df$FDR <- p.adjust(Pos_pValue_df$pValue, method="BH")

All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta%>%filter(POS%in%(filter(Pos_pValue_df,FDR<0.1)$POS))%>%droplevels()%>%
  ggplot(aes(x=RDV_Hosp,y=Sh))+
  geom_boxplot(outlier.shape = NA,varwidth = 0.2)+geom_point(aes(color=as.factor(mrn)),position = position_jitterdodge(dodge.width = 0.5),size=2)+ 
  theme_pubr()+ggtitle("Significant Positions in RDV Group")+rotate_x_text(45)+ylab("Entropy")+ guides(color=guide_legend(title="Patient"))+facet_grid(.~POS)+rremove("legend")


#Additional Tests for significant data (Positions 14960 and 15451)


PositionData=(filter(All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta,POS==14960))

PositionData$TimeBetweenSamples<-as.numeric(PositionData$TimeBetweenSamples)

model_lmer<-lmer(Sh ~ TimeBetweenSamples+ct_n1+clade+RDV_Hosp+Sample.Type..when.received.+
                   (1 | mrn), data = PositionData)

tbl_regression(model_lmer, exponentiate = TRUE)

pairs(lsmeans(model_lmer, ~ RDV_Hosp))


plot(effect(model_lmer,term=c("RDV_Hosp")),ylab="Prob")

PositionData=(filter(All_RDV_Paired_Filtered_Positions_Full_sh_persite_Meta,POS==15451))

PositionData$TimeBetweenSamples<-as.numeric(PositionData$TimeBetweenSamples)

model_lmer<-lmer(Sh ~ TimeBetweenSamples+ct_n1+clade+RDV_Hosp+Sample.Type..when.received.+
                   (1 | mrn), data = PositionData)

tbl_regression(model_lmer, exponentiate = TRUE)

pairs(lsmeans(model_lmer, ~ RDV_Hosp))


plot(effect(model_lmer,term=c("RDV_Hosp")),ylab="Prob")
