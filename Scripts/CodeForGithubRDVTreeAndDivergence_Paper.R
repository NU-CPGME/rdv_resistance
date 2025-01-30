
#Consensus Tree and Divergence


AllPaired_COVIDids<-data.frame(Ids=c(unique(All_RDV_Control_Paired_Filtered_Positions_Full_sh_persite_Meta$study_id)))

write.table(AllPaired_COVIDids,"AllPaired_COVIDids.txt",quote = F,row.names = F,col.names = F)



tree_1 <- read.tree("All_RDV_Control_Paired_Consensus_Aligned_Clean.nwk")

FastaNames_1 <- as.data.frame(tree_1$tip.label)
names(FastaNames_1)[1]<-"taxa"

FastaNames_1$study_id<-FastaNames_1$taxa


AllPaired_mrn<-data.frame(mrn=unique(Meta_All_Paired$mrn))
AllPaired_mrn$Group<-ifelse(AllPaired_mrn$mrn%in%Meta_All_RDV_Paired$mrn,"RDV","Control")

AllPaired_mrn$participant<-rownames(AllPaired_mrn)

AllPaired_mrn$participant<-as.factor(as.character(AllPaired_mrn$participant))



FastaNames_1<-dplyr::full_join(FastaNames_1,Meta_All_Paired,by="study_id")
FastaNames_1<-dplyr::full_join(FastaNames_1,AllPaired_mrn,by="mrn")

row.names(FastaNames_1) <- NULL


p <-ggtree(tree_1, layout="circular") %<+% FastaNames_1 +
  geom_tippoint(size=3,aes(color=interaction(RDV_Hosp,Group)))+scale_color_lancet()
p1<-p+ geom_tiplab(offset = 0.0001,aes(label = participant,angle=angle,show.legend=F))+ geom_treescale(fontsize=4, linesize=1, x=-0.001)#+xlim(-0.001, NA)

colourCount = length(unique(FastaNames_1$clade))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))

p1 + ggnewscale::new_scale_fill()+
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=FastaNames_1$clade),
    width=0.0005,offset = 0.2,show.legend=T)+scale_fill_manual(values = getPalette(colourCount))

p <-ggtree(tree_1) %<+% FastaNames_1 +
  geom_tippoint(size=3,aes(color=interaction(RDV_Hosp,Group)))+scale_color_lancet()
p+theme(legend.position="right")+ geom_tiplab(offset = 0.0001,aes(label = participant,show.legend=F))+ geom_treescale(fontsize=4, linesize=1, x=-0.001)

p <-ggtree(tree_1, layout="circular") %<+% FastaNames_1 +
  geom_tippoint(size=3,aes(color=interaction(RDV_Hosp,Group)))+scale_color_lancet()
p+theme(legend.position="right")+ geom_treescale(fontsize=4, linesize=1, x=-0.001)#+xlim(-0.001, NA)

p <-ggtree(tree_1) %<+% FastaNames_1 +
  geom_tippoint(size=3,aes(color=interaction(RDV_Hosp,Group)))+scale_color_lancet()
p+theme(legend.position="right")

p+geom_tiplab(as_ylab=TRUE, color='firebrick')

p<-ggtree(tree_1,layout = "circular")%<+% FastaNames_1 +
  geom_tippoint(size=3,aes(color=interaction(RDV_Hosp,Group)))+scale_color_lancet()




#Divergence

PairedDivergence<-read.csv("All_RDV_Control_Paired_Consensus_Aligned_Clean_HammingMatrix.csv")


PairedDivergence_df<-PairedDivergence%>%filter(Species.1!=Species.2)
PairedDivergence_df$Patient_1<-PairedDivergence_df$Species.1
PairedDivergence_df$Patient_2<-PairedDivergence_df$Species.2

All_Control_Paired_Samples_Per_ControlGroup<-All_Control_Paired_Samples_Per_Control
All_Control_Paired_Samples_Per_ControlGroup$Group<-"Control"
names(All_Control_Paired_Samples_Per_ControlGroup)[2]<-"RDV_Hosp"

All_RDV_Paired_Samples_Per_RDVGroup<-All_RDV_Paired_Samples_Per_RDV
All_RDV_Paired_Samples_Per_RDVGroup$Group<-"RDV"
names(All_RDV_Paired_Samples_Per_RDVGroup)[2]<-"RDV_Hosp"


PairedMetadataRDV_Control<-rbind(All_Control_Paired_Samples_Per_ControlGroup,All_RDV_Paired_Samples_Per_RDVGroup)


PairedDivergence_df$Patient_1<-plyr::mapvalues(PairedDivergence_df$Patient_1,Metadata_RDV$study_id,Metadata_RDV$mrn)
PairedDivergence_df$Patient_1<-plyr::mapvalues(PairedDivergence_df$Patient_1,Metadata_Control$study_id,Metadata_Control$mrn)
PairedDivergence_df$Patient_2<-plyr::mapvalues(PairedDivergence_df$Patient_2,Metadata_RDV$study_id,Metadata_RDV$mrn)
PairedDivergence_df$Patient_2<-plyr::mapvalues(PairedDivergence_df$Patient_2,Metadata_Control$study_id,Metadata_Control$mrn)

PairedDivergence_df_Pairs<-PairedDivergence_df%>%filter(Patient_1==Patient_2)

PairedDivergence_df_Pairs$Group<-plyr::mapvalues(PairedDivergence_df_Pairs$Patient_1,PairedMetadataRDV_Control$mrn,PairedMetadataRDV_Control$Group)

PairedDivergence_df_Pairs$RDV_Hosp_1<-plyr::mapvalues(PairedDivergence_df_Pairs$Species.1,PairedMetadataRDV_Control$study_id,as.character(PairedMetadataRDV_Control$RDV_Hosp))
PairedDivergence_df_Pairs$RDV_Hosp_2<-plyr::mapvalues(PairedDivergence_df_Pairs$Species.2,PairedMetadataRDV_Control$study_id,as.character(PairedMetadataRDV_Control$RDV_Hosp))


PairedDivergence_df_Pairs<-PairedDivergence_df_Pairs %>%
  group_by(Patient_1) %>%
  filter(row_number()==1)

PairedDivergence_df_Pairs$HammingDistance<-PairedDivergence_df_Pairs$Dist

PairedDivergence_df_Pairs$TimeBetweenSamples<-plyr::mapvalues(PairedDivergence_df_Pairs$Patient_1,All_RDV_Paired_Time$mrn,All_RDV_Paired_Time$TimeBetweenSamples)
PairedDivergence_df_Pairs$TimeBetweenSamples<-plyr::mapvalues(PairedDivergence_df_Pairs$TimeBetweenSamples,All_Control_Paired_Time$mrn,All_Control_Paired_Time$TimeBetweenSamples)
PairedDivergence_df_Pairs$TimeBetweenSamples<-as.numeric(as.character(PairedDivergence_df_Pairs$TimeBetweenSamples))



PairedDivergence_df_Pairs%>%
  ggplot(aes(x=Group,y=HammingDistance,color=Group))+geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_jco()+stat_compare_means()

PairedDivergence_df_Pairs%>%
  ggplot(aes(x=Group,y=Dist,color=Group))+geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_jco()+stat_compare_means()

PairedDivergence_df_Pairs%>%
  ggplot(aes(x=Group,y=log(Dist+0.0001),color=Group))+geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  theme_pubr()+scale_color_jco()+stat_compare_means()

PairedDivergence_df_Pairs%>%
  ggplot(aes(x=Group,y=as.factor(Dist>0),color=Group))+
  geom_count()+scale_size_area(max_size = 20)+
  scale_size(range = c(0,20))+
  theme_pubr()+scale_color_jco()

PairedDivergence_df_Pairs%>%
  ggplot(aes(x=Group,fill=as.factor(Dist>0)))+
  geom_bar(position="fill")+
  theme_pubr()+scale_fill_jco()

library(performance)
library(fitdistrplus)
plotdist(PairedDivergence_df_Pairs$HammingDistance, histo = TRUE, demp = TRUE)
descdist(PairedDivergence_df_Pairs$HammingDistance, discrete=FALSE, boot=500)
fnb<-fitdist(PairedDivergence_df_Pairs$HammingDistance, "nbinom")
fp<-fitdist(PairedDivergence_df_Pairs$HammingDistance, "pois")
plot.legend <- c( "nbinom", "pois")
denscomp(list(fnb,fp), legendtext = plot.legend)
qqcomp(list(fnb,fp), legendtext = plot.legend)
cdfcomp(list(fnb,fp), legendtext = plot.legend)
ppcomp(list(fnb,fp), legendtext = plot.legend)

gofstat(list(fnb,fp))



library(glmmTMB)
zinbm0<-glmmTMB(HammingDistance ~ Group*TimeBetweenSamples, 
                family="nbinom2",
                ziformula=~1,data = PairedDivergence_df_Pairs)

zinbm1<-glmmTMB(HammingDistance ~ Group*TimeBetweenSamples, 
                family="nbinom2",
                data = PairedDivergence_df_Pairs)


anova(zinbm1,zinbm0)
model_performance(zinbm0)
check_model(zinbm0)
compare_performance(zinbm0,zinbm1)

summary(zinbm0)
pairs(lsmeans(zinbm0, ~ Group),adjust="fdr")


check_model(zinbm0)




model <- lm(HammingDistance ~ Group*TimeBetweenSamples, PairedDivergence_df_Pairs)

summary(model)
plot(check_distribution(model))
library(lsmeans)
pairs(lsmeans(model, ~ Group))


check_model(model)
qqnorm(resid(model))
qqline(resid(model))
hist(resid(model))
plot(fitted(model),resid(model))
abline(h=0)



PairedDivergence_df_Pairs_Divergent<-PairedDivergence_df_Pairs%>%filter(HammingDistance>0)
PairedDivergence_df_Pairs_Divergent%>%
  ggplot(aes(x=TimeBetweenSamples,y=Dist,color=Group))+
  geom_point()+geom_smooth(method = "lm",se = F)+
  theme_pubr()+scale_color_jco()

PairedDivergence_df_Pairs%>%
  ggplot(aes(x=TimeBetweenSamples,y=HammingDistance,color=Group))+
  geom_jitter()+geom_smooth(method = "lm",se = F)+
  theme_pubr()+scale_color_jco()+stat_cor()


write.csv(FastaNames_1%>%dplyr::select(study_id,mrn,RDV_Hosp,Group,participant),"TreeMetadata.csv",quote = F, row.names = F)

