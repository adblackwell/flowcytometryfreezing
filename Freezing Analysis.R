library(beeswarm)
library(wesanderson)
library(lmerTest)
library(insight)
library(ICC)

#sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ICC_2.3.0         brms_2.14.0       Rcpp_1.0.5        rmcorr_0.4.1     
# [5] insight_0.9.6     lmerTest_3.1-2    lme4_1.1-23       Matrix_1.2-18    
# [9] wesanderson_0.3.6 beeswarm_0.2.3   
# 
# loaded via a namespace (and not attached):
# [1] nlme_3.1-149         matrixStats_0.57.0   xts_0.12.1           threejs_0.3.3      
# [5] rstan_2.19.3         numDeriv_2016.8-1.1  backports_1.1.10     tools_4.0.3        
# [9] R6_2.5.0             DT_0.16              colorspace_1.4-1     withr_2.3.0        
# [13] prettyunits_1.1.1    tidyselect_1.1.0     gridExtra_2.3        mnormt_2.0.2      
# [17] processx_3.4.4       Brobdingnag_1.2-6    emmeans_1.5.1        compiler_4.0.3    
# [21] cli_2.1.0            shinyjs_2.0.0        sandwich_3.0-0       colourpicker_1.1.0
# [25] scales_1.1.1         dygraphs_1.1.1.6     mvtnorm_1.1-1        psych_2.0.9       
# [29] ggridges_0.5.2       callr_3.5.1          StanHeaders_2.21.0-6 stringr_1.4.0     
# [33] digest_0.6.27        minqa_1.2.4          base64enc_0.1-3      pkgconfig_2.0.3   
# [37] htmltools_0.5.0      fastmap_1.0.1        htmlwidgets_1.5.2    rlang_0.4.8       
# [41] rstudioapi_0.11      shiny_1.5.0          generics_0.0.2       zoo_1.8-8         
# [45] crosstalk_1.1.0.1    gtools_3.8.2         dplyr_1.0.2          inline_0.3.16     
# [49] magrittr_1.5         loo_2.3.1            bayesplot_1.7.2      fansi_0.4.1       
# [53] munsell_0.5.0        abind_1.4-5          lifecycle_0.2.0      multcomp_1.4-14   
# [57] stringi_1.5.3        MASS_7.3-53          pkgbuild_1.1.0       plyr_1.8.6        
# [61] grid_4.0.3           parallel_4.0.3       promises_1.1.1       crayon_1.3.4      
# [65] miniUI_0.1.1.1       lattice_0.20-41      splines_4.0.3        tmvnsim_1.0-2     
# [69] ps_1.4.0             pillar_1.4.6         igraph_1.2.6         boot_1.3-25       
# [73] estimability_1.3     markdown_1.1         shinystan_2.5.0      codetools_0.2-16  
# [77] reshape2_1.4.4       stats4_4.0.3         rstantools_2.1.1     glue_1.4.2        
# [81] RcppParallel_5.0.2   vctrs_0.3.4          nloptr_1.2.2.2       httpuv_1.5.4      
# [85] gtable_0.3.0         purrr_0.3.4          assertthat_0.2.1     ggplot2_3.3.2     
# [89] mime_0.9             xtable_1.8-4         coda_0.19-4          later_1.1.0.1     
# [93] survival_3.2-7       rsconnect_0.8.16     tibble_3.0.4         shinythemes_1.1.2 
# [97] statmod_1.4.35       TH.data_1.0-10       ellipsis_0.3.1       bridgesampling_1.0-0


#Col2Alpha
col.alpha <- function (acol, alpha = 0.5){  # color function for plotting
  acol <- col2rgb(acol)
  acol.red <- acol["red",]/255
  acol.green <- acol["green",]/255
  acol.blue <- acol["blue",]/255
  acol <- mapply(function(red, green, blue, alphas) rgb(red, green, blue, alphas), acol.red, acol.green, acol.blue, alpha)
  as.character(acol)
}

#Load the data
allset<-read.csv("FlowCytometryFreezingData.csv",header=TRUE,stringsAsFactors=FALSE)

#Calculate lymphocyte percents from flow cytometry
allset$LymphFlowPct1<-allset$LymphCnt1/(allset$LymphCnt1+allset$NeutroCnt1)
allset$LymphFlowPct2<-allset$LymphCnt2/(allset$LymphCnt2+allset$NeutroCnt2)
allset$LymphFlowPct3<-allset$LymphCnt3/(allset$LymphCnt3+allset$NeutroCnt3)
allset$LymphFlowPct4<-allset$LymphCnt4/(allset$LymphCnt4+allset$NeutroCnt4)
allset$LymphFlowPct5<-allset$LymphCnt5/(allset$LymphCnt5+allset$NeutroCnt5)
allset$LymphFlowPct6<-allset$LymphCnt6/(allset$LymphCnt6+allset$NeutroCnt6)
allset$LymphFlowPct7<-allset$LymphCnt7/(allset$LymphCnt7+allset$NeutroCnt7)
allset$LymphFlowPct8<-allset$LymphCnt8/(allset$LymphCnt8+allset$NeutroCnt8)
allset$LymphCntM<-rowMeans(allset[,c("LymphCnt1","LymphCnt2","LymphCnt3","LymphCnt4","LymphCnt5","LymphCnt6","LymphCnt7","LymphCnt8")],na.rm=TRUE)
allset$LymphFlowPct<-rowMeans(allset[,c("LymphFlowPct1","LymphFlowPct2","LymphFlowPct3","LymphFlowPct4","LymphFlowPct5","LymphFlowPct6","LymphFlowPct7","LymphFlowPct8")],na.rm=TRUE)

#Lymphocyte recovery calculation
#Started with 10ul blood, ended up in 150ul FSB at the end. 
allset$LymphRecovery <- allset$LymphCntM/(allset$LymphQBCCnt*1000/15)*100

#Frozen samples used 20ul of rewashed cells not 10ul
allset$LymphRecovery[allset$Frozen==1]<-allset$LymphRecovery[allset$Frozen==1]/2


#summarize samples run and when
weeksummary<-table(allset$PID,allset$Weeks)
w2<-as.data.frame.matrix(weeksummary)
rbind(colSums(w2[w2$`0`==1,]),colSums(w2[w2$`2`==1,]),colSums(w2[w2$`4`==1,]),colSums(w2[w2$`14`==1,]))

#Viability Analysis
ViablebyTime<-reshape(allset[,c("PID","Weeks","Viable","LymphRecovery")],idvar="PID",timevar="Weeks",direction="wide")

ViablebyTime2<-ViablebyTime[!is.na(ViablebyTime$Viable.2) & !is.na(ViablebyTime$Viable.4),]
ViablebyTime2

LymphByTime2<-ViablebyTime[!is.na(ViablebyTime$LymphRecovery.2) & !is.na(ViablebyTime$LymphRecovery.4),]
LymphByTime2

#For Fresh samples, check if time to preparation affected results. 
summary(lm(LymphRecovery~waittime+Age+Sex,data=allset[allset$Frozen==0,]))
summary(lm(LymphFlowPct~waittime+Age+Sex,data=allset[allset$Frozen==0,]))
summary(lm(GransPct~waittime+Age+Sex,data=allset[allset$Frozen==0,]))
summary(lm(CD45PerM~waittime+Age+Sex,data=allset[allset$Frozen==0,]))
summary(lm(NKPerM~waittime+Age+Sex,data=allset[allset$Frozen==0,]))
summary(lm(CD19PerM~waittime+Age+Sex,data=allset[allset$Frozen==0,]))

#restructure lymphocytes and neutrophils to long format so each sample run is a row. There are 8 sets per sample
lmnset<-allset[,c("PID","Batch","Weeks","Frozen","LymphRecovery","LymphQBCPct","GransPct","Viable","LymphCnt1","LymphCnt2","LymphCnt3","LymphCnt4","LymphCnt5","LymphCnt6","LymphCnt7","LymphCnt8","NeutroCnt1","NeutroCnt2","NeutroCnt3","NeutroCnt4","NeutroCnt5","NeutroCnt6","NeutroCnt7","NeutroCnt8")]

lmnset2<-reshape(lmnset,varying=list(LymphCnt=c("LymphCnt1","LymphCnt2","LymphCnt3","LymphCnt4","LymphCnt5","LymphCnt6","LymphCnt7","LymphCnt8"),NeutroCnt=c("NeutroCnt1","NeutroCnt2","NeutroCnt3","NeutroCnt4","NeutroCnt5","NeutroCnt6","NeutroCnt7","NeutroCnt8")),times=1:8,direction="long")
lmnset2$LymphFlowPct<-lmnset2$LymphCnt1/(lmnset2$LymphCnt1+lmnset2$NeutroCnt1)
lmnset2$NeutroFlowPct<-lmnset2$NeutroCnt1/(lmnset2$LymphCnt1+lmnset2$NeutroCnt1)

#Calculate means (same as above, just a different format)
lmnset2Mean<-aggregate(cbind(LymphQBCPct,GransPct,LymphFlowPct,NeutroFlowPct,LymphRecovery)~PID+Frozen,data=lmnset2,FUN=mean)

lmnset2MeanTemp<-lmnset2Mean[lmnset2Mean$Frozen==1,c("PID","LymphFlowPct","NeutroFlowPct")]
names(lmnset2MeanTemp)[2:3]<- c("LymphFlowPctFrozen","NeutroFlowPctFrozen")
lmnset2Mean<-merge(lmnset2Mean,lmnset2MeanTemp,all=TRUE)

#check correlation in lymphocyte percent in flow and with QBC
cor.test(~ LymphQBCPct + LymphFlowPct,data=lmnset2Mean[lmnset2Mean$Frozen==0,])
cor.test(~ LymphQBCPct + LymphFlowPct,data=lmnset2Mean[lmnset2Mean$Frozen==1,])
cor.test(~ LymphFlowPct + LymphFlowPctFrozen,data=lmnset2Mean[lmnset2Mean$Frozen==0,])

t.test(lmnset2Mean[lmnset2Mean$Frozen==0,]$LymphQBCPct,lmnset2Mean[lmnset2Mean$Frozen==0,]$LymphFlowPct*100,paired=TRUE)
t.test(lmnset2Mean[lmnset2Mean$Frozen==1,]$LymphQBCPct,lmnset2Mean[lmnset2Mean$Frozen==1,]$LymphFlowPct*100,paired=TRUE)

#Aggregate by week and batch
lmnset2Mean2<-aggregate(cbind(LymphQBCPct,GransPct,LymphFlowPct,NeutroFlowPct,LymphRecovery)~PID+Frozen+Weeks+Batch,data=lmnset2,FUN=mean,na.rm=TRUE)

lmnset2Mean2$change <- (lmnset2Mean2$LymphFlowPct*100)-lmnset2Mean2$LymphQBCPct
#Check Batch effects
mod<-lm(change~as.factor(Batch), data=lmnset2Mean2[lmnset2Mean2$Frozen==1,])
anova(mod)

mod<-lm(LymphRecovery~as.factor(Batch), data=lmnset2Mean2[lmnset2Mean2$Frozen==1,])
anova(mod)

#Further examine lymphocyte recovery
allset$change <- (allset$LymphFlowPct*100)-allset$LymphQBCPct
cor.test(allset$LymphRecovery[allset$Frozen==1],allset$change[allset$Frozen==1])

mean(allset$LymphRecovery[allset$Frozen==0],na.rm=TRUE)
sd(allset$LymphRecovery[allset$Frozen==0],na.rm=TRUE)
mean(allset$LymphRecovery[allset$Frozen==1],na.rm=TRUE)
sd(allset$LymphRecovery[allset$Frozen==1],na.rm=TRUE)

mod<-lm(LymphRecovery~as.factor(Weeks), data=lmnset2Mean2[lmnset2Mean2$Frozen==1,])
anova(mod)

#Examine variance contributed by Week, PID, and Batch
mod<-lmer(LymphRecovery~WeekExact+(1|PID)+(1|Batch)+(1|FreezeBatch),data=allset[allset$Frozen==1,])
Recoveryvar<-c(get_variance(mod)$var.intercept,get_variance(mod)$var.fixed)/(get_variance(mod)$var.random+get_variance(mod)$var.residual+get_variance(mod)$var.fixed)
Recoveryvar

#Same thing for Change
mod<-lmer(change~WeekExact+(1|PID)+(1|Batch)+(1|FreezeBatch),data=allset[allset$Frozen==1,])
Changevar<-c(get_variance(mod)$var.intercept,get_variance(mod)$var.fixed)/(get_variance(mod)$var.random+get_variance(mod)$var.residual+get_variance(mod)$var.fixed)
Changevar

#Check correlation between Week 2 and Week 4 values
Week2<-lmnset2Mean2[lmnset2Mean2$Weeks==2,c("PID","LymphRecovery","LymphFlowPct","LymphQBCPct","change")]
Week4<-lmnset2Mean2[lmnset2Mean2$Weeks==4,c("PID","LymphRecovery","LymphFlowPct","change")]
names(Week4)<-c("PID","LymphRecovery.4","LymphFlowPct.4","change.4")
all24<-merge(Week2,Week4)

cor.test(all24$LymphRecovery,all24$LymphRecovery.4)
cor.test(all24$LymphFlowPct,all24$LymphFlowPct.4)
cor.test(all24$change,all24$change.4)

t.test(lmnset2Mean[lmnset2Mean$Frozen==0,]$LymphFlowPct,lmnset2Mean[lmnset2Mean$Frozen==0,]$LymphFlowPctFrozen,paired=TRUE)

t.test(lmnset2Mean[lmnset2Mean$Frozen==0,]$NeutroFlowPct,lmnset2Mean[lmnset2Mean$Frozen==0,]$NeutroFlowPctFrozen,paired=TRUE)


#Figure 1 Lymphocyte plot
library(png)
#load example scatters exported from the Guava software
frozen<-readPNG("FrozenScatter.png")
fresh<-readPNG("FreshScatter.png")

cols<-c("black",rev(wes_palette("Chevalier1",4))[c(1,4,3,2)])
tiff("Lymphocytes2.tif",width=3500,height=2400,units="px",compression="lzw",pointsize=14,res=350)
layout(matrix(c(1,1, 2,2, 3,3, 4,4, 
                5,5,5 ,6,6,6, 7,7),byrow=TRUE,ncol=8),heights=c(1,0.9))


#1
par(mar=c(8,4,4,0))
plot(c(0,10000),c(0,10000),type="n",xaxs="i",yaxs="i",xlab="Forward Scatter",ylab="Side Scatter",xaxt="n",yaxt="n")
axis(1,at=c(0,2000,4000,6000,8000,10000),labels=c(0,NA,4000,NA,8000,NA))
axis(2,at=c(0,2000,4000,6000,8000,10000),labels=c(0,NA,4000,NA,8000,NA))
rasterImage(fresh,0,0,10000,10000)
polygon(c(6800,10000,10000,6800),c(0,0,4000,4000),col="white",border=NA)
box()
text(4000,2000,"Lymphocytes",adj=0)
text(1000,9000,"Granulocytes",adj=0)
at<-par("usr")[1]- (par("usr")[2]-par("usr")[1])*0.25
mtext("A",3,line=1,at=at,cex=2)
mtext("Fresh",3,adj=0.5,line=0.5,cex=0.75)

#2
par(mar=c(8,1,4,3))
plot(c(0,10000),c(0,10000),type="n",xaxs="i",yaxs="i",xlab="Forward Scatter",ylab=NA,yaxt="n",xaxt="n")
axis(1,at=c(0,2000,4000,6000,8000,10000),labels=c(0,NA,4000,NA,8000,NA))
axis(2,labels=FALSE)
rasterImage(frozen,0,0,10000,10000)
polygon(c(6800,10000,10000,6800),c(0,0,4000,4000),col="white",border=NA)
box()
text(4000,2000,"Lymphocytes",adj=0)
text(1000,9000,"Granulocytes",adj=0)
mtext("Frozen",3,adj=0.5,line=0.5,cex=0.75)

#3
par(mar=c(5,5,2,2))
boxplot(I(LymphFlowPct*100)~Frozen,col=cols[2:3],data=lmnset2,names=c("Fresh","Frozen"),ylab="Lymphocyte % (Flow Cytometry)",xlab=NA,ylim=c(0,100))
beeswarm(I(LymphFlowPct*100)~Frozen,data=lmnset2,col="black",pch=19,spacing=0.8,cex=0.3,add=TRUE)
at<-par("usr")[1]- (par("usr")[2]-par("usr")[1])*0.37
mtext("B",3,line=-1,at=at,cex=2)

#4
par(mar=c(6,4,2,2))
plot(I(LymphFlowPct*100)~I(LymphQBCPct+runif(nrow(lmnset2),-0.25,0.25)),data=lmnset2,col=col.alpha(cols[2:3],0.75)[lmnset2$Frozen+1],ylab="Lymphocyte % (Flow Cytometry)",xlab="Lymphocyte % (QBC)",ylim=c(0,100))
abline(0,1,lty=2,lwd=2)
abline(lm(I(LymphFlowPct*100)~LymphQBCPct,data=lmnset2[lmnset2$Frozen==0,]),col=cols[2],lwd=2)
abline(lm(I(LymphFlowPct*100)~LymphQBCPct,data=lmnset2[lmnset2$Frozen==1,]),col=cols[3],lwd=2)

at<-par("usr")[1]- (par("usr")[2]-par("usr")[1])*0.35
mtext("C",3,line=-1,at=at,cex=2)

#5
par(mar=c(4,5,1,2))
boxplot(allset$LymphRecovery~allset$Weeks,col=cols[2:4],names=c("0 W","2 W","4 W","14 W"),ylab="Lymphocyte Recovery (%)",xlab=NA,log="y")
beeswarm(allset$LymphRecovery~allset$Weeks,col="black",pch=19,spacing=0.9,cex=0.6,add=TRUE)
at<-par("usr")[1]- (par("usr")[2]-par("usr")[1])*0.2
mtext("D",3,line=-1,at=at,cex=2)

#6
par(mar=c(4,5,1,2))
boxplot(allset$change~allset$Weeks,col=cols[2:4],names=c("0 W","2 W","4 W","14 W"),ylab="Change in Lymphocyte (%)",xlab=NA)
beeswarm(allset$change~allset$Weeks,col="black",pch=19,spacing=0.9,cex=0.6,add=TRUE)
at<-par("usr")[1]- (par("usr")[2]-par("usr")[1])*0.2

mtext("E",3,line=-1,at=at,cex=2)

#7
par(mar=c(4,4,1,2))
barplot(cbind(Recoveryvar[c(1,3,2,4)],Changevar[c(1,3,2,4)])*100,ylim=c(0,100),col=cols[-1],ylab="Proportion of Variance",names.arg=c("R","C"))
legend(x=0,y=100,legend=rev(c("ID","Freeze","Thaw","Storage")),fill=rev(cols[-1]),bty="n")
at<-par("usr")[1]- (par("usr")[2]-par("usr")[1])*0.35
mtext("F",3,line=-1,at=at,cex=2)
dev.off()

#Examine viability
m<-lm(Viable~as.factor(Weeks),data=allset)
t.test(Viable~Weeks,data=allset[allset$Weeks %in% c(2,4),])
t.test(Viable~Weeks,data=allset[allset$Weeks %in% c(2,14),])
t.test(Viable~Weeks,data=allset[allset$Weeks %in% c(4,14),])
mean(allset$Viable,na.rm=TRUE)
sd(allset$Viable,na.rm=TRUE)

cor.test(allset$Viable,allset$LymphRecovery)
cor.test(allset$Viable,allset$LymphFlowPct)
cor.test(allset$TotalViable,allset$LymphRecovery)
cor.test(allset$TotalViable,allset$LymphFlowPct)
cor.test(allset$LymphViable,allset$LymphRecovery)
cor.test(allset$LymphViable,allset$LymphFlowPct)
cor.test(allset$GranViable,allset$LymphRecovery)
cor.test(allset$GranViable,allset$LymphFlowPct)
cor.test(allset$TotalViable,allset$LymphViable)
plot(allset$LymphViable,allset$LymphRecovery)

#Variablility in duplicates from different sets

repeatsets<-list(
Lymph=c("LymphFlowPct1","LymphFlowPct2","LymphFlowPct3","LymphFlowPct4","LymphFlowPct5","LymphFlowPct6","LymphFlowPct7","LymphFlowPct8"),
CD45=c("CD45Per1","CD45Per2"),
CD3=c("CD3Per1","CD3Per2","CD3Per7"),
CD3CD45=c("CD3CD45Per1","CD3CD45Per2"),
CD4=c("CD4Per4","CD4Per5","CD4Per6"),
CD4.2=c("CD4Per2","CD4Per4","CD4Per5","CD4Per6"),
CD8=c("CD8Per3","CD8Per7"),
CD19=c("CD19Per1","CD19Per8"),
CD28=c("CD28Per3","CD28Per4"),
CD45RA=c("CD45RAPer3","CD45RAPer4","CD45RAPer6"),
CD57=c("CD57Per3","CD57Per4"))

##Table 2 CVs
CVs<-lapply(repeatsets,function(s) apply(allset[,s],1,function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)))
CVtab<-Reduce(rbind,lapply(CVs,function(x) t(aggregate(x,by=list(Frozen=allset$Frozen),FUN=mean,na.rm=TRUE)[,2])))
row.names(CVtab)<-names(repeatsets)  
CVtab

CVtabsd<-Reduce(rbind,lapply(CVs,function(x) t(aggregate(x,by=list(Frozen=allset$Frozen),FUN=sd,na.rm=TRUE)[,2])))
row.names(CVtabsd)<-names(repeatsets)  
CVtabsd

ICCs<-Reduce(rbind,lapply(repeatsets,function(s){
  rs<-reshape(allset[,c("PID","Frozen","Weeks",s)],varying=list(var=s),idvar=c("PID","Frozen","Weeks"),direction="long",v.names="var")
  IC1<-ICCest(PID,var,rs[rs$Frozen==0,])
  IC2<-ICCest(PID,var,rs[rs$Frozen==1,])
  return(c(IC1$ICC,IC2$ICC))
} ))
row.names(ICCs)<-names(repeatsets)  
ICCs

plot(CD19Per1~CD19Per8,data=allset,col=c("black","red")[Frozen+1])
plot(CD4Per2~CD4Per6,data=allset,col=c("black","red")[Frozen+1])
plot(CD4Per2~CD4Per4,data=allset,col=c("black","red")[Frozen+1])
plot(CD4Per2~CD4Per5,data=allset,col=c("black","red")[Frozen+1])
plot(NKPer1~NKPer7,data=allset,col=c("black","red")[Frozen+1])
abline(0,1)
#on frozen NKPer1 looks elevated. Might be defined differently though
#NKPer1 is CD16/CD56pCD3nCD19nCD45p / CD45p
#NKPer7 is CD56 or CD16 CD3n / Lymph

allset$CD4.CPPerM<-rowMeans(allset[,c("CD4Per4","CD4Per5","CD4Per6")],na.rm=TRUE)
allset$CD8.CPPerM<-rowMeans(allset[,c("CD8Per3","CD8Per7")],na.rm=TRUE)

originals<-names(allset)[grep("Per[1-8M]?$",names(allset))]
#In this case just use the fresh value twice if someone was run twice.
Fresh<-allset[allset$Frozen==0,c("PID","LymphRecovery","Weeks","Batch","FreezeBatch","WeekExact",originals)]
Frozen<-allset[allset$Frozen==1,c("PID","LymphRecovery","Weeks","Batch","FreezeBatch","WeekExact",originals)]
names(Frozen)[-1]<-paste0(names(Frozen)[-1],".F")
all<-merge(Fresh,Frozen)
#check sample sizes
data.frame(names(Fresh),colSums(!is.na(Fresh)),colSums(!is.na(Frozen)))

#Function to use to summarize stats.
getstats<-function(variables,data){
  t(sapply(variables,function(v){
    if(sum(!is.na(data[,v]))>2){
      cc<-cor.test(data[,v],data[,paste0(v,".F")],alternative="greater")
      tt<-t.test(data[,paste0(v,".F")],data[,v],paired=TRUE)
      out<-c(r=cc$estimate,pr=cc$p.value,M1=mean(data[,v],na.rm=TRUE),M2=mean(data[,paste0(v,".F")],na.rm=TRUE),t=tt$estimate,pt=tt$p.value,c1=tt$conf.int[1],c2=tt$conf.int[2])
    }
    else c(NA,NA,NA,NA)
  }))
}



#various statistical tables, some of which are reported in the paper. See code further below where these get formatted into the Tables for the paper.
ffcor<-getstats(originals,all)
ffcor

allM<-aggregate(all[,-1],by=list(PID=all$PID),mean,na.rm=TRUE)
ffcorM<-getstats(originals,allM)
ffcorM

all2W<-all[all$Weeks.F==2,]
ffcor.2W<-getstats(originals,all2W)
ffcor.2W
  
all4W<-all[all$Weeks.F==4,]
ffcor.4W<-getstats(originals,all4W)
ffcor.4W

all14W<-all[all$Weeks.F==14,]
ffcor.14W<-getstats(originals,all14W)
all14W

#Some are percents of some other subset, as follows:
LPercents<-c("CD127Per6","CD16CD56Per1","CD16Per7","CD19Per1","CD19Per8","CD19PerM","CD20Per","CD25Per6","CD27Per","CD28Per3","CD28Per4","CD28PerM","CD3Per1","CD3Per2","CD3Per7","CD3PerM","CD45Per1","CD45Per2","CD45PerM","CD45RAPer3","CD45RAPer4","CD45RAPer6", "CD45RAPerM","CD4Per2","CD4Per4","CD4Per5","CD4Per6","CD4PerM","CD4.CPPerM","CD56Per7","CD57Per3","CD57Per4","CD57PerM","CD8Per2","CD8Per3","CD8Per7","CD8.CPPerM","CD8PerM","IgDPer","CD183Per5","CD196Per5","CD294Per5")

#Now adjust for LPercents for CD45 component from set 1 and 2 and then recheck
Fresh.Ladj<-cbind(Fresh[,c("PID","Weeks")],Fresh[,LPercents]/Fresh$CD45PerM*100)
Frozen.Ladj<-cbind(Frozen[,c("PID","Weeks.F")],Frozen[,paste0(LPercents,".F")]/Frozen$CD45PerM*100)
all.Ladj<-merge(Fresh.Ladj,Frozen.Ladj)

ffcoradj<-getstats(LPercents,all.Ladj)
ffcoradj

all.LadjM<-aggregate(all.Ladj[,-1],by=list(PID=all.Ladj$PID),mean,na.rm=TRUE)
ffcoradjM<-getstats(LPercents,all.LadjM)
ffcoradjM

#just 2 week
all2Wadj<-all.Ladj[all.Ladj$Weeks.F==2,]
ffcor.2Wadj<-getstats(LPercents,all2Wadj)
ffcor.2Wadj
  
all4Wadj<-all.Ladj[all.Ladj$Weeks.F==4,]
ffcor.4Wadj<-getstats(LPercents,all4Wadj)
ffcor.4Wadj

all14Wadj<-all.Ladj[all.Ladj$Weeks.F==14,]
ffcor.14Wadj<-getstats(LPercents,all14Wadj)
ffcor.14Wadj

#Remove low lymphrecovery and check again
#LymphCntM best I've got
allset2<-allset[allset$LymphRecovery>10,]
Fresh2<-allset2[allset2$Frozen==0,c("PID","Weeks",originals)]
Frozen2<-allset2[allset2$Frozen==1,c("PID","Weeks",originals)]
names(Frozen2)[-1]<-paste0(names(Frozen2)[-1],".F")
all2<-merge(Fresh2,Frozen2)

ffcor2<-getstats(originals,all2)
ffcor2

allM2<-aggregate(all2[,-1],by=list(PID=all2$PID),mean,na.rm=TRUE)
ffcorM2<-getstats(originals,allM2)
ffcorM2

#adjusted and low conc eliminated
Fresh2.Ladj<-cbind(Fresh2[,c("PID","Weeks")],Fresh2[,LPercents]/Fresh2$CD45PerM*100)
Frozen2.Ladj<-cbind(Frozen2[,c("PID","Weeks.F")],Frozen2[,paste0(LPercents,".F")]/Frozen2$CD45PerM*100)
all2.Ladj<-merge(Fresh2.Ladj,Frozen2.Ladj)
ffcoradj2<-getstats(LPercents,all2.Ladj)
ffcoradj2  

#counts
LymphPercents<-c("CD3PerM","CD4Per2","CD4.CPPerM","CD8Per2","CD8.CPPerM","CD16Per7","CD19Per1","CD19Per8","CD20Per","CD25Per6","CD27Per","CD28PerM","CD45PerM","CD45RAPerM","CD56Per7","CD57PerM","CD127Per6","CD183Per5","CD196Per5","CD294Per5","IgDPer")
CD45Percents<-c("CD19CD45Per1","CD3CD45PerM","CD4CD3CD45Per2","CD8CD3CD45Per2","NKPer1")
CD4Percents<-c("CD4NaivePer","CD4CD45RAPer4","CD4SenPer","TH17Per","TH1Per","TH2Per","TregPer")
CD8Percents<-c("CD8NaivePer","CD8CD45RAPer3","CD8SenPer")
CD19Percents<-c("BMemoryPer","BNaivePer","BNonSwitcherPer","BPlasmaPer")

#Converting to counts, but note, not adjusting variable names to reflect using counts
allsub<-cbind(allset[,c("PID","Weeks","LymphQBCCnt","Frozen",originals)])
allsub[,LymphPercents] <- allsub[,LymphPercents]*allsub$LymphQBCCnt/100*1000
allsub[,CD45Percents] <- allsub[,CD45Percents]*allsub$CD45PerM/100
allsub[,CD4Percents] <- allsub[,CD4Percents]*allsub$CD4CD3CD45Per2/100
allsub[,CD8Percents] <- allsub[,CD8Percents]*allsub$CD8CD3CD45Per2/100
allsub[,CD19Percents] <- allsub[,CD19Percents]*allsub$CD19CD45Per1/100

Freshcnt<-allsub[allsub$Frozen==0,]
Frozencnt<-allsub[allsub$Frozen==1,]
names(Frozencnt)[-1]<-paste0(names(Frozencnt)[-1],".F")
allcnt<-merge(Freshcnt,Frozencnt)

ffcorcnt<-getstats(originals,allcnt)
ffcorcnt

#Make Tables
##Formatting function for Tables
formcor<-function(r){
  paste0(format(round(r[,1],2),nsmall=2),ifelse(r[,2]>0.05,"ns",""))
}


#Table 3. Mean percent lymphocytes positive for antibody tags in fresh and frozen samples
cbind(ffcorM[match(LymphPercents, rownames(ffcorM)),c(3,4,7,8)],
      ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),c(3,4,7,8)])

ffcorM2[match(LymphPercents, rownames(ffcorM2)),c(3,4,6)]

#Figure 2 Fresh Frozen comparison
cols<-c("black",rev(wes_palette("Chevalier1",4))[c(1,4,3,2)])
tiff("FreshFrozen.tif",width=4000,height=2000,units="px",compression="lzw",pointsize=14,res=350)
par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
plot(ffcorM[match(LymphPercents, rownames(ffcorM)),3],ffcorM[match(LymphPercents, rownames(ffcorM)),4],xlab="Fresh %",ylab="Frozen %",xlim=c(4,95),ylim=c(4,95))
polygon(c(1,1,100,100),c(0.8,1.2,120,80),border=NA,col=cols[4])
points(ffcorM[match(LymphPercents, rownames(ffcorM)),3],ffcorM[match(LymphPercents, rownames(ffcorM)),4],pch=19)
abline(0,1,lwd=2)
apply(ffcorM[match(LymphPercents, rownames(ffcorM)),c(3,4)],1,function(x) lines(c(x[1],x[1]),c(x[1],x[2])))
box()
LymphLabels<-c("CD3","CD4-PE","CD4-PerCP-Cy5.5","CD8-APC","CD8-PerCP-Cy5.5","CD16","CD19-APC","CD19-PerCP-Cy5.5","CD20","CD25","CD27","CD28","CD45","CD45RA","CD56","CD57","CD127","CD183","CD196","CD294","IgD"
)
delta<-(ffcorM[match(LymphPercents, rownames(ffcorM)),4]-ffcorM[match(LymphPercents, rownames(ffcorM)),3])/ffcorM[match(LymphPercents, rownames(ffcorM)),3]
showlabel<-abs(delta)>0.20
text(ffcorM[match(LymphPercents, rownames(ffcorM)),3][showlabel]+c(0,0,0,-0,0,0,-10,0,-5),ffcorM[match(LymphPercents, rownames(ffcorM)),4][showlabel]+c(5,-2,0,0,0,0,0,0,5),LymphLabels[showlabel],pos=4)
points(ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),3],ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),4])
apply(cbind(ffcorM[match(LymphPercents, rownames(ffcorM)),c(3,4)],ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),c(3,4)]),1,function(x) lines(c(x[1],x[3]),c(x[2],x[4]),lty=2))

plot(ffcorM[match(LymphPercents, rownames(ffcorM)),3],ffcorM[match(LymphPercents, rownames(ffcorM)),4],xlab="Fresh %",ylab="Frozen %",xlim=c(2,95),ylim=c(2,95),log="xy")
polygon(c(1,1,110,110),c(0.8,1.2,132,88),border=NA,col=cols[4])
points(ffcorM[match(LymphPercents, rownames(ffcorM)),3],ffcorM[match(LymphPercents, rownames(ffcorM)),4],pch=19)
abline(0,1,lwd=2)
apply(ffcorM[match(LymphPercents, rownames(ffcorM)),c(3,4)],1,function(x) lines(c(x[1],x[1]),c(x[1],x[2])))
box()
text(ffcorM[match(LymphPercents, rownames(ffcorM)),3][showlabel]+c(-7,0,0,-12,0,-18,-12,0,0),ffcorM[match(LymphPercents, rownames(ffcorM)),4][showlabel]+c(5,-1,0,-5,0,-4,0,0,0),LymphLabels[showlabel],pos=4)
points(ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),3],ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),4])
apply(cbind(ffcorM[match(LymphPercents, rownames(ffcorM)),c(3,4)],ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),c(3,4)]),1,function(x) lines(c(x[1],x[3]),c(x[2],x[4]),lty=2))
dev.off()

#Correlation tables
#Table 4 Pearson’s correlations between percentages measured in fresh and frozen samples part 1
cbind(LymphPercents,
      formcor(ffcorM[match(LymphPercents, rownames(ffcorM)),c(1,2)]),
      formcor(ffcor.2W[match(LymphPercents, rownames(ffcor.2W)),c(1,2)]),
      formcor(ffcor.4W[match(LymphPercents, rownames(ffcor.4W)),c(1,2)]),
      formcor(ffcor.14W[match(LymphPercents, rownames(ffcor.14W)),c(1,2)])
      )
# #Mean correlations by week
c(mean(ffcorM[match(LymphPercents, rownames(ffcorM)),1]),mean(ffcor.2W[match(LymphPercents, rownames(ffcor.2W)),1]),mean(ffcor.4W[match(LymphPercents, rownames(ffcor.4W)),1]),mean(ffcor.14W[match(LymphPercents, rownames(ffcor.14W)),1]))

#Table 4 Part 2
cbind(LymphPercents,
      formcor(ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),c(1,2)]),
      formcor(ffcor.2Wadj[match(LymphPercents, rownames(ffcor.2Wadj)),c(1,2)]),
      formcor(ffcor.4Wadj[match(LymphPercents, rownames(ffcor.4Wadj)),c(1,2)]),
      formcor(ffcor.14Wadj[match(LymphPercents, rownames(ffcor.14Wadj)),c(1,2)]))

c(mean(ffcoradjM[match(LymphPercents, rownames(ffcoradjM)),1],na.rm=TRUE),mean(ffcor.2Wadj[match(LymphPercents, rownames(ffcor.2Wadj)),1],na.rm=TRUE),mean(ffcor.4Wadj[match(LymphPercents, rownames(ffcor.4Wadj)),1],na.rm=TRUE),mean(ffcor.14Wadj[match(LymphPercents, rownames(ffcor.14Wadj)),1],na.rm=TRUE))

CD45Percents<-c("CD19CD45Per1","CD3CD45PerM","CD4CD3CD45Per2","CD8CD3CD45Per2","NKPer1")
CD4Percents<-c("CD4NaivePer","CD4CD45RAPer4","CD4SenPer","TH17Per","TH1Per","TH2Per","TregPer")
CD8Percents<-c("CD8NaivePer","CD8CD45RAPer3","CD8SenPer")
CD19Percents<-c("BMemoryPer","BNaivePer","BNonSwitcherPer","BPlasmaPer")
#Table 5. Mean percent lymphocyte subsets in fresh and frozen samples
rbind(ffcorM[match(CD45Percents, rownames(ffcorM)),c(3,4,5,7,8)],
      ffcorM[match(CD4Percents, rownames(ffcorM)),c(3,4,5,7,8)],
      ffcorM[match(CD8Percents, rownames(ffcorM)),c(3,4,5,7,8)],
      ffcorM[match(CD19Percents, rownames(ffcorM)),c(3,4,5,7,8)])

#table using counts instead of percents, not in paper      
rbind(ffcorcnt[match(CD45Percents, rownames(ffcorcnt)),c(1,3,4,6,2)],
      ffcorcnt[match(CD4Percents, rownames(ffcorcnt)),c(1,3,4,6,2)],
      ffcorcnt[match(CD8Percents, rownames(ffcorcnt)),c(1,3,4,6,2)],
      ffcorcnt[match(CD19Percents, rownames(ffcorcnt)),c(1,3,4,6,2)])


#Table 6. Pearson’s correlations between percentages obtained for fresh and frozen samples      
rbind(
cbind(ffcorM[match(CD45Percents, rownames(ffcorM)),1],
      ffcor.2W[match(CD45Percents, rownames(ffcor.2W)),1],
      ffcor.4W[match(CD45Percents, rownames(ffcor.4W)),1],
      ffcor.14W[match(CD45Percents, rownames(ffcor.14W)),1]),

cbind(ffcorM[match(CD4Percents, rownames(ffcorM)),1],
      ffcor.2W[match(CD4Percents, rownames(ffcor.2W)),1],
      ffcor.4W[match(CD4Percents, rownames(ffcor.4W)),1],
      ffcor.14W[match(CD4Percents, rownames(ffcor.14W)),1]),

cbind(ffcorM[match(CD8Percents, rownames(ffcorM)),1],
      ffcor.2W[match(CD8Percents, rownames(ffcor.2W)),1],
      ffcor.4W[match(CD8Percents, rownames(ffcor.4W)),1],
      ffcor.14W[match(CD8Percents, rownames(ffcor.14W)),1]),

cbind(ffcorM[match(CD19Percents, rownames(ffcorM)),1],
      ffcor.2W[match(CD19Percents, rownames(ffcor.2W)),1],
      ffcor.4W[match(CD19Percents, rownames(ffcor.4W)),1],
      ffcor.14W[match(CD19Percents, rownames(ffcor.14W)),1]))

#Table of p-values for correlations by week, used for Table 6
rbind(
  cbind(ffcorM[match(CD45Percents, rownames(ffcorM)),2],
        ffcor.2W[match(CD45Percents, rownames(ffcor.2W)),2],
        ffcor.4W[match(CD45Percents, rownames(ffcor.4W)),2],
        ffcor.14W[match(CD45Percents, rownames(ffcor.14W)),2]),
  
  cbind(ffcorM[match(CD4Percents, rownames(ffcorM)),2],
        ffcor.2W[match(CD4Percents, rownames(ffcor.2W)),2],
        ffcor.4W[match(CD4Percents, rownames(ffcor.4W)),2],
        ffcor.14W[match(CD4Percents, rownames(ffcor.14W)),2]),
  
  cbind(ffcorM[match(CD8Percents, rownames(ffcorM)),2],
        ffcor.2W[match(CD8Percents, rownames(ffcor.2W)),2],
        ffcor.4W[match(CD8Percents, rownames(ffcor.4W)),2],
        ffcor.14W[match(CD8Percents, rownames(ffcor.14W)),2]),
  
  cbind(ffcorM[match(CD19Percents, rownames(ffcorM)),2],
        ffcor.2W[match(CD19Percents, rownames(ffcor.2W)),2],
        ffcor.4W[match(CD19Percents, rownames(ffcor.4W)),2],
        ffcor.14W[match(CD19Percents, rownames(ffcor.14W)),2]))


#batch effects for table 6
BatchTests<-c("CD19CD45Per1","CD3CD45PerM","CD4CD3CD45Per2","CD8CD3CD45Per2","NKPer1","CD4NaivePer","CD4CD45RAPer4","CD4SenPer","TH17Per","TH1Per","TH2Per","TregPer","CD8NaivePer","CD8CD45RAPer3","CD8SenPer","BMemoryPer","BNaivePer","BNonSwitcherPer","BPlasmaPer")
library(MCMCglmm)
batch<-t(sapply(BatchTests, function(v){
  print(v)
  mod<-lmer(formula(paste0("I(",v,".F-",v,")~(1|Batch.F)+(1|FreezeBatch.F)")),data=all)
  Changevar<-c(get_variance(mod)$var.intercept)/(sum(get_variance(mod)$var.intercept)+get_variance(mod)$var.residual)
  Changevar
}))


corem<-c("CD19CD45Per1","CD3CD45PerM","NKPer1","CD4CD3CD45Per2","CD4NaivePer","CD4SenPer","CD8CD3CD45Per2","CD8NaivePer","CD8SenPer")


#Figure 3 Scatters with Percents
rndRec<-1 + (all$LymphRecovery.F>5) + (all$LymphRecovery.F>7) + (all$LymphRecovery.F>10) + (all$LymphRecovery.F>15)
lymphcols<-rev(wes_palette("Zissou1",5))
tiff("Fig 3 Scatters.tif",width=2000,height=2000,units="px",compression="lzw",pointsize=14,res=300)
layout(matrix(c(10,1,2,3, 10,4,5,6, 10,7,8,9, 12,11,11,11),byrow=T,ncol=4),widths=c(0.05,0.3,0.3,0.3),heights=c(0.3,0.3,0.3,0.05))
par(mar=c(2,2,1,1))
titles<-c("B cells","T cells","Natural Killer","CD4 T","Naive CD4","Senescent CD4 T","CD8 T","Naive CD8 T","Senescent CD8 T")
for(i in 1:length(corem)){
  minmax<-range(c(all[,corem[i]],all[,paste0(corem[i],".F")]),na.rm=TRUE)
  plot(all[,corem[i]],all[,paste0(corem[i],".F")],xlab=NA,ylab=NA,pch=19,col=lymphcols[rndRec],ylim=minmax,xlim=minmax)
  rng <- par("usr")
  text(rng[1]+1*(rng[2]-rng[1]),rng[3]+0.05*(rng[4]-rng[3]),titles[i],pos=2,col="red")
  text(rng[1]+0.05*(rng[2]-rng[1]),rng[3]+0.90*(rng[4]-rng[3]),paste0("r=",format(round(ffcorM[corem[i],1],2),nsmall=2)),pos=4,col="red")
  abline(0,1)
  abline(lm(all[,paste0(corem[i],".F")]~all[,corem[i]]),lty=2)
}
par(mar=c(0,0,0,0))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,"Frozen (%)",srt=90,cex=1.7)
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,"Fresh (%)",srt=0,cex=1.7)
dev.off()

#Figure 4 Scatters with Counts
tiff("Figure 5 ScattersCNT.tif",width=2000,height=2000,units="px",compression="lzw",pointsize=14,res=300)
layout(matrix(c(10,1,2,3, 10,4,5,6, 10,7,8,9, 12,11,11,11),byrow=T,ncol=4),widths=c(0.05,0.3,0.3,0.3),heights=c(0.3,0.3,0.3,0.05))
par(mar=c(2,2,1,1))
titles<-c("B cells","T cells","Natural Killer","CD4 T","Naive CD4","Senescent CD4 T","CD8 T","Naive CD8 T","Senescent CD8 T")
for(i in 1:length(corem)){
  minmax<-range(c(allcnt[,corem[i]],allcnt[,paste0(corem[i],".F")]),na.rm=TRUE)
  plot(allcnt[,corem[i]],allcnt[,paste0(corem[i],".F")],xlab=NA,ylab=NA,pch=19,col=lymphcols[rndRec],ylim=minmax,xlim=minmax)
  rng <- par("usr")
  text(rng[1]+1*(rng[2]-rng[1]),rng[3]+0.05*(rng[4]-rng[3]),titles[i],pos=2,col="red")
  text(rng[1]+0.05*(rng[2]-rng[1]),rng[3]+0.90*(rng[4]-rng[3]),paste0("r=",format(round(ffcorcnt[corem[i],1],2),nsmall=2)),pos=4,col="red")
  abline(0,1)
  abline(lm(allcnt[,paste0(corem[i],".F")]~allcnt[,corem[i]]),lty=2)
}
par(mar=c(0,0,0,0))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,expression(paste("Frozen (cells/",mu,"l)")),srt=90,cex=1.7)
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,expression(paste("Fresh (cells/",mu,"l)")),srt=0,cex=1.7)
dev.off()

#Figure 7 mean values by week
out<-rbind(
  cbind(ffcor.2W[match(corem, rownames(ffcor.2W)),c(5,7,8)],
        ffcor.4W[match(corem, rownames(ffcor.4W)),c(5,7,8)],
        ffcor.14W[match(corem, rownames(ffcor.14W)),c(5,7,8)]))
tiff("Figure 6 WeekChange.tif",width=2000,height=2000,units="px",compression="lzw",pointsize=14,res=300)
titles<-c("B cells","T cells","Natural Killer","CD4 T","Naive CD4","Senescent CD4 T","CD8 T","Naive CD8 T","Senescent CD8 T")
layout(matrix(c(10,1,2,3, 10,4,5,6, 10,7,8,9, 12,11,11,11),byrow=T,ncol=4),widths=c(0.05,0.3,0.3,0.3),heights=c(0.3,0.3,0.3,0.05))
par(mar=c(2,2,1,1))
cols<-wes_palette("Zissou1",5)[c(4,1,5,1,3,3,1,3,3)]
cols<-c(wes_palette("GrandBudapest1"),wes_palette("GrandBudapest2"))
cols<-c(cols[-3],cols[1:2])
for(i in 1:nrow(out)){
  plot(0,0,type="n",xlim=c(2,4),ylim=c(-15,15),xaxt="n",xlab=NA,ylab=NA,bg="black")
  axis(1,2:4,labels=c("2","4","14"))
  polygon(c(2:4,4:2),c(out[i,c(2,5,8,9,6,3)]),border=NA,col=cols[i])
  lines(2:4,c(out[i,c(1,4,7)]),col="black")
  abline(0,0)
  text(3,15,titles[i], pos=1,cex=1.5)
}  
par(mar=c(0,0,0,0))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,"Change in %",srt=90,cex=1.7)
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,"Weeks Frozen",srt=0,cex=1.7)
dev.off()

#Bland-Altman
tiff("Figure 4 Bland-Altman.tif",width=2000,height=2000,units="px",compression="lzw",pointsize=14,res=300)
layout(matrix(c(10,1,2,3, 10,4,5,6, 10,7,8,9, 12,11,11,11),byrow=T,ncol=4),widths=c(0.05,0.3,0.3,0.3),heights=c(0.3,0.3,0.3,0.05))
par(mar=c(2,2,1,1))
titles<-c("B cells","T cells","Natural Killer","CD4 T","Naive CD4","Senescent CD4 T","CD8 T","Naive CD8 T","Senescent CD8 T")
for(i in 1:length(corem)){
  avg<- (all[,paste0(corem[i],".F")] + all[,corem[i]])/2
  diff<- all[,paste0(corem[i],".F")] - all[,corem[i]]
  minmax<-max(abs(diff),na.rm=TRUE)*1.2
  plot(avg,diff,xlab=NA,ylab=NA,pch=19,col=lymphcols[rndRec],ylim=c(-45,45))

  abline(h=mean(diff,na.rm=TRUE),lty=2)
  sddiff<-sd(diff,na.rm=TRUE)
  abline(h=mean(diff,na.rm=TRUE)+2*sddiff,lty=3)
  abline(h=mean(diff,na.rm=TRUE)-2*sddiff,lty=3)
  rng <- par("usr")
  text(rng[1]+1*(rng[2]-rng[1]),rng[3]+0.05*(rng[4]-rng[3]),titles[i],pos=2,col="red")
}
par(mar=c(0,0,0,0))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,"Frozen (%) - Fresh (%)",srt=90,cex=1.7)
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,"(Frozen (%) + Fresh (%))/2",srt=0,cex=1.7)
dev.off()

#Count. Not included in paper because difficult to interpret, because relative count of different groups is so different.
tiff("Figure X Bland-Altman Count.tif",width=2000,height=2000,units="px",compression="lzw",pointsize=14,res=300)
layout(matrix(c(10,1,2,3, 10,4,5,6, 10,7,8,9, 12,11,11,11),byrow=T,ncol=4),widths=c(0.05,0.3,0.3,0.3),heights=c(0.3,0.3,0.3,0.05))
par(mar=c(2,2,1,1))
titles<-c("B cells","T cells","Natural Killer","CD4 T","Naive CD4","Senescent CD4 T","CD8 T","Naive CD8 T","Senescent CD8 T")
for(i in 1:length(corem)){
  avg<- (allcnt[,paste0(corem[i],".F")] + allcnt[,corem[i]])/2/1000
  diff<- (allcnt[,paste0(corem[i],".F")] - allcnt[,corem[i]])/1000
  minmax<-max(abs(diff),na.rm=TRUE)*1.2
  plot(avg,diff,xlab=NA,ylab=NA,pch=19,col=lymphcols[rndRec],ylim=c(-1.5,1.5))
  
  abline(h=mean(diff,na.rm=TRUE),lty=2)
  sddiff<-sd(diff,na.rm=TRUE)
  abline(h=mean(diff,na.rm=TRUE)+2*sddiff,lty=3)
  abline(h=mean(diff,na.rm=TRUE)-2*sddiff,lty=3)
  rng <- par("usr")
  text(rng[1]+1*(rng[2]-rng[1]),rng[3]+0.05*(rng[4]-rng[3]),titles[i],pos=2,col="red")
}
par(mar=c(0,0,0,0))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,expression(paste("Frozen - Fresh (x1000 cells/",mu,"l)")),srt=90,cex=1.7)
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(0.5,0.5,expression(paste("(Frozen + Fresh (x1000 cells/",mu,"l))/2")),srt=0,cex=1.7)
dev.off()