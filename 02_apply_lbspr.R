# apply LBSPR to different simulated data sets
# need to calculate the 'true' SBPR for each
# scenarios: constant recruitment, single month
# scenarios: sporadic recruitment, single month
# scenarios: constant recruitment, double month
# scenarios: sporadic recruitment, double month

# Pull data from mantis_om, populate below
# compare the actual Bzero at a given time to the estimated for different scenarios
# and 

#devtools::install_github("AdrianHordyk/LBSPR")
library(LBSPR)

# example
MyPars <- new("LB_pars")
MyPars@Species <- "MySpecies"
MyPars@Linf <- 100 
MyPars@L50 <- 66 
MyPars@L95 <- 70
MyPars@MK <- 1.5 
MyPars@L_units <- "mm"
datdir <- DataDir()
Len1 <- new("LB_lengths", LB_pars=MyPars, file=paste0(datdir, "/LFreq_MultiYr.csv"), 
            dataType="freq")

#==quote Prince 2015 for average crustacean MK ratio
#====simple model 
use_it<-(dcast(filter(use_df,model=="Simple"),Month~Size,value.var="Var"))

MyPars2 <- new("LB_pars")
MyPars2@Species <- "MySpecies"
MyPars2@Linf <- 120
MyPars2@L50 <- 75 
MyPars2@L95 <- 80
MyPars2@MK <- 1.55 
MyPars2@L_units <- "mm"
Len2 <- Len1
Len2@LMids<-as.numeric(colnames(use_it)[-1])
Len2@LData<-t(use_it[,-1])
Len2@Years<-seq(1,12)
Len2@NYears<-12

plotSize(Len2)

myFit_simp <- LBSPRfit(MyPars2, Len2)
myFit_simp@Ests

#====complex model 
use_it<-(dcast(filter(use_df,model=="Complex"),Month~Size,value.var="Var"))
#==project under no fishing to get Linf for each OM
MyPars2 <- new("LB_pars")
MyPars2@Species <- "MySpecies"
MyPars2@Linf <- 165 
MyPars2@L50 <- 75 
MyPars2@L95 <- 80
MyPars2@MK <- 1.55 
MyPars2@L_units <- "mm"
Len2 <- Len1
Len2@LMids<-as.numeric(colnames(use_it)[-1])
trmp<-t(use_it[,-1])
Len2@LData<-trmp[,which(apply(trmp,2,sum)>0)]
Len2@Years<-seq(1,12)[which(apply(trmp,2,sum)>0)]
Len2@NYears<-length(Len2@Years)

plotSize(Len2)

apply(Len2@LData,2,sum)
ugh<-Len2@LData
for(x in 2:nrow(ugh))
  ugh[x,]<-ugh[x-1,]+ugh[x,]

sweep(ugh,2,apply(Len2@LData,2,sum),"/")

myFit_comp <- LBSPRfit(MyPars2, Len2)
myFit_comp@Ests



plot_spr<-data.frame(SPR=as.numeric(c(myFit_simp@Ests[,4],myFit_comp@Ests[,4])),
           OM=c(rep("Simple",12),rep("Complex",Len2@NYears)),
           ests="LBSPR")


#==calculate true SPR in each month
ber<-apply(outs_complex$eq_spawner_n,1,sum)
yer<-apply(sweep(outs_complex$n_matrix,2,params$maturity_at_size,"*"),1,sum)
yer_c<-yer[(length(yer)-11):length(yer)]
true_spr_com<-yer_c/ber

berr<-apply(outs_simple$eq_spawner_n,1,sum)
yerr<-apply(sweep(outs_simple$n_matrix,2,params$maturity_at_size,"*"),1,sum)
yerr_cc<-yerr[(length(yerr)-11):length(yerr)]
true_spr_sim<-yerr_cc/berr

add_spr<-data.frame(SPR=as.numeric(c(true_spr_sim,true_spr_com)),
                    OM=c(rep("Simple",12),rep("Complex",12)),
                    ests="True")
in_spr<-rbind(plot_spr,add_spr)

png('plots/figure_lbspr_comp.png',height=6,width=6,res=350,units='in')
ggplot(in_spr)+
  geom_boxplot(aes(x=OM,y=SPR,fill=ests),position=position_dodge(0))+
  theme_bw()+
  expand_limits(y=0)+
  theme(legend.position=c(0.2,.2),
        legend.title = element_blank())+
  scale_fill_manual(values=c("#f8766d","#619cff","#00ba38"))
dev.off()
        

add_spr2<-data.frame(SPR=as.numeric(c(em_spr)),
                    OM=c(rep("Complex",12)),
                    ests="EM")
in_spr2<-rbind(in_spr,add_spr2)

median(filter(in_spr,))

png('plots/figure_lbspr_comp_w_em.png',height=5,width=5,res=350,units='in')
ggplot(in_spr2)+
  geom_boxplot(aes(x=OM,y=SPR,fill=ests),position=position_dodge(0))+
  theme_bw()+
  expand_limits(y=0)+
  theme(legend.position=c(0.2,.2),
        legend.title = element_blank())+
  scale_fill_manual(values=c("#00ba38","#f8766d","#619cff"))
dev.off()

        