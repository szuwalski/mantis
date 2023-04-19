#===============================
# plot output from mantis model
#===============================
library(PBSmodelling)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggridges)

#  read in EM output
rep_files<-c("em/complex/mantis.rep")
#rep_files<-c("admb_devs_dat/mantis.rep")
outs<-list(list())

for(x in 1:length(rep_files))
  outs[[x]]<-readList(rep_files[x])

png('plots/figure_mod_fits.png',height=8,width=5,res=350,units='in')

#=======estimated processes vs. true================
par(mfrow=c(4,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
#================fits to data=================
div_n<-100000000
plot(outs[[1]]$survey_obs/div_n~outs[[1]]$survey_months,
     ylim=c(0,max(outs[[1]]$survey_ob/div_n * exp(1.96*sqrt(log(1+params$surv_sd^2))))),pch=16,
     las=1,ylab="Survey numbers (hundred millions)",
     xaxt='n')
for(j in 1:length(outs[[1]]$survey_months))
{
  segments(x0=outs[[1]]$survey_months[j],x1=outs[[1]]$survey_months[j],
           y0=outs[[1]]$survey_ob[j]/div_n  /   exp(1.96*sqrt(log(1+params$surv_sd^2))),
           y1=outs[[1]]$survey_ob[j]/div_n * exp(1.96*sqrt(log(1+params$surv_sd^2))))
}
lines(outs[[1]]$survey_pred_comp/div_n~outs[[1]]$survey_months,lwd=3,col='blue')
legend("topleft",bty='n',legend="Survey")
legend("top",bty='n',lty=c(NA,1),pch=c(16,NA),legend=c("Data","Model fit"),lwd=c(NA,2))
plot(outs[[1]]$cat_obs/div_n~outs[[1]]$catch_months,
     ylim=c(0,max(outs[[1]]$cat_obs/div_n,outs[[1]]$catch_pred_comp/div_n)),pch=16,las=1,
     ylab="Catch (hundred millions)",
     xlab="Month",xaxt='n')
for(j in 1:length(outs[[1]]$cat_obs))
{
  segments(x0=outs[[1]]$catch_months[j],x1=outs[[1]]$catch_months[j],
           y0=outs[[1]]$cat_obs[j]/div_n  /   exp(1.96*sqrt(log(1+params$cat_sd^2))),
           y1=outs[[1]]$cat_obs[j]/div_n * exp(1.96*sqrt(log(1+params$cat_sd^2))))
}
lines(outs[[1]]$catch_pred_comp/div_n~outs[[1]]$catch_months,lwd=2,col='seagreen')
legend("topleft",bty='n',legend="Catch")
mtext(side=2,line=2.5,"Numbers (10000000)",outer=T,adj=0.75)
mtext(side=1,line=2.5,"Month",outer=T)

plot_rec_est<-rep(0,length(params$all_months))
plot_rec_est[params$rec_months]<-log(outs[[1]]$recruits)
plot(plot_rec_est~params$all_months,type='l',ylim=c(0,max(plot_rec_est)),las=1,
     ylab='',xaxt='n',lwd=2)

plot_rec_true<-(outs_complex$recruitment)
lines(plot_rec_true,col=2,lty=2)

mtext(side=2,outer=F,line=2.75,"Recruits (ln(100,000,000))")

plot(outs[[1]]$applied_f_matrix[,8],type='l',ylim=c(0,max(outs[[1]]$applied_f_matrix)),
     las=1,ylab="Fishing mortality",xlab="Month",xlim=c(0,length(params$all_months)),lwd=2)
lines(f_mort1[,1],col=2,lty=2)
mtext(side=2,outer=F,line=2.75,"Fishing mortality")
legend("topleft",bty='n',col=c(1,2),lty=c(1,2),legend=c("Estimated","True"))
dev.off()


#==survey size comps
par_in<-ceiling(sqrt(nrow(outs[[1]]$survey_size_comp_obs)))
par(mfcol=c(par_in,par_in),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(outs[[1]]$survey_size_comp_pred_comp))
{
  plot(outs[[1]]$survey_size_comp_pred_comp[x,],yaxt='n',xaxt='n',ylim=c(0,0.5),pch=16,col='grey')
  lines(outs[[1]]$survey_size_comp_obs[x,],col=2,lty=1)
}
plot.new()
legend('center',col=c('grey','red'),pch=c(16,NA),lty=c(NA,1),legend=c("Obs","Est"),bty='n')

png('plots/figure_mod_fits_size.png',height=6,width=6,res=350,units='in')
par(mfrow=c(2,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
boxplot(outs[[1]]$survey_size_comp_obs,ylim=c(0,.3),xaxt='n',las=1)
lines(apply(outs[[1]]$survey_size_comp_pred_comp,2,median),col=2,lwd=2)
legend("topleft",bty='n',legend="Survey")
boxplot(outs[[1]]$catch_size_comp_obs,names=params$mid_pts,ylim=c(0,.3),las=1)
lines(apply(outs[[1]]$catch_size_comp_pred_comp,2,median),col=2,lwd=2)
legend("topleft",bty='n',legend="Fishery")

mtext(side=2,outer=T,line=2.5,"Proportion")
mtext(side=1,outer=T,line=2.5,"Size (mm)")
dev.off()

#==catch size comps
par_in<-ceiling(sqrt(nrow(outs[[1]]$catch_size_comp_pred_comp)))
par(mfrow=c(par_in,par_in),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(outs[[1]]$catch_size_comp_pred_comp))
{
  plot(outs[[1]]$catch_size_comp_pred_comp[x,],yaxt='n',xaxt='n',ylim=c(0,0.6),pch=16,col='grey')
  lines(outs[[1]]$catch_size_comp_obs[x,],col=2,lty=3)
}


#==size transition matrix
outs_complex$size_trans

in_mat<-outs[[1]]$size_trans_matrix
colnames(in_mat)<-outs[[1]]$sizes
rownames(in_mat)<-outs[[1]]$sizes
in_g<-data.frame(melt(t(in_mat)))
colnames(in_g)<-c("Postmolt","Premolt","Density")

p <- ggplot(in_g)
p <- p + geom_density_ridges(aes(x=Postmolt, y=Premolt, height = Density, group=Premolt,
                                 fill=stat(y),alpha=.9999),stat = "identity",scale=3) +
  scale_fill_viridis_c()+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90)) 
print(p)

##==plot true vs. estimate growth and selectivities#==growth
png('plots/figure_sel_grow.png',height=6,width=4,res=350,units='in')

par(mfcol=c(2,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
plot(outs[[1]]$survey_selectivity~params$mid_pts,type='l',las=1,xaxt='n',lwd=2,ylim=c(0,1))
lines(params$surv_sel~params$mid_pts,lty=2,col=2,lwd=2)
text(x=65,y=0.75,"Survey")
text(x=85,y=0.55,"Fishery")
lines(outs[[1]]$fish_selectivity~params$mid_pts,type='l',las=1,xaxt='n',lwd=2)
lines(params$fish_sel~params$mid_pts,lty=2,col=2,lwd=2)
mtext(side=2,line=2.5,"Selectivity")
legend('right',lty=c(1,2),col=c(1,2),legend=c("Estimated","True"),bty='n')

true_growth_inc <- params$alpha + params$beta*params$mid_pts
plot(outs[[1]]$growth_inc~params$mid_pts,type='l',las=1,lwd=2,xlab='Size (mm)',ylab='Post-molt size (mm)')
lines(true_growth_inc~params$mid_pts,lty=2,col=2,lwd=2)
mtext(side=1,"Size (mm)",line=2.)
mtext(side=2,"Post-molt size (mm)",line=2.5)
dev.off()
#==================================
# plot equilibrium size comps
#=====================================
offit<-8 # this is dumb and because of the number of projected months in bzero calcs

n_at_len<-outs[[1]]$n_at_len[(nrow(outs[[1]]$n_at_len)-11-offit):(nrow(outs[[1]]$n_at_len)-offit),]
survey_n<-sweep(n_at_len,2,outs[[1]]$survey_selectivity,"*")
est_eq_spawner_n<-sweep(survey_n,2,params$maturity_at_size,FUN="*")
est_eq_spawner_bio<-sweep(est_eq_spawner_n,2,params$wt_at_size,FUN="*")

colnames(est_eq_spawner_bio)<-outs[[1]]$sizes
in_g<-data.frame(melt((est_eq_spawner_bio)))
colnames(in_g)<-c("Month","Size","Numbers")

size_yr <- ggplot(dat=in_g) 
size_yr <- size_yr + geom_density_ridges(aes(x=Size, y=Month, height = Numbers,
                                             group = Month, 
                                             alpha=.9999),stat = "identity",scale=3,fill='royalblue3') +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  labs(x="Carapace width (mm)") +
  #scale_fill_manual(values=c("blue4"))+
  scale_y_continuous(name='Month',position='right')
size_yr

#============================================================
# calculate estimated SPR
#============================================================
# last year of spawning biomass
comp_n_at_len<-outs[[1]]$n_at_len[(params$all_months[length(params$all_months)]-11):params$all_months[length(params$all_months)],]
comp_survey_n<-sweep(comp_n_at_len,2,outs[[1]]$survey_selectivity,"*")
est_comp_spawner_n<-sweep(comp_survey_n,2,params$maturity_at_size,FUN="*")
est_comp_spawner_bio<-sweep(est_comp_spawner_n,2,params$wt_at_size,FUN="*")

eq_spbio_tot_month<-apply(est_eq_spawner_bio,1,sum)
comp_spbio_tot_month<-apply(est_comp_spawner_bio,1,sum)

em_spr<-comp_spbio_tot_month/eq_spbio_tot_month
boxplot(em_spr)

#============================================================
# compare the survey size composition data to the equilibrium
#============================================================


use_months<-outs[[1]]$survey_months
unq_yr<-use_months[length(use_months)]/12
for(x in 0:(unq_yr-1))
{
  #==select months
  mon_rng<-seq((x)*12+1,(x+1)*12)
  #==pull months inyear
  use_ind<-which(!is.na(match(use_months,mon_rng)))
  in_yr_mon<-use_months[use_ind]
  plot_months<-which(!is.na(match(mon_rng,in_yr_mon)))
  
  use_size_surv<-outs[[1]]$survey_size_comp_obs[use_ind,]
  use_n_surv<-outs[[1]]$survey_obs[use_ind]
  plot_n_surv<-sweep(use_size_surv,1,use_n_surv,"*")
  colnames(plot_n_surv)<-outs[[1]]$sizes
  rownames(plot_n_surv)<-plot_months
  
  in_f<-data.frame(melt((plot_n_surv)))
  colnames(in_f)<-c("Month","Size","Numbers")
  base_plot<-ggplot()+
    geom_line(data=in_g,aes(x=Size,y=Numbers))+
    theme_bw()+
    facet_wrap(~Month)+  
    geom_line(data=in_f,aes(x=Size,y=Numbers),col=2)
  png(paste('plots/size_obs_vs_equilibrium/fig_mon_',x,'.png',sep=''),height=8,width=8,res=350,units='in') 
  print(base_plot)
  dev.off()
 }
