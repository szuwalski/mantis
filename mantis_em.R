#===============================
# plot output from mantis model
#===============================
library(PBSmodelling)
#  read in EM output
rep_files<-c("admb/mantis.rep")
outs<-list(list())

for(x in 1:length(rep_files))
  outs[[x]]<-readList(rep_files[x])

#=======estimated processes vs. true================
par(mfrow=c(2,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
div_n<-100000000
plot_rec_est<-rep(0,length(all_months))
plot_rec_est[rec_months]<-outs[[1]]$recruits/div_n
plot(plot_rec~all_months,type='l',ylim=c(0,max(outs[[1]]$recruits,exp(in_log_rec))/div_n),las=1,
     ylab='',xaxt='n')
plot_rec_true<-rep(0,length(all_months))
plot_rec_true[rec_months]<-exp(in_log_rec)/div_n
lines(plot_rec_true,col=2)
legend("topleft",bty='n',col=c(1,2),lty=1,legend=c("Estimated","True"))
mtext(side=2,outer=F,line=2.75,"Recruits (hundred million)")

plot(outs[[1]]$applied_f_matrix[,4],type='l',ylim=c(0,max(outs[[1]]$applied_f_matrix)),
     las=1,ylab="Fishing mortality",xlab="Month")
lines(f_mort1[,1],col=2)
mtext(side=2,outer=F,line=2.75,"Fishing mortality")

#================fits to data=================
par(mfrow=c(2,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
plot(outs[[1]]$survey_obs/div_n~outs[[1]]$survey_months,
     ylim=c(0,max(outs[[1]]$survey_ob/div_n * exp(1.96*sqrt(log(1+surv_sd^2))))),pch=16,
     las=1,ylab="Survey numbers (hundred millions)",
     xaxt='n')
for(j in 1:length(outs[[1]]$survey_months))
{
  segments(x0=outs[[1]]$survey_months[j],x1=outs[[1]]$survey_months[j],
           y0=outs[[1]]$survey_ob[j]/div_n  /   exp(1.96*sqrt(log(1+surv_sd^2))),
           y1=outs[[1]]$survey_ob[j]/div_n * exp(1.96*sqrt(log(1+surv_sd^2))))
}
lines(outs[[1]]$survey_pred_comp/div_n~outs[[1]]$survey_months,lwd=3,col='blue')
legend("topleft",bty='n',legend="Survey")

plot(outs[[1]]$cat_obs/div_n~outs[[1]]$catch_months,
     ylim=c(0,max(outs[[1]]$cat_obs/div_n,outs[[1]]$catch_pred_comp/div_n)),pch=16,las=1,
     ylab="Catch (hundred millions)",
     xlab="Month")
for(j in 1:length(outs[[1]]$cat_obs))
{
  segments(x0=outs[[1]]$catch_months[j],x1=outs[[1]]$catch_months[j],
           y0=outs[[1]]$cat_obs[j]/div_n  /   exp(1.96*sqrt(log(1+cat_sd^2))),
           y1=outs[[1]]$cat_obs[j]/div_n * exp(1.96*sqrt(log(1+cat_sd^2))))
}
lines(outs[[1]]$catch_pred_comp/div_n~outs[[1]]$catch_months,lwd=2,col='red')
legend("topleft",bty='n',legend="Catch")
mtext(side=2,line=2.5,"Numbers (10000000)",outer=T)
mtext(side=1,line=2.5,"Month",outer=T)
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

#==catch size comps
par_in<-ceiling(sqrt(nrow(outs[[1]]$catch_size_comp_pred_comp)))
par(mfrow=c(par_in,par_in),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(outs[[1]]$catch_size_comp_pred_comp))
{
  plot(outs[[1]]$catch_size_comp_pred_comp[x,],yaxt='n',xaxt='n',ylim=c(0,0.6))
  lines(outs[[1]]$catch_size_comp_obs[x,],col=2,lty=3)
}


#==size transition matrix
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggridges)

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
