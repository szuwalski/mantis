#===============================
# plot output from mantis model
#===============================
library(PBSmodelling)
#  read in EM output
rep_files<-c("C:/mantis/admb/mantis.rep")
outs<-list(list())

for(x in 1:length(rep_files))
  outs[[x]]<-readList(rep_files[x])

names(outs[[1]])

plot(outs[[1]]$recruits,type='l',ylim=c(0,max(outs[[1]]$recruits)))
plot(outs[[1]]$applied_f_matrix[,4],type='l',ylim=c(0,max(outs[[1]]$applied_f_matrix)))


plot(outs[[1]]$survey_obs~outs[[1]]$survey_months,
     ylim=c(0,max(outs[[1]]$survey_pred_comp,outs[[1]]$survey_obs)))
lines(outs[[1]]$survey_pred_comp~outs[[1]]$survey_months,type='b',pch=16)

plot(outs[[1]]$cat_obs~outs[[1]]$catch_months,
     ylim=c(0,max(outs[[1]]$cat_obs,outs[[1]]$catch_pred_comp)))
lines(outs[[1]]$catch_pred_comp~outs[[1]]$catch_months,type='b',pch=16)

#==survey size comps
par_in<-ceiling(sqrt(nrow(outs[[1]]$survey_size_comp_obs)))
par(mfrow=c(par_in,par_in),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(outs[[1]]$survey_size_comp_pred_comp))
{
  plot(outs[[1]]$survey_size_comp_pred_comp[x,],yaxt='n',xaxt='n')
  lines(outs[[1]]$survey_size_comp_obs[x,],col=2,lty=3)
}

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
