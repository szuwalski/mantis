#========================================
# implement OM across uncertainty axes
#=======================================
source("functions/om_function.R")
source("functions/scenario_plot_2.R")
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(patchwork)
library(cowplot)

params<-NULL
#############################
# common model characteristics
#############################
params$sizes         <-seq(50,250,5) # choose start sizes so they have roughly the same intermolt duration
params$years         <-seq(2000,2021) # simulated year span
params$tot_time_steps<-length(params$years)*12
all_months           <-seq(1,params$tot_time_steps)

#==mid points to the size bins from 'sizes'
in_break<-params$sizes
params$mid_pts<-rep(0,length(in_break)-1)
for(x in 1:length(params$mid_pts))
  params$mid_pts[x] <- (in_break[x]+in_break[x+1])/2

#==fishing and natural mortality parameters
params$nat_mort         <-1

#==growth
params$growth_sd       <-4
params$alpha           <-4
params$beta            <-1.15
#==maturity at size
params$mat_50           <-75
params$mat_slope        <-1
params$maturity_at_size <-1/(1+exp(-params$mat_slope*(params$mid_pts-params$mat_50)))

#==weight at size
params$wt_at_size    <- exp(2.85 * log(params$mid_pts) - 10.43)

#==fishery selectivity
params$fish_50       <-85
params$fish_slope    <-1
params$fish_sel      <-1/(1+exp(-params$fish_slope*(params$mid_pts-params$fish_50)))

#==molting probability
params$molt_50       <-150
params$molt_slope    <- -.05
params$molt_sel      <-1/(1+exp(-params$molt_slope*(params$mid_pts-params$molt_50)))

#==for simulated data and projections
params$size_comp_weight_surv <-150
params$size_comp_weight_catch<-150

#==when the surveys occur for simulated data
params$all_months<-seq(1,(params$tot_time_steps))
survey_sampling_1<-all_months[which(params$all_months%%8==0)]
survey_sampling_2<-all_months[which(params$all_months%%5==0)]
survey_sampling_3<-seq(96,108)
survey_sampling_3<-1

params$use_surv_samp<-sort(union(union(survey_sampling_1,survey_sampling_2),survey_sampling_3))

#==observation error for survey and catches
params$surv_sd  <-0.2
params$cat_sd   <-0.05

#==survey selectivity
params$surv_50   <-60
params$surv_slope<-.3
params$surv_sel  <-1/(1+exp(-params$surv_slope*(params$mid_pts-params$surv_50)))

#====================================
# parameters for complex life history
#====================================

#==timing of life history processes
                         #J,F,M,A,M,J,J,A,S,O,N,D
params$grow_month    <-c(0,0,0,1,0,0,1,0,0,0,0,0)
params$rec_month     <-c(1,0,0,0,1,0,0,0,1,0,0,0)
params$spawn_month   <-c(0,0,0,0,0,0,1,0,0,0,0,1)
params$fish_prop     <-c(0,0,.1,.1,.1,.1,0,0,0,.2,.2,.15) # this is the proportion of fish_mort that occurs in a month
params$fish_months   <-all_months[which(all_months*rep(params$fish_prop,length(params$years))!=0)]

#==calculate fishing months
params$fish_months   <-all_months[which(all_months*rep(params$fish_prop ,length(params$years))!=0)]

params$growth_months <-all_months[which(all_months*rep(params$grow_month,length(params$years))!=0)]
params$rec_months    <-all_months[which(all_months*rep(params$rec_month,length(params$years))!=0)]
params$spawn_months  <-all_months[which(all_months*rep(params$spawn_month,length(params$years))!=0)]

#==specify fishing mortality
f_mort1              <-matrix(0,ncol=length(params$mid_pts),nrow=params$tot_time_steps)
in_fmort             <-rnorm(length(params$fish_months),0.20001,0.000)
in_fmort[in_fmort<0.001] <-0
f_mort1[params$fish_months,]<-in_fmort
params$f_mort        <-sweep(f_mort1,2,params$fish_sel,FUN="*")

#==recruitment
#==rec_mu can be a vector if recruitment strength varies by recruitment month 
params$rec_mu<-c(20,20,20)
params$rec_sd<-0

#==add proportion of recruitment occuring in each season
params$proj_yr<-100

params$surv_sd<-0.2
params$cat_sd<-0.05
params$surv_size_samp<-10000
params$cat_size_samp<-10000

outs_complex_init<-operating_model(params,
                              sim_data=FALSE,
                              input_init=0)

outs_complex<-operating_model(params,
                              sim_data=TRUE,
                              dat_file="C:/Users/cody.szuwalski/Work/mantis/em/complex/mantis.DAT",
                              pin_file="C:/Users/cody.szuwalski/Work/mantis/em/complex/mantis.PIN",
                              input_init=1,input_init_N=outs_complex_init$eq_numbers[1,])

file.copy(from="admb/mantis.exe",to='em/complex')

png('plots/figure_2_complex.png',height=6,width=8.5,res=350,units='in')
scenario_plot(params,outs_complex)
dev.off()
#====================================
# parameters for simple life history
#====================================

#==timing of life history processes
#J,F,M,A,M,J,J,A,S,O,N,D
params_simp<-params
params_simp$grow_month    <-c(0,0,0,0,0,0,1,0,0,0,0,0)
params_simp$rec_month     <-c(0,0,0,0,0,0,1,0,0,0,0,0)
params_simp$spawn_month   <-c(0,0,0,0,0,0,1,0,0,0,0,0)
params_simp$fish_prop     <-rep(1/12,12) # this is the proportion of fish_mort that occurs in a month
params_simp$fish_months   <-all_months[which(all_months*rep(params_simp$fish_prop,length(params_simp$years))!=0)]

#==calculate fishing months
params_simp$fish_months   <-all_months[which(all_months*rep(params_simp$fish_prop ,length(params_simp$years))!=0)]
params_simp$growth_months <-all_months[which(all_months*rep(params_simp$grow_month,length(params_simp$years))!=0)]
params_simp$rec_months    <-all_months[which(all_months*rep(params_simp$rec_month,length(params_simp$years))!=0)]
params_simp$spawn_months  <-all_months[which(all_months*rep(params_simp$spawn_month,length(params_simp$years))!=0)]

#==specify fishing mortality
f_mort1              <-matrix(0,ncol=length(params_simp$mid_pts),nrow=params_simp$tot_time_steps)
in_fmort             <-rnorm(length(params_simp$fish_months),0.20001,0.0000)
in_fmort[in_fmort<0] <-0
f_mort1[params_simp$fish_months,]<-in_fmort
params_simp$f_mort        <-sweep(f_mort1,2,params_simp$fish_sel,FUN="*")

#==recruitment
#==rec_mu can be a vector if recruitment strength varies by recruitment month 
params_simp$rec_mu<-c(20)
params_simp$rec_sd<-0

params_simp$proj_yr<-100

outs_simple_init<-operating_model(params_simp,
                                   sim_data=FALSE,
                                   input_init=0)

outs_simple<-operating_model(params_simp,
                              sim_data=FALSE,
                              dat_file="C:/Users/cody.szuwalski/Work/mantis/em/simple/mantis.DAT",
                              pin_file="C:/Users/cody.szuwalski/Work/mantis/em/simple/mantis.PIN",
                              input_init=1,input_init_N=outs_simple_init$eq_numbers[1,])

png('plots/figure_2_simple.png',height=6,width=8.5,res=350,units='in')
scenario_plot(params_simp,outs_simple)
dev.off()
cbind(outs_simple$recruitment,
      outs_complex$recruitment)
#=======================================
# plot equilibrium size comps by model
#=======================================
colnames(outs_simple$eq_spawner_bio)<-params_simp$mid_pts
colnames(outs_complex$eq_spawner_bio)<-params$mid_pts
df1<-melt(outs_complex$eq_spawner_bio)
df2<-melt(outs_simple$eq_spawner_bio)
colnames(df1)<-c("Month","Size","Var")
colnames(df2)<-c("Month","Size","Var")
df1$model<-"Complex"
df2$model<-"Simple"
use_df<-rbind(df1,df2)

##==ggridges
p <- ggplot(data=use_df) 
eq_spbio <- p + geom_density_ridges(aes(x=Size, y=Month, height = Var, group = as.factor(Month), 
                                      fill=stat(y),alpha=.9999),stat = "identity",scale=1.25) +
  #scale_fill_viridis_c()+
  theme_bw() +
  facet_wrap(~model)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0)) +
  #text = element_text(size=20)) +
  #scale_fill_manual(values=c(rep(1,36),2,2,2))+
  labs(x="Body length (mm)")

png('plots/figure_equilibrium_sb.png',height=6,width=5,res=350,units='in')
eq_spbio
dev.off()

yurp<-ggplot(data=use_df)+
  geom_line(aes(x=Size,y=Var,group=Month))+
  facet_wrap(~model)+theme_bw()


#=======================================
# plot catch size comps by model
#=======================================
colnames(outs_simple$c_matrix)<-params_simp$mid_pts
colnames(outs_complex$c_matrix)<-params$mid_pts
df1<-melt(outs_complex$c_matrix)
df2<-melt(outs_simple$c_matrix)
colnames(df1)<-c("Month","Size","Var")
colnames(df2)<-c("Month","Size","Var")
df1$model<-"Complex"
df2$model<-"Simple"
use_df_cat<-rbind(df1,df2)

##==ggridges
p <- ggplot(data=use_df_cat) 
eq_catch <- p + geom_density_ridges(aes(x=Size, y=Month, height = Var, group = as.factor(Month), 
                                        fill=stat(y),alpha=.9999),stat = "identity",scale=1.25) +
  #scale_fill_viridis_c()+
  theme_bw() +
  facet_wrap(~model)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0)) +
  #text = element_text(size=20)) +
  #scale_fill_manual(values=c(rep(1,36),2,2,2))+
  labs(x="Body length (mm)")

eq_catch

yurp<-ggplot(data=use_df_cat)+
  geom_line(aes(x=Size,y=Var,group=Month))+
  facet_wrap(~model)+theme_bw()

#=======================================
# plot final year size comps by mdoel
#=======================================
use_simp<-outs_simple$c_matrix[(nrow(outs_simple$n_matrix)-11):nrow(outs_simple$n_matrix),]
#use_simp<-sweep(use_simp,2,params$surv_sel,FUN="*")
use_comp<-outs_complex$c_matrix[(nrow(outs_complex$n_matrix)-11):nrow(outs_complex$n_matrix),]
#use_comp<-sweep(use_comp,2,params$surv_sel,FUN="*")
colnames(use_simp)<-params$mid_pts
colnames(use_comp)<-params$mid_pts
df1<-melt(use_comp)
df2<-melt(use_simp)
colnames(df1)<-c("Month","Size","Var")
colnames(df2)<-c("Month","Size","Var")
df1$model<-"Complex"
df2$model<-"Simple"
use_df<-rbind(df1,df2)

##==ggridges
p <- ggplot(data=use_df) 
surv_comp <- p + geom_density_ridges(aes(x=Size, y=Month, height = Var, group = as.factor(Month), 
                                        fill=stat(y),alpha=.9999),stat = "identity",scale=1.25) +
  #scale_fill_viridis_c()+
  theme_bw() +
  facet_wrap(~model)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0)) +
  #text = element_text(size=20)) +
  #scale_fill_manual(values=c(rep(1,36),2,2,2))+
  labs(x="Body length (mm)")

surv_comp

yurp<-ggplot(data=use_df)+
  geom_line(aes(x=Size,y=Var,group=Month))+
  facet_wrap(~model)+theme_bw()



