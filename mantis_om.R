# Operating model for mantis shrimp
# key characteristics: 
#     two recruitment events a year
#     two growth periods a year
#     fishery doesn't operate in the winter or summer
#     closed season changed length (need to implement)
#     monthly time step
#     start at 50mm
#     mature at 80mm
#     live 3 to 4 years

####################
# TO INCLUDE AT A LATER DATE
# include sex
# include probability of molting by size
# include a distribution of recruits to size bins
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggridges)
library(gridExtra)

#############################
# model characteristics
#############################
sizes         <-seq(50,190,10) # choose start sizes so they have roughly the same intermolt duration
years         <-seq(2000,2021) # simulated year span
tot_time_steps<-length(years)*12
all_months    <-seq(1,tot_time_steps)

#==mid points to the size bins from 'sizes'
in_break<-sizes
mid_pts<-rep(0,length(in_break)-1)
for(x in 1:length(mid_pts))
  mid_pts[x] <- (in_break[x]+in_break[x+1])/2

nat_mort      <-1

                 #J,F,M,A,M,J,J,A,S,O,N,D
grow_month    <-c(0,0,1,0,0,0,1,0,0,0,0,0)
rec_month     <-c(0,1,0,0,0,0,0,1,0,0,0,0)
spawn_month   <-c(0,0,0,0,0,0,1,0,0,0,0,1)
fish_prop     <-c(0,0,.1,.1,.1,.1,0,0,0,.2,.2,.15) # this is the proportion of fish_mort that occurs in a month
fish_months   <-all_months[which(all_months*rep(fish_prop,length(years))!=0)]

grow_month    <-c(0,0,0,0,0,0,1,0,0,0,0,0)
rec_month     <-c(0,0,0,0,0,0,0,1,0,0,0,0)
spawn_month   <-c(0,0,0,0,0,0,1,0,0,0,0,1)
fish_prop     <-c(.1,.1,.1,.1,.1,.1,.10,.10,.10,.2,.2,.15) # this is the proportion of fish_mort that occurs in a month
fish_months   <-all_months[which(all_months*rep(fish_prop,length(years))!=0)]


growth_months   <-all_months[which(all_months*rep(grow_month,length(years))!=0)]
rec_months    <-all_months[which(all_months*rep(rec_month,length(years))!=0)]
spawn_months  <-all_months[which(all_months*rep(spawn_month,length(years))!=0)]

#==maturity at size
mat_50<-75
mat_slope<-1
maturity_at_size<-1/(1+exp(-mat_slope*(mid_pts-mat_50)))
plot(maturity_at_size~mid_pts,type='b',las=1)

#==weight at size
wt_at_size <- exp(2.85 * log(mid_pts) - 10.43)
plot(wt_at_size~mid_pts)

#==fishery selectivity
fish_50<-85
fish_slope<-1
fish_sel<-1/(1+exp(-fish_slope*(mid_pts-fish_50)))
plot(fish_sel~mid_pts,type='b')

f_mort1              <-matrix(0,ncol=length(mid_pts),nrow=tot_time_steps)
in_fmort             <-rnorm(length(fish_months),0.15,0.05)
in_fmort[in_fmort<0]<-0
#in_fmort<-0
f_mort1[fish_months,]<-in_fmort
f_mort              <-sweep(f_mort1,2,fish_sel,FUN="*")

###################################
#==build a size transition matrix==
make_size_trans<-function(alpha,beta,growth_sd,out_plot=FALSE,no_x=FALSE,
                          input_size)
{
  avg_post_molt<-alpha + beta*input_size
  size_trans_mat<-matrix(ncol=length(input_size),nrow=length(input_size))
  
  for(x in 1:nrow(size_trans_mat))
  {
    tmp<-dnorm(input_size,avg_post_molt[x],growth_sd[x])  
    tmp[seq(1,nrow(size_trans_mat))<x]<-0
    size_trans_mat[x,]<-round(tmp/sum(tmp) ,3)
    
  }
  
  rownames(size_trans_mat)<-input_size
  colnames(size_trans_mat)<-input_size
  in_g<-data.frame(melt(t(size_trans_mat)))
  colnames(in_g)<-c("Postmolt","Premolt","Density")
  
  p <- ggplot(in_g)
  p <- p + geom_density_ridges(aes(x=Postmolt, y=Premolt, height = Density, group=Premolt,
                                   fill=stat(y),alpha=.9999),stat = "identity",scale=3) +
    scale_fill_viridis_c()+
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90)) 
  
  if(no_x==TRUE)
    p<-p+theme(axis.text.x=element_blank(),
               axis.title.x=element_blank(),
               axis.ticks.x = element_blank())
  
  if(out_plot==TRUE)
    print(p)
  
  list(size_trans_mat,p)
}

size_trans_1<-make_size_trans(alpha=4,beta=1.2,
                              growth_sd=rep(5,length(mid_pts)),out_plot=F,no_x=FALSE,input_size=mid_pts)
# plot size transition matrix
size_trans_1[[2]]


##########################################
#==make up initial numbers at size at time
#==by iterating the size transition matrix
#==on a recruitment
init_numbers<-matrix(ncol=length(mid_pts),nrow=4)
dummy<-c(exp(10),rep(0,length(mid_pts)-1))
counter<-0
for(x in 1:4)
{
 counter<-0
 while(counter<x)
 {
   dummy<-(dummy%*%size_trans_1[[1]])*exp(-nat_mort/length(grow_month))
   counter<-counter+1
 }
 init_numbers[x,]<-dummy
   
}
#==check they look sensible
 plot(init_numbers[1,]) 
 for(x in 1:nrow(init_numbers))
   lines(init_numbers[x,])

 ###############################
 #==project population
 ###############################
 n_matrix<-matrix(ncol=length(mid_pts),nrow=c(tot_time_steps))
 n_matrix[1,]<-apply(init_numbers,2,sum)
 c_matrix<-matrix(ncol=length(mid_pts),nrow=c(tot_time_steps))
 log_recruits<-rep(0,tot_time_steps)
 in_log_rec<-rnorm(length(rec_months),20,1)
 log_recruits[rec_months] <- in_log_rec
   
for(x in 1:(tot_time_steps))
{
 #==recruitment(this could be spread over several size bins)
  #==this can also be linked to spawning biomass later
  if(any(rec_months==x))
    n_matrix[x,1]<-n_matrix[x,1] + exp(log_recruits[x])
  
 #==growth
  if(any(growth_months==x))
    n_matrix[x,]<-n_matrix[x,]%*%size_trans_1[[1]]
  
  #==fishing
  #==need a better way to do this...
    n_matrix[x,]<-n_matrix[x,]*(exp(-f_mort[x,]))
    c_matrix[x,]<-n_matrix[x,]*(1-exp(-f_mort[x,]))   

  #==natural mortality
    if(x<tot_time_steps)
     n_matrix[x+1,]<-n_matrix[x,]*exp(-nat_mort/12)

}
 ################################################################### 
 # Rerun with no fishing and constant recruitment set at the average
 ################################################################### 
 proj_yr<-100
 proj_t_step<-proj_yr*12
 proj_months<-seq(1,proj_t_step)
 n_matrix_b0<-matrix(ncol=length(mid_pts),nrow=c(proj_t_step))
 n_matrix_b0[1,]<-apply(init_numbers,2,sum)
 
 proj_growth_months <-proj_months[which(proj_months*rep(grow_month,length(proj_yr))!=0)]
 proj_rec_months    <-proj_months[which(proj_months*rep(rec_month,length(proj_yr))!=0)]
 proj_spawn_months  <-proj_months[which(proj_months*rep(spawn_month,length(proj_yr))!=0)]
 
 for(x in 1:(proj_t_step))
 {
   #==recruitment(this could be spread over several size bins)
   #==this can also be linked to spawning biomass later
   if(any(proj_rec_months==x))
     n_matrix_b0[x,1]<-n_matrix_b0[x,1] + exp(mean(log_recruits[log_recruits!=0]))
   
   #==growth
   if(any(proj_growth_months==x))
     n_matrix_b0[x,]<-n_matrix_b0[x,]%*%size_trans_1[[1]]
   
   #==natural mortality
   if(x<proj_t_step)
     n_matrix_b0[x+1,]<-n_matrix_b0[x,]*exp(-nat_mort/12)
   
 }
 
 eq_numbers<-n_matrix_b0[(proj_t_step-11):proj_t_step,]
 eq_spawner_n<-sweep(eq_numbers,2,maturity_at_size,FUN="*")
 eq_spawner_bio<-sweep(eq_spawner_n,2,wt_at_size,FUN="*")

 ####################################
 # visualize equilibrium population
 ################################
 rownames(eq_numbers)<-seq(1,12)
 colnames(eq_numbers)<-mid_pts
 in_g<-data.frame(melt(t(eq_numbers)))
 colnames(in_g)<-c("Size","Month","Numbers")
 in_g$year<-floor(in_g$Month)
 
 eq_plot_n <- ggplot(in_g)
 eq_plot_n <- eq_plot_n + geom_density_ridges(aes(x=Size, y=Month, height = Numbers, group=Month,
                                  fill=stat(y),alpha=.9999),stat = "identity",scale=5) +
   scale_fill_viridis_c()+
   theme_bw() +
   theme(legend.position = "none",
         axis.text.x = element_text(angle = 90)) 
 
 print(eq_plot_n)
 
 ####################################
 # visualize equilibrium spawners
 ################################
 rownames(eq_spawner_n)<-seq(1,12)
 colnames(eq_spawner_n)<-mid_pts
 in_g<-data.frame(melt(t(eq_spawner_n)))
 colnames(in_g)<-c("Size","Month","Numbers")
 in_g$year<-floor(in_g$Month)
 
 eq_plot <- ggplot(in_g)
 eq_plot <- eq_plot + geom_density_ridges(aes(x=Size, y=Month, height = Numbers, group=Month,
                                              fill=stat(y),alpha=.9999),stat = "identity",scale=5) +
   scale_fill_viridis_c()+
   theme_bw() +
   theme(legend.position = "none",
         axis.text.x = element_text(angle = 90)) 
 
 print(eq_plot)
 
 
###########################################
# compare equilibrium spawners to observed
##########################################
 
 n_spawner<-sweep(n_matrix,2,maturity_at_size,FUN="*")
 bio_spawner<-sweep(n_spawner,2,wt_at_size,FUN="*")
 spr_calc<-bio_spawner
 spr_calc_tot<-rep(0,nrow(bio_spawner))
 loop_seq<-seq(1,nrow(bio_spawner))%%12
 loop_seq[loop_seq==0]<-12
 for(x in 1:12)
 {
   inds<-which(loop_seq==x)
   for(y in 1:length(inds))
   {
    spr_calc[inds[y],]<-bio_spawner[inds[y],]/eq_spawner_bio[x,]
    spr_calc_tot[inds[y]]<-sum(bio_spawner[inds[y],])/sum(eq_spawner_bio[x,])
   }
 }
 hist(spr_calc) 
 hist(spr_calc_tot)   
 
 
###################################################################
# Rerun with no fishing to find dynamic unfished biomass
###################################################################
 n_matrix_dynb0<-matrix(ncol=length(mid_pts),nrow=c(tot_time_steps))
 n_matrix_dynb0[1,]<-apply(init_numbers,2,sum)

 for(x in 1:(tot_time_steps))
 {
   #==recruitment(this could be spread over several size bins)
   #==this can also be linked to spawning biomass later
   if(any(rec_months==x))
     n_matrix_dynb0[x,1]<-n_matrix_dynb0[x,1] + exp(log_recruits[x])
   
   #==growth
   if(any(growth_months==x))
     n_matrix_dynb0[x,]<-n_matrix_dynb0[x,]%*%size_trans_1[[1]]
   
   #==natural mortality
   if(x<tot_time_steps)
     n_matrix_dynb0[x+1,]<-n_matrix_dynb0[x,]*exp(-nat_mort/12)
   
 }


 
####################################
# visualize population
################################
 rownames(n_matrix)<-seq(1,tot_time_steps)/12-0.0001
 colnames(n_matrix)<-mid_pts
 in_g<-data.frame(melt(t(n_matrix)))
 colnames(in_g)<-c("Size","Time","Numbers")
 in_g$year<-floor(in_g$Time)
 
 p <- ggplot(in_g)
 p <- p + geom_density_ridges(aes(x=Size, y=Time, height = Numbers, group=Time,
                                  fill=stat(y),alpha=.9999),stat = "identity",scale=15) +
   scale_fill_viridis_c()+
   theme_bw() +
   theme(legend.position = "none",
         axis.text.x = element_text(angle = 90)) 

   print(p)
   
   num_yr<-in_g%>%
     group_by(year)%>%
     summarize(tot_num = sum(Numbers,na.rm=T))
   
  q<- ggplot(data=num_yr,aes(x=year,y=tot_num))+
     geom_line(lwd=1.4)+
     theme_bw()+expand_limits(y=0)
   print(q)
   grid.arrange(p,q,layout_matrix=matrix(c(1,2,2),ncol=3))
   
   
####################################
# visualize catch
################################
rownames(c_matrix)<-seq(1,tot_time_steps)/12
colnames(c_matrix)<-mid_pts
in_g<-data.frame(melt(t(c_matrix)))
colnames(in_g)<-c("Size","Time","Numbers")
in_g$year<-floor(in_g$Time)
   
   p <- ggplot(in_g)
   p <- p + geom_density_ridges(aes(x=Size, y=Time, height = Numbers, group=Time,
                                    fill=stat(y),alpha=.9999),stat = "identity",scale=5) +
     scale_fill_viridis_c()+
     theme_bw() +
     theme(legend.position = "none",
           axis.text.x = element_text(angle = 90)) 
   
   print(p)
 
cat_yr<-in_g%>%
  group_by(year)%>%
  summarize(tot_catch = sum(Numbers,na.rm=T))
 
q<- ggplot(data=cat_yr,aes(x=year,y=tot_catch))+
  geom_line(lwd=1.4)+
  theme_bw()+expand_limits(y=0)
print(q)
grid.arrange(p,q,layout_matrix=matrix(c(1,2,2),ncol=3))

###########################
# write a .DAT file
###########################  

dat_file<-'admb/mantis_pindat.DAT'
file.create(dat_file)

cat("# log_rec",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(in_log_rec,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat("# fmort",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(in_fmort,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat("# init_n_at_l",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(log(n_matrix[1,]),file=dat_file,append=TRUE) 



all_months<-seq(1,(tot_time_steps))
survey_sampling_1<-all_months[which(all_months%%8==0)]
survey_sampling_2<-all_months[which(all_months%%5==0)]
survey_sampling_3<-seq(96,108)
survey_sampling_3<-1
use_surv_samp<-sort(union(union(survey_sampling_1,survey_sampling_2),survey_sampling_3))
surv_sd<-0.2
cat_sd<-0.05

#==survey selectivity
surv_50<-60
surv_slope<-.3
surv_sel<-1/(1+exp(-surv_slope*(mid_pts-surv_50)))
#plot(surv_sel,type='b')

#==calculate fishing months
fish_months<-all_months[which(all_months*rep(fish_prop,length(years))!=0)]

dat_file<-'admb/mantis.DAT'
file.create(dat_file)

cat("# simulated mantis shrimp data file",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat("# start_mo",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(1,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat("# end_mo",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(tot_time_steps,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of sizes",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(mid_pts),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# sizes",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(mid_pts,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of survey months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(use_surv_samp),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(use_surv_samp,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey numbers",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
surv_obs<-sweep(n_matrix,2,surv_sel,FUN="*")[use_surv_samp,]
surv_obs_error<-apply(surv_obs,1,sum)*exp(rnorm(nrow(surv_obs),0,surv_sd))
cat(surv_obs_error,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey size composition",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
for(x in 1:nrow(surv_obs))
{
  size_comp<-hist(sample(mid_pts,size=10000,prob=surv_obs[x,],replace=TRUE),plot=FALSE,breaks=sizes)$density
  cat(size_comp/sum(size_comp),file=dat_file,append=TRUE)
  cat("\n",file=dat_file,append=TRUE)
}

cat("# number of catch months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(fish_months),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

# cat("# which catch months",file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# cat(which((fish_prop!=0)),file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)

cat("# all catch months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(fish_months,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# catch numbers",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
catch_obs<-c_matrix[fish_months,]
catch_obs_error<-apply(c_matrix[fish_months,],1,sum)*exp(rnorm(length(fish_months),0,cat_sd))
cat(catch_obs_error,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# catch size composition",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
for(x in 1:nrow(catch_obs))
{
  size_comp<-hist(sample(mid_pts,size=10000,prob=catch_obs[x,],replace=TRUE),plot=FALSE,breaks=sizes)$density
  cat(size_comp/sum(size_comp),file=dat_file,append=TRUE)
  cat("\n",file=dat_file,append=TRUE)
}

# cat("# number of recruitment months",file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# cat(sum(rec_month),file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# 
# cat("# which recruitment months",file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# cat(which(rec_month!=0),file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)

cat("# number of total recruitment months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(rec_months),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# all recruitment months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(rec_months,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of growth months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(sum(grow_month),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# which growth months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(which(grow_month!=0),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# catch cv",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(cat_sd,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey cv",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(surv_sd,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# catch eff N",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(150,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey eff N",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(150,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
# cat("# number of total growth months",file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# cat(length(growth_months),file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# 
# cat("# all growth months",file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# cat(growth_months,file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
