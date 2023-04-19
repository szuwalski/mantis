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

operating_model<-function(params,sim_data=FALSE,dat_file=NA,pin_file=NA)
{
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
          axis.text.x = element_text(angle = 90)) +
    ylab("Premolt size")+xlab("Postmolt size")
  
  if(no_x==TRUE)
    p<-p+theme(axis.text.x=element_blank(),
               axis.title.x=element_blank(),
               axis.ticks.x = element_blank())
  
  if(out_plot==TRUE)
    print(p)
  
  list(size_trans_mat,p)
}

size_trans_1<-make_size_trans(alpha=params$alpha,
                              beta=params$beta,
                              growth_sd=rep(params$growth_sd,length(params$mid_pts)),
                              out_plot=F,no_x=FALSE,input_size=params$mid_pts)

##########################################
#==make up initial numbers at size at time
#==by iterating the size transition matrix
#==on a recruitment
init_numbers<-matrix(ncol=length(params$mid_pts),nrow=4)
dummy<-c(exp(10),rep(0,length(params$mid_pts)-1))
counter<-0
for(x in 1:4)
{
  counter<-0
  while(counter<x)
  {
    dummy<-(dummy%*%size_trans_1[[1]])*exp(-params$nat_mort/length(params$grow_month))
    counter<-counter+1
  }
  init_numbers[x,]<-dummy
  
}

###############################
#==project population
###############################
n_matrix<-matrix(ncol=length(params$mid_pts),nrow=c(params$tot_time_steps))
n_matrix[1,]<-apply(init_numbers,2,sum)
c_matrix<-matrix(ncol=length(params$mid_pts),nrow=c(params$tot_time_steps))

#==set up for autocorrelated too and differences in intensity per
log_recruits<-rep(0,params$tot_time_steps)
in_log_rec<-rnorm(length(params$rec_months),params$rec_mu,params$rec_sd)
log_recruits[params$rec_months] <- in_log_rec

for(x in 1:(params$tot_time_steps))
{
  #==recruitment(this could be spread over several size bins)
  #==this can also be linked to spawning biomass later
  if(any(params$rec_months==x))
    n_matrix[x,1]<-n_matrix[x,1] + exp(log_recruits[x])
  
  #==growth
  if(any(params$growth_months==x))
  {
    #==molters
    temp_molt<-n_matrix[x,]*params$molt_sel
    #==non-molters
    temp_no_molt<-n_matrix[x,]*(1-params$molt_sel)
    
    n_matrix[x,]<-temp_molt%*%size_trans_1[[1]]  +  temp_no_molt
  }
  
  #==fishing
  #==need a better way to do this...
  n_matrix[x,]<-n_matrix[x,]*(exp(-params$f_mort[x,]))
  c_matrix[x,]<-n_matrix[x,]*(1-exp(-params$f_mort[x,]))   
  
  #==natural mortality
  if(x<params$tot_time_steps)
    n_matrix[x+1,]<-n_matrix[x,]*exp(-params$nat_mort/12)
  
}

################################################################### 
# Rerun with no fishing and constant recruitment set at the average
################################################################### 

proj_t_step<-params$proj_yr*12
proj_months<-seq(1,proj_t_step)
n_matrix_b0<-matrix(ncol=length(params$mid_pts),nrow=c(proj_t_step))
n_matrix_b0[1,]<-apply(init_numbers,2,sum)

proj_growth_months        <-proj_months[which(proj_months*rep(params$grow_month,params$proj_yr)!=0)]
proj_rec_months           <-proj_months[which(proj_months*rep(params$rec_month,params$proj_yr)!=0)]

for(x in 1:(proj_t_step))
{
  #==recruitment(this could be spread over several size bins)
  #==this can also be linked to spawning biomass later
  if(any(proj_rec_months==x))
    n_matrix_b0[x,1]<-n_matrix_b0[x,1] + exp(mean(log_recruits[ which(proj_months%%12==x%%12)],na.rm=T))
  
  #==growth
  if(any(proj_growth_months==x))
  {
    #==molters
    temp_molt<-n_matrix_b0[x,]*params$molt_sel
    #==non-molters
    temp_no_molt<-n_matrix_b0[x,]*(1-params$molt_sel)
    
    n_matrix_b0[x,]<-temp_molt%*%size_trans_1[[1]]  +  temp_no_molt
  }
  
  #==natural mortality
  if(x<proj_t_step)
    n_matrix_b0[x+1,]<-n_matrix_b0[x,]*exp(-params$nat_mort/12)
  
}

eq_numbers<-n_matrix_b0[(proj_t_step-11):proj_t_step,]
eq_spawner_n<-sweep(eq_numbers,2,params$maturity_at_size,FUN="*")
eq_spawner_bio<-sweep(eq_spawner_n,2,params$wt_at_size,FUN="*")


###################################################################
# Rerun with no fishing to find dynamic unfished biomass
###################################################################
n_matrix_dynb0<-matrix(ncol=length(params$mid_pts),nrow=c(params$tot_time_steps))
n_matrix_dynb0[1,]<-apply(init_numbers,2,sum)

for(x in 1:(params$tot_time_steps))
{
  #==recruitment(this could be spread over several size bins)
  #==this can also be linked to spawning biomass later
  if(any(params$rec_months==x))
    n_matrix_dynb0[x,1]<-n_matrix_dynb0[x,1] + exp(log_recruits[x])
  
  #==growth
  if(any(params$growth_months==x))
  {
    #==molters
    temp_molt<-n_matrix_dynb0[x,]*params$molt_sel
    #==non-molters
    temp_no_molt<-n_matrix_dynb0[x,]*(1-params$molt_sel)
    
    n_matrix_dynb0[x,]<-temp_molt%*%size_trans_1[[1]]  +  temp_no_molt
  }
    
  #==natural mortality
  if(x<params$tot_time_steps)
    n_matrix_dynb0[x+1,]<-n_matrix_dynb0[x,]*exp(-params$nat_mort/12)
  
}

if(sim_data==TRUE)
{
###########################
# write a .DAT file
###########################  

file.create(pin_file)
  
cat("# init_n_at_l",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(log(n_matrix[1,]),file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat("# log_rec",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(in_log_rec,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# nat mort",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$nat_mort,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# surv q",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(0.5,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# surv sel 50",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$surv_50,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# surv sel slope",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$surv_slope,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# fish sel 50",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$fish_50,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# fish sel slope",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$fish_slope,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# fmort",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(in_fmort,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
 
cat("# growth alpha",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$alpha,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# growth beta",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(0.4,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# growth slope",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$beta,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)


cat("# molt 50",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$molt_50,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# molt slope",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(params$molt_slope,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

cat("# dummy",file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)
cat(1,file=pin_file,append=TRUE)
cat("\n",file=pin_file,append=TRUE)

file.create(dat_file)

cat("# simulated mantis shrimp data file",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat("# start_mo",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(1,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat("# end_mo",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$tot_time_steps,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of sizes",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(params$mid_pts),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# sizes",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$mid_pts,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of survey months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(params$use_surv_samp),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$use_surv_samp,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey numbers",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
surv_obs<-sweep(n_matrix,2,params$surv_sel,FUN="*")[params$use_surv_samp,]
surv_obs_error<-apply(surv_obs,1,sum)*exp(rnorm(nrow(surv_obs),0,params$surv_sd))
cat(surv_obs_error,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey size composition",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
for(x in 1:nrow(surv_obs))
{
  size_comp<-hist(sample(params$mid_pts,size=params$surv_size_samp,prob=surv_obs[x,],replace=TRUE),plot=FALSE,breaks=params$sizes)$density
  cat(size_comp/sum(size_comp),file=dat_file,append=TRUE)
  cat("\n",file=dat_file,append=TRUE)
}

cat("# number of catch months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(params$fish_months),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

# cat("# which catch months",file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# cat(which((fish_prop!=0)),file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)

cat("# all catch months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$fish_months,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# catch numbers",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
catch_obs<-c_matrix[params$fish_months,]
catch_obs_error<-apply(c_matrix[params$fish_months,],1,sum)*exp(rnorm(length(params$fish_months),0,params$cat_sd))
cat(catch_obs_error,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# catch size composition",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
for(x in 1:nrow(catch_obs))
{
  size_comp<-hist(sample(params$mid_pts,size=params$cat_size_samp,prob=catch_obs[x,],replace=TRUE),plot=FALSE,breaks=params$sizes)$density
  cat(size_comp/sum(size_comp),file=dat_file,append=TRUE)
  cat("\n",file=dat_file,append=TRUE)
}

# cat("# number of recruitment months",file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# cat(sum(params$rec_month),file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# 
# cat("# which recruitment months",file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)
# cat(which(params$rec_month!=0),file=dat_file,append=TRUE)
# cat("\n",file=dat_file,append=TRUE)

cat("# number of total recruitment months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(params$rec_months),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# all recruitment months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$rec_months,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of rec months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(sum(params$rec_month),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# which rec months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(which(params$rec_month!=0),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of growth months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(sum(params$grow_month),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# which growth months",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(which(params$grow_month!=0),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# catch cv",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$cat_sd,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey cv",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$surv_sd,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# catch eff N",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$size_comp_weight_catch,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey eff N",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$size_comp_weight_surv,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
}

list(size_trans=size_trans_1,n_matrix=n_matrix,c_matrix=c_matrix,
     eq_numbers=eq_numbers,eq_spawner_n=eq_spawner_n,eq_spawner_bio=eq_spawner_bio,
     n_matrix_dynb0=n_matrix_dynb0,recruitment=log_recruits)
}