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
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggridges)

#############################
# model characteristics
#############################
sizes         <-seq(50,190,10) # choose start sizes so they have roughly the same intermolt duration
growth_month  <-c(7,3)         # month in which growth occurs
rec_month     <-c(8,2)         # month in which recruitment occurs
spawning_month<-c(6,12)        # months in which spawning occurs (not actually used yet)
years         <-seq(2009,2021) # simulated year span
total_time_steps<-length(years)*12

nat_mort      <-1
fmort         <-rep(1,length(years))
fish_month    <-c(3,4,5,6,10,11,12)                 # natural mortality
fish_prop     <-c(0,0,.1,.1,.1,.05,0,0,0,.3,.2,.15) # this is the proportion of fish_mort that occurs in a month
fish_mort     <-rnorm(length(years),0.05,0.5)       # this is the fishing mortality by year...have to think about best way to do this

#==mid points to the size bins from 'sizes'
in_break<-sizes
mid_pts<-rep(0,length(in_break)-1)
for(x in 1:length(mid_pts))
  mid_pts[x] <- (in_break[x]+in_break[x+1])/2

#==fishery selectivity
fish_50<-65
fish_slope<-1
fish_sel<-1/(1+exp(-fish_slope*(mid_pts-fish_50)))
#==check it
plot(fish_sel,type='b')

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

size_trans_1<-make_size_trans(alpha=4,beta=1.1,
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
   dummy<-(dummy%*%size_trans_1[[1]])*exp(-nat_mort/length(growth_month))
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
 n_matrix<-matrix(ncol=length(mid_pts),nrow=c(total_time_steps))
 n_matrix[1,]<-apply(init_numbers,2,sum)
 c_matrix<-matrix(ncol=length(mid_pts),nrow=c(total_time_steps))
   
for(x in 1:(total_time_steps-1))
{
 #==recruitment(this could be spread over several size bins)
  #==this can also be linked to spawning biomass later
  if(any(x%%rec_month==0))
    n_matrix[x,1]<-n_matrix[x,1] + exp(rnorm(1,10,1))
  
 #==growth
  if(any(x%%growth_month==0))
    n_matrix[x,]<-n_matrix[x,]%*%size_trans_1[[1]]
  
  #==fishing
  #==need a better way to do this...
  if(any(x%%fish_month==0))
  {
    n_matrix[x,]<-n_matrix[x,]*(exp(-fish_sel*fmort[x%%12+1]*fish_prop[x%%12+1]))
    c_matrix[x,]<-n_matrix[x,]*(1-exp(-fish_sel*fmort[x%%12+1]*fish_prop[x%%12+1]))   
  }
  
  #==natural mortality
  n_matrix[x+1,]<-n_matrix[x,]*exp(-nat_mort/12)

}

####################################
# visualize population
################################
 rownames(n_matrix)<-seq(1,total_time_steps)/12
 colnames(n_matrix)<-mid_pts
 in_g<-data.frame(melt(t(n_matrix)))
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
   
   num_yr<-in_g%>%
     group_by(year)%>%
     summarize(tot_num = sum(Numbers,na.rm=T))
   
   ggplot(data=num_yr,aes(x=year,y=tot_num))+
     geom_line()+
     theme_bw()
   
####################################
# visualize catch
################################
rownames(c_matrix)<-seq(1,total_time_steps)/12
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
 
ggplot()+
  geom_line(data=num_yr,aes(x=year,y=tot_num))+
  geom_line(data=cat_yr,aes(x=year,y=tot_catch),col='red')+
  theme_bw()  
 
###########################
# next steps...
###########################  
