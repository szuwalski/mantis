library(ggplot2)
library(dplyr)
library(reshape)
library(ggridges)
library(gghighlight)
library(gridExtra)

in_dat<-read.csv("data/size_comps.csv")

#==total numbers at size
melted<-melt(in_dat, id.vars = "size")
colnames(melted)<-c("Size","Month","Numbers")

##==ggridges
xlab <- paste0("\n", xlab)
p <- ggplot(data=melted) 
p_natl <- p + geom_density_ridges(aes(x=Size, y=Month, height = Numbers, group = as.factor(Month), 
                                      fill=stat(y),alpha=.9999),stat = "identity",scale=1.5) +
  #scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0)) +
  #text = element_text(size=20)) +
  #scale_fill_manual(values=c(rep(1,36),2,2,2))+
  labs(x="Body length (mm)")
print(p_natl)



######################################################
# size compositions data
#############################################
size_comp<-in_dat
tot_n<-apply(in_dat[,2:ncol(in_dat)],2,sum,na.rm=T)
size_comp[,2:ncol(in_dat)]<-sweep(in_dat[,2:ncol(in_dat)],2,tot_n,"/")
melted<-melt(size_comp, id.vars = "size")
colnames(melted)<-c("Size","Month","Numbers")

##==ggridges
xlab <- paste0("\n", xlab)
p <- ggplot(data=melted) 
p_natl_2 <- p + geom_density_ridges(aes(x=Size, y=Month, height = Numbers, group = as.factor(Month), 
                                      fill=stat(y),alpha=.9999),stat = "identity",scale=1.5) +
  #scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0)) +
  #text = element_text(size=20)) +
  #scale_fill_manual(values=c(rep(1,36),2,2,2))+
  labs(x="Body length (mm)")
print(p_natl_2)

library(patchwork)
patchwork <- p_natl_2 + p_natl 
patchwork[[2]] = patchwork[[2]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank() )
patchwork


#####################################
#==sum over month
#==total numbers at size
melted<-melt(in_dat, id.vars = "size")
colnames(melted)<-c("Size","Month","Numbers")

##==ggridgess
xlab <- paste0("\n", xlab)
p <- ggplot(data=melted) 
p_natl <- p + geom_density_ridges(aes(x=Size, y=Month, height = Numbers, group = as.factor(Month), 
                                      fill=stat(y),alpha=.9999),stat = "identity",scale=1.5) +
  #scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0)) +
  #text = element_text(size=20)) +
  #scale_fill_manual(values=c(rep(1,36),2,2,2))+
  labs(x="Body length (mm)")
print(p_natl)
