#==make maps of study area with catch history and survey timings
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(autoimage)
library(akima)  #install.packages("akima")
library(viridis)
library(broom)
library(maps)
library("rnaturalearth")
library(interp)
library(RColorBrewer)
library(reshape2) # for melt
library(mgcv)  
library(PBSmapping)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(metR)
library(ggplot2)
#install.packages("rnaturalearthdata")
library(reshape)
library(dplyr)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(png)
library(grid)
library(PBSmodelling)
library(patchwork)
library(cowplot)
library(dplyr)
library(ggplot2)
library(dplyr)
library(reshape)
library(ggridges)
library(gghighlight)
library(gridExtra)
library(sf)

annotation_custom2 <-   function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, data){ layer(data = data, 
                                                                                                      stat = StatIdentity, 
                                                                                                      position = PositionIdentity,
                                                                                                      geom = ggplot2:::GeomCustomAnn,
                                                                                                      inherit.aes = TRUE, 
                                                                                                      params = list(grob = grob,xmin = xmin, xmax = xmax,
                                                                                                                    ymin = ymin, ymax = ymax))}


bohai<-read.csv("data/bohai.csv")
world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<- 117.51
lon_2<- 122.501
lat_1<-min(bohai$Lat,na.rm=T)
lat_2<-41.01
in_fin<-melt(filter(bohai,Species=='Oratosquilla oratoria')[,c(2,3,7,9,10)],id=c('Lat','Lon','Year',"Month"))
in_fin<-in_fin[complete.cases(in_fin),]

sf_use_s2(FALSE)
china_map<-ggplot() +
  geom_sf(data=world) +
  coord_sf(ylim = c(10,50), xlim = c(90,140), expand = FALSE) +
  geom_rect(aes(xmin=117,xmax=123,ymin=37,ymax=41),color='red',fill=NA)+
  theme_map()

inset_CH<-ggplotGrob(china_map)
#==plot distribution in Bohai Sea
op_map<-ggplot() + 
  #geom_tile(data=in_fin, aes(x = Lon, y = Lat, fill = value),width=.5,height=.25) +
  geom_point(data=in_fin, aes(x = Lon, y = Lat)) +
  # geom_point(data = gam_dat, 
  #            aes(x = lon, y = lat,fill=log_abund_101), 
  #            shape = 16,size=.45) +
  # scale_fill_distiller(palette="Spectral", na.value="grey") +
  # facet_wrap(~Year) +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10))+
  theme(legend.position=c(.88,.82),
        legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")+
  annotation_custom2(rasterGrob(in_png_mantis, interpolate=TRUE), 
                     xmin=120, xmax=123, ymin=37, ymax=39.25,data=in_dat[1:1280,])

study_area <- op_map +
  annotation_custom(grob = inset_CH, xmin = 117.11, xmax = 120,
                    ymin = 39.5, ymax = 41) 
study_area

in_dat<-read.csv("data/size_comps.csv")
in_png_mantis<-readPNG('presentations/mantis_grey.png')
#==total numbers at size
melted<-melt(in_dat, id.vars = "size")
colnames(melted)<-c("Size","Month","Numbers")

##==ggridges
xlab <- paste0("\n", xlab)
p <- ggplot(data=melted) 
p_natl <- p + geom_density_ridges(aes(x=Size, y=Month, height = Numbers, group = as.factor(Month), 
                                      fill=stat(y),alpha=.9999),stat = "identity",scale=1.25) +
  #scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0)) +
  #text = element_text(size=20)) +
  #scale_fill_manual(values=c(rep(1,36),2,2,2))+
  labs(x="Body length (mm)")
  
#==plot recent time series of catch (+ timing of fishery?)
catches<-data.frame(year=rep(seq(1998,2022),2),
                    data=c(rep('catch',25),rep('survey',25)),
                    value=c(rnorm(length(seq(1998,2022)),20,4),rnorm(length(seq(1998,2022)),20,6)))
cat_plot<-ggplot(catches)+
  geom_line(aes(x=year,y=value))+theme_bw()+
  facet_wrap(~data,ncol=1)+expand_limits(y=0)+ylab('tons')

png('plots/figure_1.png',height=6,width=8,res=350,units='in')
(study_area | p_natl) + plot_layout(widths=c(3.5,1)) +plot_annotation(tag_levels='a')
dev.off()

png('plots/figure_1_alt.png',height=6,width=10,res=350,units='in')
(study_area | cat_plot | p_natl) + plot_layout(widths=c(3.5,1,1)) +plot_annotation(tag_levels='a')
dev.off()