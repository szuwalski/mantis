#==============================
# plot scenario characteristics
#==============================

scenario_plot<-function(params,outs)
{
trm<-data.frame(Fishing=ceiling(params$fish_prop),
                Growth=params$grow_month,
                Recruitment=params$rec_month,
                month=month.name)
trm$month <- factor(trm$month, levels = trm$month)

tryit<-melt(trm,id.vars='month')
timing_plot<-ggplot(tryit,aes(x=month,y=variable,col=as.factor(value)))+
  geom_point(pch=15,size=10)+theme_bw()+
  theme(legend.position='none')+ylab('')+xlab('')+
  scale_x_discrete(position = "top") + theme(axis.text.x=element_text(angle=45,hjust=0)) +
  scale_color_manual(values=c("light grey","black"))

tlo<-data.frame(Maturity=params$maturity_at_size,
                Molting=params$molt_sel,
                size=params$mid_pts)
tlp<-melt(tlo,id.vars='size')
mat_plot<-ggplot(tlp,aes(x=size,y=value,color=variable))+
  geom_line()+
  theme_bw()+
  theme(legend.position=c(.7,.5),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))+
  ylab("Probability")

tlo<-data.frame(Fishery=params$fish_sel,
                Survey=params$surv_sel,
                size=params$mid_pts)
tlp<-melt(tlo,id.vars='size')
sel_plot<-ggplot(tlp,aes(x=size,y=value,color=variable))+
  geom_line()+
  theme_bw()+
  theme(legend.position=c(.7,.5),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))+
  ylab("Selectivity")

wt_plot<-ggplot(data.frame(Weight=params$wt_at_size,
                           size=params$mid_pts),
                aes(x=size,y=Weight))+
  geom_line()+theme_bw()+
  theme(legend.position='none')

fmort_plot<-ggplot(data.frame(fmort=params$f_mort[,ncol(params$f_mort)],
                              time=params$all_months),
                   aes(x=time,y=fmort))+
  geom_line()+theme_bw()+
  theme(legend.position='none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_continuous(position = "right") +
  ylab("Fishing mortality")


rec_plot<-ggplot(data.frame(recruitment=outs$recruitment,
                            time=params$all_months),
                 aes(x=time,y=recruitment))+
  geom_line()+geom_point()+theme_bw()+
  theme(legend.position='none')+
  scale_y_continuous(position = "right")+
  ylab("Recruitment")



p<-((timing_plot/outs$size_trans[[2]]) | (mat_plot/sel_plot/wt_plot) | (fmort_plot/rec_plot)) + 
  plot_layout(widths=c(2.5,1,1)) +
  plot_annotation(tag_levels='a')

return(p)
}