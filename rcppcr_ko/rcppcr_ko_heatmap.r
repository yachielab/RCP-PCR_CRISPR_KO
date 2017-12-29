#!/usr/bin/env Rscript
library(ggplot2)
library(scales) # for muted function
library(cowplot)
library(gridExtra)
library(grid)

args = commandArgs(trailingOnly=TRUE)

group.colors <- c("-" = "#dbdada", "gRNA" = "#80c0e5", "PAM" ="#f79999", "Del" = "#000000", "Ins" = "#7f0114","Mut"="#f5f99f","Ref"="#FFFFFF" ,"- (WT)" = "#dbdada","+ (Indel)"= "#f5f99f","++ (Frameshift)"="#FFD700")

b = read.csv(file=args[1],header=TRUE)
h = read.csv(file=args[2],header=TRUE)

H = h[ h$base_index == 0, ]
h = h[ h$base_index > 0, ]


t = ggplot(H, aes(base_index,ORDER)) +   geom_tile(aes(fill = factor(H$Char_stat)))+geom_text(data=H,aes(y=ORDER,label=SampleID))  + theme_bw()  + theme(plot.margin= unit(c(1, 1, -1, 1), "lines")) + scale_x_continuous(expand=c(0,0),labels=NULL,breaks=c(),position="top" ) +  scale_y_continuous(expand=c(0,0),labels = NULL,breaks=c()) +scale_fill_manual(values=group.colors) + theme(legend.position='none')+ theme(axis.title.y = element_blank())+labs(x= "KO Status")
p = ggplot(b, aes(x=ORDER, y=Ratio)) + geom_bar(stat="identity",  position=position_dodge()) + geom_errorbar(aes(ymin=Ratio-SE, ymax=Ratio+SE), width=.2, position=position_dodge(2.5))+theme_classic() + theme(axis.title.y = element_blank(),axis.text.x =element_text(size=12,color = "#000000"))+ scale_x_continuous(expand=c(0,0),labels=NULL)  +  scale_y_continuous(expand=c(0,0),position="top",breaks=c(0,0.5,1))+coord_flip()+labs(y= "Relative read out within well")+  geom_hline(aes(yintercept=0.25), colour="#D3D3D3", linetype="dashed",alpha=0.3)+  geom_hline(aes(yintercept=0.5), colour="#D3D3D3", linetype="dashed",alpha=0.3)+  geom_hline(aes(yintercept=0.75), colour="#D3D3D3", linetype="dashed",alpha=0.3)+  geom_hline(aes(yintercept=1.0), colour="#D3D3D3", linetype="dashed",alpha=0.3)
T = ggplot(h, aes(base_index,ORDER)) +   geom_tile(aes(fill = factor(h$Char_stat)))  + theme_bw()   + scale_x_continuous(expand=c(0,0),position="top",breaks=c(25,50,75))  +  scale_y_continuous(expand=c(0,0),labels = NULL,breaks=c())+scale_fill_manual(values=group.colors) + theme(legend.position='none')+ theme(axis.title.y = element_blank(),axis.text.x =element_text(size=12,color = "#000000"))+labs(x = "Target locus position (bp)")


g1 = t+ theme(plot.margin = unit(c(1, 0, 1, 1), "cm"))
g2 = T+ theme(plot.margin = unit(c(1, 0, 1, 0), "cm"))
g3 = p+ theme(plot.margin = unit(c(1, 1, 1, 0), "cm"))

#Arrange them in a grid
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg3 <- ggplot_gtable(ggplot_build(g3))

P = plot_grid(gg1, gg2, gg3, align = "h", ncol = 3, rel_widths = c(2.5/10,5/10,2.5/10))

pH =  nrow(H) *0.2
pW = 10

ggsave(args[3],width=pW,height=pH,limitsize=FALSE)
