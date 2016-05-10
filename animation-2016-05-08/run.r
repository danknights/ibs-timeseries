source('../code/load.r')
source('../code/animate.r')
library('vegan')
library('car')

bd <- vegdist(xsp)

visit.counts <- sapply(split(m$Visit_no., m$Patient_no.),length)
multi.ix <- m$Patient_no. %in% names(visit.counts)[visit.counts > 1]

pc <- cmdscale(as.matrix(bd)[multi.ix,multi.ix],3)
m.multi <- droplevels(m[multi.ix,])


for(step in seq(0,max(m$Visit_no.)+3,.025)){
    cat(step,'\n')
    png(sprintf('images/step-%0.5f.png',step),width=800,height=800)
    plot.animation.frame(pc, m.multi$Patient_no., m.multi$Visit_no., 
                         step=step, c('#FF0000','#00FF00','#0000FF')[m.multi$Cohort],
                         ellipse.alpha=.1)
    dev.off()
}
