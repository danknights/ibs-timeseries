library('RColorBrewer')
library('vegan')
ALPHA <- 0.05

# from https://chitchatr.wordpress.com/2011/12/30/convex-hull-around-scatter-plot-in-r/
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  

# test species and genus differences
test.xs <- list(species=xspc, genus=xgnc)
# different group comparisons
test.ixs <- list('HC v. IBS'=ixc.hc | ixc.ibs,
                 'HC v. IBSD'=ixc.hc | ixc.ibsd,
                 'HC v. IBSC'=ixc.hc | ixc.ibsc,
                 'IBSC v. IBSD'=ixc.ibsc | ixc.ibsd
                )



# run all combinations of tests
for(i in 1:length(test.xs)){
    x.name <- names(test.xs)[i]
    test.x <- test.xs[[i]]

    for(j in 1:length(test.ixs)){
        test.name <- names(test.ixs)[j]
        test.ix <- test.ixs[[j]]

        # ANOSIM test for clustering
        aa <- anosim(as.dist(as.matrix(vegdist(test.x))[test.ix,test.ix]), droplevels(mc$Cohort[test.ix]))
        if(aa$signif <= ALPHA){
            cat(sprintf('Significant ANOSIM clustering for %s, %s: ',x.name, test.name))
            cat(paste('R=',aa$statistic,', p=', aa$signif,'\n',sep=''))
        }

        # PERMANOVA test for clustering
        aa <- adonis(as.dist(as.matrix(vegdist(test.x))[test.ix,test.ix]) ~ droplevels(mc$Cohort[test.ix]))
        pval <- aa$aov.tab[1,'Pr(>F)']
        if(pval <= ALPHA){
            cat(sprintf('Significant ADONIS clustering for %s, %s: ',x.name, test.name))
            cat(paste('Partial R2=',round(aa$aov.tab[1,'R2'],3),', p=', pval,'\n',sep=''))
            pc <- cmdscale(vegdist(test.x))

            pdf(sprintf('clustering_%s_%s.pdf',x.name, gsub(' ','_',test.name)),width=4,height=4)
            # plot(pc[,1], pc[,2], col=brewer.pal(9,'Set1')[-6][as.numeric(droplevels(mc$Cohort[test.ix]))])
            # for(lev in as.numeric(droplevels(mc$Cohort[test.ix]))){
            #     ix <- as.numeric(droplevels(mc$Cohort[test.ix])) == lev
            #     Plot_ConvexHull(pc[ix,1], pc[ix,2], lcolor=brewer.pal(9,'Set1')[-6][lev])
            # }

            plot(pc[,1], pc[,2], col=brewer.pal(9,'Set1')[-6][as.numeric(mc$Cohort[test.ix])])
            for(lev in as.numeric(mc$Cohort[test.ix])){
                ix <- as.numeric(mc$Cohort[test.ix]) == lev
                Plot_ConvexHull(pc[ix,1], pc[ix,2], lcolor=brewer.pal(9,'Set1')[-6][lev])
            }
            dev.off()
        }
    }
}