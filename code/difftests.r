library('beeswarm')
source('code/stats.r')
source('code/util.r')
source('code/risk.index.r')

# source('load.r')
ALPHA <- 0.25

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

        cat(sprintf('Significant q-values for %s, %s:\n',x.name, test.name))
        difftest <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=TRUE)
        difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=FALSE)

        if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
            pdf(sprintf('differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
            for(k in signif.ix){
                if(!is.null(difftest$norm.test.pvals)){
                    norm.test <- difftest$norm.test.pvals[k]
                } else {
                    norm.test <- '0'
                }
                if(norm.test < 0.05){
                    qval <- difftest.np$qvalues[k]
                } else {
                    qval <- difftest$qvalues[k]
                }
                
                cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',norm.test,'\n',sep=''))
                beeswarm(data.transform(test.x)[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),
                # beeswarm(test.x[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),
                        xlab='',ylab='Relative abundance', main=colnames(test.x)[k])
                bxplot(data.transform(test.x)[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
                # bxplot(test.x[test.ix,col.ix] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
            }
            dev.off()
        } else {
            cat("None\n")
        }
    }
}
