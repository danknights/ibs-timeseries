library('vegan')
ALPHA <- .05

dsp <- as.matrix(vegdist(xsp))
dgn <- as.matrix(vegdist(xgn))

# get all within-patient distances
patient.nos <- sort(unique(m$Patient_no.))
wds.sp <- numeric(length(patient.nos))
names(wds.sp) <- sprintf('Subject_%03d',patient.nos)
wds.gn <- numeric(length(patient.nos))
names(wds.gn) <- sprintf('Subject_%03d',patient.nos)
for(i in 1:length(unique(m$Patient_no.))){
    patient.no <- patient.nos[i]
    # within-patient distances
    patient.ix <- m$Patient_no. == patient.no
    if(sum(patient.ix) == 1){
        wds.sp[i] <- NA
        wds.gn[i] <- NA
    } else {
        # calculates all within-cloud distances
        tmp <- dsp[patient.ix,patient.ix]
        tmp <- tmp[upper.tri(tmp)]

        # calculates only distance to previous timepoint
        # tmp <- sapply(1:(sum(patient.ix)-1), function(ixx) dsp[which(patient.ix)[ixx],which(patient.ix)[ixx+1]])
        wds.sp[i] <- mean(tmp)

        # calculates all within-cloud distances
        tmp <- dgn[patient.ix,patient.ix]
        tmp <- tmp[upper.tri(tmp)]
        # calculates only distance to previous timepoint
        # tmp <- sapply(1:(sum(patient.ix)-1), function(ixx) dgn[which(patient.ix)[ixx],which(patient.ix)[ixx+1]])
        wds.gn[i] <- mean(tmp)
    }
}

# test species and genus distances
test.xs <- list(species=wds.sp, genus=wds.gn)
# different group comparisons
test.ixs <- list('HC v. IBS'=ixc.hc | ixc.ibs,
                 'HC v. IBSD'=ixc.hc | ixc.ibsd,
                 'HC v. IBSC'=ixc.hc | ixc.ibsc,
                 'IBSC v. IBSD'=ixc.ibsc | ixc.ibsd
                )
# list of "reference" group names
compare.to <- list('HC v. IBS'='Healthy',
                 'HC v. IBSD'='Healthy',
                 'HC v. IBSC'='Healthy',
                 'IBSC v. IBSD'='IBS-C'
                )

# run all combinations of tests
for(i in 1:length(test.xs)){
    x.name <- names(test.xs)[i]
    test.x <- test.xs[[i]]

    for(j in 1:length(test.ixs)){
        test.name <- names(test.ixs)[j]
        compare.to.name <- compare.to[[test.name]]
        test.ix <- test.ixs[[j]] & !is.na(test.x)

        tt <- t.test(test.x[test.ix] ~ mc$Cohort[test.ix] != compare.to.name) 

        if(tt$p.value < ALPHA) {
            cat(sprintf('t-test variability, %s, %s: ', x.name, test.name))
            cat('p=',round(tt$p.value,4), ', t statistic=',round(tt$statistic,4),'\n',sep='')

            pdf(sprintf('variability_%s_%s.pdf',x.name, gsub(' ','_',test.name)),width=4,height=4)
            beeswarm(test.x[test.ix] ~ mc$Cohort[test.ix] != compare.to.name,
                     xlab=paste(test.name, ' (FALSE is ', compare.to.name,')',sep=''),
                     ylab='Variability (mean within-subject distance)')
            bxplot(test.x[test.ix] ~ mc$Cohort[test.ix] != compare.to.name,add=TRUE)
            dev.off()
        }
    }
}
