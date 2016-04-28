ix <- ixc.hc | ixc.ibsd; disease.case <- "IBS-D"
ri <- cv.risk.index(sqrt(xgnc[ix,]), droplevels(mc$Cohort[ix]) == disease.case, 
                    alpha=.25, correct.pvalues=FALSE, nfolds=-1)
cat('t-test: hold-out genus-level risk index IBS-D vs. healthy:\n')
cat('p=',t.test(ri$risk.index ~ droplevels(mc$Cohort[ix]) != "Healthy")$p.value, '\n', sep='')

pdf('risk_index_genus_IBSD_v_Healthy.pdf',width=4,height=4)
beeswarm(ri$risk.index ~ droplevels(mc$Cohort[ix]) == disease.case,
        xlab='IBSD', ylab='Bacterial Risk Index')
bxplot(ri$risk.index ~ droplevels(mc$Cohort[ix]) == disease.case,add=TRUE)
dev.off()

ix <- ixc.hc | ixc.ibs
ri <- cv.risk.index(sqrt(xgnc[ix,]), droplevels(mc$Cohort[ix]) != "Healthy", 
                    alpha=.25, correct.pvalues=FALSE, nfolds=-1)
cat('t-test: hold-out genus-level risk index IBS vs. healthy:\n')
cat('p=',t.test(ri$risk.index ~ droplevels(mc$Cohort[ix]) != "Healthy")$p.value, '\n', sep='')
