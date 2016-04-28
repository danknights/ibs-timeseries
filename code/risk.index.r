# Derives an additive risk index: bad bugs minus good bugs
# for y == TRUE
# x is taxon/OTU/feature table
# y is TRUE/FALSE for disease or other phenotype
# alpha determines FDR threshold for inclusion of bugs
# eps is small value to add to all zero entries in feature table
# 
# return value includes risk.index.f, a function that takes a matrix with the same
# named columns as x and returns the risk index for each row.
# diff.tests can be included as a prior result to avoid recomputation
# this should be output from differentiation.test(x,y)
"get.risk.index" <- function(x, y,
		alpha=0.05,
		threshold.method=c('acc','spec')[1],
		threshold.spec=.80,
		transform.type=c('none','asin-sqrt','sqrt')[1],
		threshold=NULL,
		eps=NULL,
		parametric=FALSE, verbose=FALSE, correct.pvalues=FALSE,
		diff.tests=NULL
	){
	if(any(is.na(y))) stop('y must not contain NA values')
	if(any(is.na(x))) stop('x must not contain NA values')
	require('ROCR')

	x <- data.transform(x,transform.type)
	if(is.null(diff.tests)){
		diff.tests <- differentiation.test(x, y, alpha=alpha, parametric=parametric, correct.pvalues=correct.pvalues)
	}

	# if there are a range of alpha values, build multiple lists of hits
	hit.ix <- which(diff.tests$qvalues <= alpha)
	
	if(length(hit.ix) == 0){
		warning('No hits found, predicting risk of 0 for all patients\n')
		good.bugs <- NULL
		bad.bugs <- NULL
		sum.good.bugs <- NA
		sum.bad.bugs <- NA
		risk.index <- rep(0,nrow(x))
		risk.index.f <- function(x) {return (rep(0,nrow(x)))}
		threshold <- NA	
	} else {
		if (verbose) {
		cat("There were",length(hit.ix),"taxa significant at FDR of",alpha,'\n')
		}

		# get rid of zeros (if eps > 0)
		if(!is.null(eps)) x[x==0] <- eps

		# identify good/bad bugs
		good.bugs <- diff.tests$classwise.means[,1] > diff.tests$classwise.means[,2] & diff.tests$qvalues < alpha
		bad.bugs <- diff.tests$classwise.means[,1] < diff.tests$classwise.means[,2] & diff.tests$qvalues < alpha

		sum.good.bugs <- numeric(nrow(x))
		sum.bad.bugs <- numeric(nrow(x))
		if(sum(good.bugs) > 0) sum.good.bugs <- rowSums(x[,good.bugs,drop=F])
		if(sum(bad.bugs) > 0) sum.bad.bugs <- rowSums(x[,bad.bugs,drop=F])

		# risk index is simply additive bad bugs vs good bugs
		risk.index <- sum.bad.bugs - sum.good.bugs

		# risk.index.f 
		"risk.index.f" <- function(x){
			sum.good.bugs.i <- numeric(nrow(x))
			sum.bad.bugs.i <- numeric(nrow(x))
			if(sum(good.bugs) > 0) sum.good.bugs.i <- rowSums(x[,good.bugs,drop=F])
			if(sum(bad.bugs) > 0) sum.bad.bugs.i <- rowSums(x[,bad.bugs,drop=F])
			return(sum.bad.bugs.i - sum.good.bugs.i)
		}

		# choose threshold based on requested threshold method
		if(is.null(threshold)){
			if(threshold.method=='acc'){
				# get threshold that is most accurate for this risk index
				pred <- prediction(risk.index, y)
				perf <- performance(pred, 'acc')
				threshold.ix <- which.max(perf@y.values[[1]])
				threshold <- perf@x.values[[1]][threshold.ix]
			} else if(threshold.method == 'spec'){
				# get threshold that is most accurate for this risk index
				pred <- prediction(risk.index, y)
				perf <- performance(pred, 'sens', 'spec')
				threshold.ix <- max(which(perf@x.values[[1]] >= threshold.spec))
				threshold <- perf@alpha.values[[1]][threshold.ix]
			} else {
				stop(paste('Unknown threshold method:', threshold.method))
			}
		}
	}
			
	res <- list(
		risk.index=risk.index,
		risk.index.f=risk.index.f,
		diff.tests=diff.tests,
		good.bugs=good.bugs,
		bad.bugs=bad.bugs,
		hit.ix=hit.ix,
		alpha=alpha,
		threshold=threshold
	)
	class(res) <- "risk.index"
	return(res)
}

"predict.risk.index" <- function(model, x){
	yhat <- model$risk.index.f(x)
	yhat <- yhat >= model$threshold
	return(yhat)
}

"best.risk.index" <- function(...) {
	res <- tune.risk.index(...)
	return(res$best.model)
}
"tune.risk.index" <- function(x, y, ranges=NULL, 
		alpha=0.05,
		threshold.method=c('acc','spec')[1],
		threshold.spec=.80,
		nfolds=10
	) {
	
	alphas <- alpha
	threshold.specs <- threshold.spec
	if(!is.null(ranges)){
		if('alpha' %in% names(ranges)) alphas <- ranges$alpha
		if('threshold.spec' %in% names(ranges)) threshold.specs <- ranges$threshold.spec
	}
	
	performances <- expand.grid(alpha=alphas, threshold.spec=threshold.specs)
	performances <- cbind(performances, error=rep(0,nrow(performances)), dispersion=rep(0,nrow(performances)))
	
	# performance matrix containing the error rate per fold
	per.fold.errors <- matrix(0, nrow=nrow(performances), ncol=nfolds)
	
	folds <- balanced.folds(y,nfolds=nfolds)
	
	cv.yhat <- matrix(FALSE,nrow=nrow(performances), ncol=nrow(x))
	cv.risk.index <- matrix(0,nrow=nrow(performances), ncol=nrow(x))
	cv.thresholds <- matrix(0,nrow=nrow(performances), ncol=nrow(x))
	cv.models <- list()
	for(i in 1:nrow(performances)) cv.models[[i]] <- list()
	
	for(fold in 1:length(unique(folds))){
		fold.ix <- which(fold == folds)
		# precalculate differentiation.tests
		diff.tests <- differentiation.test(x[-fold.ix,,drop=F], y[-fold.ix], alpha=1, parametric=FALSE, correct.pvalues=FALSE)

		for(i in 1:nrow(performances)){
			alpha <- performances$alpha[i]
			threshold.spec <- performances$threshold.spec[i]
			res <- get.risk.index(x[-fold.ix,,drop=F],y[-fold.ix], alpha=alpha, threshold.spec=threshold.spec, threshold.method=threshold.method, diff.tests=diff.tests)
			cv.risk.index[i,fold.ix] <- res$risk.index.f(x[fold.ix,,drop=F])
			cv.yhat[i,fold.ix] <- cv.risk.index[i,fold.ix] >= res$threshold
			cv.thresholds[i,fold] <- res$threshold
			cv.models[[i]][[fold]] <- res
			per.fold.errors[i, fold] <- mean(y[fold.ix] != cv.yhat[i, fold.ix])
		}
	}
	performances[,'error'] <- rowMeans(per.fold.errors)
	performances[,'dispersion'] <- apply(per.fold.errors,1,sd)
	
	best.ix <- which.min(performances[,'error'])
	best.alpha <- performances[best.ix,'alpha']
	best.threshold.spec <- performances[best.ix,'threshold.spec']
	
	ri <- get.risk.index(x,y,alpha=best.alpha, threshold.spec=best.threshold.spec, threshold.method='spec')	

	res <- list()	
	res$best.model <- ri
	res$best.parameters <- performances[best.ix,1:2]
	res$best.parameters.index <- best.ix
	res$best.performance <- performances[best.ix,3:4]
	res$method <- "get.risk.index"
	res$nparcomb <- nrow(performances)
	res$performances <- performances
	res$sampling <- sprintf('%d-fold cross validation', nfolds)
	res$train.ind <- sapply(1:max(folds), function(xx) which(folds != xx))
	res$cv.thresholds <- cv.thresholds
	res$cv.risk.index <- cv.risk.index
	res$cv.yhat <- cv.yhat
	res$cv.models <- cv.models
	return(res)	
}

# performs cross-validated assessment of risk index performance
# returns hold-out risk indices, chosen thresholds,
# and classification of each sample at its threshold
#
#
# y must be boolean (a logical)
# 
# returns AUC, maximum f1, fprs/tprs for a sequence of thresholds
#
# nfolds=-1 means leave-one-out
"cv.risk.index" <- function(x, y, alpha=.25, nfolds=10, threshold.spec, 
		transform.type='none', correct.pvalues=FALSE, parametric=FALSE){
	x <- data.transform(x, transform.type)
	folds <- balanced.folds(y,nfolds=nfolds)
	if(!is.logical(y)) stop('cv.risk.index: y must be a logical')

	yhat <- y

	risk.index <- numeric(nrow(x))
	thresholds <- numeric(length(unique(folds)))
	models <- list()
	
	for(fold in 1:length(unique(folds))){
		fold.ix <- which(fold == folds)
		res <- get.risk.index(x[-fold.ix,,drop=F],y[-fold.ix], alpha=alpha, correct.pvalues=correct.pvalues, parametric=parametric)
		risk.index[fold.ix] <- res$risk.index.f(x[fold.ix,,drop=F])
		yhat[fold.ix] <- risk.index[fold.ix] >= res$threshold
		thresholds[fold] <- res$threshold
		models[[fold]] <- res
	}
	return(list(risk.index=risk.index, yhat=yhat, thresholds=thresholds, models=models))
}