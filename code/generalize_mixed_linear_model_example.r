library('lme4')

# generate log-normal microbiome data; 100 samples with 10 species
x <- matrix(exp(rnorm(1000)),100); colnames(x) <- sprintf('Species%03d',1:ncol(x)); rownames(x) <- sprintf('Sample%03d',1:nrow(x)); x <- sweep(x,1,rowSums(x),'/')

# fake metadata with symptoms coming and going within patients
# 10 samples per patient, with 10 patients
map <- data.frame(symptoms=sample(rep(c(0,1),each=50)), subject=factor(rep(sprintf('Subject%03d',1:10),each=10)))

# generalized linear mixed models -- test only the first species
glm1 <- glmer(symptoms ~ x[,1] + (1 | subject), family=binomial, data=map)

# print results summary
summary(glm1)

# extract p-value
summary(glm1)[[10]][2,]