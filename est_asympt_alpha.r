#Script for estimating asymptotic alpha
#Based on asymptoticMK by Ben Haller (http://benhaller.com/messerlab/asymptoticMK.html)
#Tuomas Hämälä, January 2022

library(nls2)

#Function for testing the deviation between candidate genes and the genome-wide average
lrt <- function(N,S){
	gw.pN <- gw.N[-c(1,length(gw.N))] / sum(gw.N)
	gw.pS <- gw.S[-c(1,length(gw.S))] / sum(gw.S)
	gw.dN <- (-3/4 * log(1 - (gw.N[length(gw.N)]/sum(gw.N)) * 4/3)) #JC correction for substitutions
	gw.dS <- (-3/4 * log(1 - (gw.S[length(gw.S)]/sum(gw.S)) * 4/3))
	pN <- N[-c(1,length(N))] / sum(N)
	pS <- S[-c(1,length(S))] / sum(S)
	dN <- (-3/4 * log(1 - (N[length(N)]/sum(N)) * 4/3))
	dS <- (-3/4 * log(1 - (S[length(S)]/sum(S)) * 4/3))
	gw.alpha <- 1 - (gw.pN/gw.pS)/(gw.dN/gw.dS)
	test.alpha <- 1 - (pN/pS)/(dN/dS)
	alpha <- c(gw.alpha, test.alpha)
	f <- c(1:(length(N)-2)/(length(N)-1), 1:(length(N)-2)/(length(N)-1))
	x <- c(rep(0,length(N)-2),rep(1,length(N)-2))
	st1 <- expand.grid(a=seq(-1,1,length.out=res+1), b=seq(-1,1,length.out=res), c=seq(1,10,length.out=res+1), d=seq(1,10,length.out=res))
	st2 <- expand.grid(a=seq(-1,1,length.out=res+1), b=seq(-1,1,length.out=res), c=seq(1,10,length.out=res+1))
	fit0 <- nls2(alpha ~ a + b * exp(-c * f) + d * x, start=st1, algorithm="brute-force", control=nls.control(maxiter=nrow(st)))
	fit1 <- nls2(alpha ~ a + b * exp(-c * f) + d * x, start=fit0, control=nls.control(maxiter=200))
	fit2 <- nls2(alpha ~ a + b * exp(-c * f), start=st2, algorithm="brute-force", control=nls.control(maxiter=nrow(st)))
	fit3 <- nls2(alpha ~ a + b * exp(-c * f), start=fit2, control=nls.control(maxiter=200))
	anova(fit1,fit3)[2,"Pr(>F)"]
}

#Function for generating bootstraps
boot <- function(N,S, n=500, verbose=T){
	sel <- rmultinom(n, sum(N), N)
	neut <- rmultinom(n, sum(S), S)
	f <- 1:(length(N)-2)/(length(N)-1)
	res <- sapply(1:n, function(i){
		tryCatch({
			if(verbose == T & (i == 1 | i == n | i %% 10 == 0)) cat("Cycle",i,"/",n,"\n")
			pN <- sel[-c(1,length(N)),i] / sum(sel[,i])
			pS <- neut[-c(1,length(S)),i] / sum(neut[,i])
			dN <- (-3/4 * log(1 - (sel[length(N),i]/sum(sel[,i])) * 4/3))
			dS <- (-3/4 * log(1 - (neut[length(S),i]/sum(neut[,i])) * 4/3))
			alpha <- 1 - (pN/pS)/(dN/dS)
			res <- 10
			st <- expand.grid(a=seq(-1,1,length.out=res+1), b=seq(-1,1,length.out=res), c=seq(1,10,length.out=res+1))
			fit0 <- nls2(alpha ~ a + b * exp(-c * f), start=st, algorithm="brute-force", control=nls.control(maxiter=nrow(st)))
			fit1 <- nls2(alpha ~ a + b * exp(-c * f), start=fit0, control=nls.control(maxiter=200))
			alpha_asympt <- predict(fit1, newdata=data.frame(f=1))
			return(alpha_asympt)
		},error=function(e){
			cat(conditionMessage(e),"\n")
			return(NA)
		})
	})
	return(res)
}

gw.N <- read.table("J1_gw_0fold.sfs",sep=",") #Genome-wide 0fold unfolded sfs
gw.S <- read.table("J1_gw_4fold.sfs",sep=",") #Genome-wide 4fold unfolded sfs
gw.N <- apply(gw.N,2,sum) #If multiple sfs have been concatenated, take a sum of each bin
gw.S <- apply(gw.S,2,sum) 

N <- read.table("J1_deg_field_0fold.sfs",sep=",") #Candidate gene 0fold unfolded sfs
S <- read.table("J1_deg_field_4fold.sfs",sep=",") #Candidate gene 4fold unfolded sfs
N <- apply(N,2,sum)
S <- apply(S,2,sum)

pN <- N[-c(1,length(N))] / sum(N)
pS <- S[-c(1,length(S))] / sum(S)
dN <- (-3/4 * log(1 - (N[length(N)]/sum(N)) * 4/3))
dS <- (-3/4 * log(1 - (S[length(S)]/sum(S)) * 4/3))
alpha <- 1 - (pN/pS)/(dN/dS)
f <- 1:(length(N)-2)/(length(N)-1)
res <- 10
st <- expand.grid(a=seq(-1,1,length.out=res+1), b=seq(-1,1,length.out=res), c=seq(1,10,length.out=res+1))
fit0 <- nls2(alpha ~ a + b * exp(-c * f), start=st, algorithm="brute-force", control=nls.control(maxiter=nrow(st)))
fit1 <- nls2(alpha ~ a + b * exp(-c * f), start=fit0, control=nls.control(maxiter=200))
alpha_asympt <- predict(fit1, newdata=data.frame(f=1))
cat("R2:",1 - sum(resid(fit1)^2)/sum((alpha - mean(alpha))^2), "alpha:", alpha_asympt)
lrt(N,S)
boots <- boot(N,S)
cat(alpha_asympt, quantile(boots, c(0.025,0.975), na.rm=T), sep="\t")
