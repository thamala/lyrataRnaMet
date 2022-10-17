#Script for finding differentially methylated genes
#Tuomas Hämälä, October 2022

library(data.table)
library(qvalue)

#generate data for glm, shown here for the cg context
#this was done for all three contexts and the results merged
df <- read.table(gzfile("CG_all.bg.gz"),header=T,na.string=".") #dataframe from bedtools unionbedg (example: bedtools unionbedg -header -filler . -names ind0 ind1 ind2 -i met0.bg met1.bg met2.bg > out.bg)
df <- df[apply(df[,4:ncol(df)],1,function(x){sum(is.na(x))/(ncol(df)-3)}) < 0.5,] #filter out sites with >= 50% missing data
genes <- read.table("lyrata_genes_loc.txt") #locations and names of genes
info <- read.table("lyrata_met_info.txt", header=T) #info file with field and pop for each individual
n <- nrow(genes)
out <- rbindlist(lapply(1:n, function(i){
	if(i==1|i==n|i%%10==0) cat("Line",i,"/",n,"\n")
	gene <- df[df$chrom == genes$V1[i] & df$end >= genes$V2[i] & df$end <= genes$V3[i],]
	if(nrow(gene) == 0) return()
	met <- apply(gene[,4:ncol(gene)], 2, function(x){round(sum(x/100,na.rm=T))})
	unmet <- apply(gene[,4:ncol(gene)], 2, function(x){round(sum(1-x/100,na.rm=T))})
	data.frame(gene=genes$V5[i],cg_met=as.numeric(met),cg_unmet=as.numeric(unmet),cg_sites=nrow(gene),id=names(met))
}))
out <- merge(out, info, by="id")
write.table(out, file="CG_met_glm.txt", row.names=F, quote=F, sep="\t")

#run glm for genes with >= 10 cytosines
df <- read.table("all_contx_met_glm.txt", header=T)
df <- subset(df, (cg_sites + chg_sites + chh_sites) >= 10)
df$ch_prop <- (df$chg_met + df$chh_met) / (df$chg_met + df$chh_met + df$chg_unmet + df$chh_unmet)
df$ch_met <- df$chg_met + df$chh_met
df$ch_unmet <- df$chg_unmet + df$chh_unmet

n <- length(unique(df$gene))
res <- rbindlist(lapply(1:n, function(i){
	tryCatch({
		if(i==1|i==n|i%%100==0) cat("Gene",i,"/",n,"\n")
		gene <- subset(df, gene == unique(df$gene)[i])
		if(sum(gene$cg_met + df$ch_met) < 3) return()
		m0 <- glm(cbind(cg_met, cg_unmet) ~ ch_prop + pop + field + pop:field, family=binomial(link="logit"), data=gene)
		m1 <- glm(cbind(cg_met, cg_unmet) ~ ch_prop + pop + field, family=binomial(link="logit"), data=gene)
		m2 <- glm(cbind(cg_met, cg_unmet) ~ ch_prop + pop, family=binomial(link="logit"), data=gene)
		m3 <- glm(cbind(cg_met, cg_unmet) ~ ch_prop + field, family=binomial(link="logit"), data=gene)
		cg_p_inter <- anova(m0,m1,test="LRT")[2,"Pr(>Chi)"]
		cg_p_field <- anova(m1,m2,test="LRT")[2,"Pr(>Chi)"]
		cg_p_pop <- anova(m1,m3,test="LRT")[2,"Pr(>Chi)"]
		m0 <- glm(cbind(ch_met, ch_unmet) ~ pop + field + pop:field, family=binomial(link="logit"), data=gene)
		m1 <- glm(cbind(ch_met, ch_unmet) ~ pop + field, family=binomial(link="logit"), data=gene)
		m2 <- glm(cbind(ch_met, ch_unmet) ~ pop, family=binomial(link="logit"), data=gene)
		m3 <- glm(cbind(ch_met, ch_unmet) ~ field, family=binomial(link="logit"), data=gene)
		ch_p_inter <- anova(m0,m1,test="LRT")[2,"Pr(>Chi)"]
		ch_p_field <- anova(m1,m2,test="LRT")[2,"Pr(>Chi)"]
		ch_p_pop <- anova(m1,m3,test="LRT")[2,"Pr(>Chi)"]
		data.frame(gene=unique(df$gene)[i],cg_p_inter=cg_p_inter,cg_p_field=cg_p_field,cg_p_pop=cg_p_pop,ch_p_inter=ch_p_inter,ch_p_field=ch_p_field,ch_p_pop=ch_p_pop)
	},warning=function(w){cat("Warning in gene",as.character(unique(df$gene)[i]),":",conditionMessage(w),"\n")})
}))

#add q-values
cg_inter <- na.omit(data.frame(gene=res$gene, p=res$cg_p_inter))
cg_inter$cg_q_inter <- qvalue(cg_inter$p)$qvalues
cg_inter <- cg_inter[,-2]
cg_field <- na.omit(data.frame(gene=res$gene, p=res$cg_p_field))
cg_field$cg_q_field <- qvalue(cg_field$p)$qvalues
cg_field <- cg_field[,-2]
cg_pop <- na.omit(data.frame(gene=res$gene, p=res$cg_p_pop))
cg_pop$cg_q_pop <- qvalue(cg_pop$p)$qvalues
chisq <- qchisq(cg_pop$p,1,lower.tail=F)
lambda <- median(chisq)/qchisq(0.5,1)
newchisq <- chisq/lambda
cg_pop$cg_p_pop_gif <- pchisq(newchisq,1,lower.tail=F)
cg_pop$cg_q_pop_gif <- qvalue(cg_pop$cg_p_pop_gif)$qvalues
cg_pop <- cg_pop[,-2]
res <- merge(res, cg_inter, by="gene", all=T)
res <- merge(res, cg_field, by="gene", all=T)
res <- merge(res, cg_pop, by="gene", all=T)
ch_inter <- na.omit(data.frame(gene=res$gene, p=res$ch_p_inter))
ch_inter$ch_q_inter <- qvalue(ch_inter$p)$qvalues
ch_inter <- ch_inter[,-2]
ch_field <- na.omit(data.frame(gene=res$gene, p=res$ch_p_field))
ch_field$ch_q_field <- qvalue(ch_field$p)$qvalues
ch_field <- ch_field[,-2]
ch_pop <- na.omit(data.frame(gene=res$gene, p=res$ch_p_pop))
ch_pop$ch_q_pop <- qvalue(ch_pop$p)$qvalues
chisq <- qchisq(ch_pop$p,1,lower.tail=F)
lambda <- median(chisq)/qchisq(0.5,1)
newchisq <- chisq/lambda
ch_pop$ch_p_pop_gif <- pchisq(newchisq,1,lower.tail=F)
ch_pop$ch_q_pop_gif <- qvalue(ch_pop$ch_p_pop_gif)$qvalues
ch_pop <- ch_pop[,-2]
res <- merge(res, ch_inter, by="gene", all=T)
res <- merge(res, ch_field, by="gene", all=T)
res <- merge(res, ch_pop, by="gene", all=T)

write.table(res, file="all_contx_glm_res.txt", row.names=F, quote=F, sep="\t")
