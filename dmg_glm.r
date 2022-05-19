#Script for finding differentially methylated genes
#Tuomas Hämälä, January 2022

library(data.table)
library(qvalue)

df <- read.table(gzfile("CG_all.bg.gz"),header=T,na.string=".") #dataframe from bedtools unionbedg (example: bedtools unionbedg -header -filler . -names ind0 ind1 ind2 -i met0.bg met1.bg met2.bg > out.bg)
df <- df[apply(df[,4:ncol(df)],1,function(x){sum(is.na(x))/(ncol(df)-3)}) < 0.5,] #filter out sites with >= 50% missing data
genes <- read.table("lyrata_genes_loc.txt") #locations and names of genes
info <- read.table("lyrata_met_info.txt", header=T) #info file with field and pop for each individual

#generate data for glm
n <- nrow(genes)
out <- rbindlist(lapply(1:n, function(i){
	if(i==1|i==n|i%%10==0) cat("Line",i,"/",n,"\n")
	gene <- df[df$chrom == genes$V1[i] & df$end >= genes$V2[i] & df$end <= genes$V3[i],]
	if(nrow(gene) == 0) return()
	met <- apply(gene[,4:ncol(gene)], 2, function(x){round(2*sum(x/100,na.rm=T))})
	unmet <- apply(gene[,4:ncol(gene)], 2, function(x){round(2*sum(1-x/100,na.rm=T))})
	data.frame(gene=genes$V5[i],met=as.numeric(met),unmet=as.numeric(unmet),sites=nrow(gene),id=names(met))
}))
out <- merge(out, info, by="id")
write.table(out, file="CG_met_glm.txt", row.names=F, quote=F, sep="\t")

#run glm for genes with >= 10 cytosines
df <- read.table("CG_met_glm.txt", header=T)
df <- df[df$sites >= 10,]

n <- length(unique(df$gene))
res <- rbindlist(lapply(1:n, function(i){
	tryCatch({
		if(i==1|i==n|i%%100==0) cat("Gene",i,"/",n,"\n")
		gene <- df[df$gene == unique(df$gene)[i],]
		if(sum(gene$met) < 3) return()
		m0 <- glm(cbind(met, unmet) ~ pop + field + pop:field, family=binomial(link="logit"), data=gene)
		m1 <- glm(cbind(met, unmet) ~ pop + field, family=binomial(link="logit"), data=gene)
		m2 <- glm(cbind(met, unmet) ~ pop, family=binomial(link="logit"), data=gene)
		m3 <- glm(cbind(met, unmet) ~ field, family=binomial(link="logit"), data=gene)
		p_inter <- anova(m0,m1,test="LRT")[2,"Pr(>Chi)"]
		p_field <- anova(m1,m2,test="LRT")[2,"Pr(>Chi)"]
		p_pop <- anova(m1,m3,test="LRT")[2,"Pr(>Chi)"]
		low <- subset(gene, field=="LomSite")
		high <- subset(gene, field=="SpSite")
		lfc <- log2((sum(low$met)/sum(c(low$met,low$unmet)))/((sum(high$met)/sum(c(high$met,high$unmet)))))
		data.frame(gene=unique(df$gene)[i],p_inter=p_inter,p_field=p_field,p_pop=p_pop, sites=gene$sites[1], mean=sum(gene$met)/sum(c(gene$met,gene$unmet)), lfc=lfc)
	},warning=function(w){cat("Warning in gene",as.character(unique(df$gene)[i]),":",conditionMessage(w),"\n")})
}))
inter <- na.omit(data.frame(gene=res$gene, p=res$p_inter))
inter$q_inter <- qvalue(inter$p)$qvalues
inter <- inter[,-2]
field <- na.omit(data.frame(gene=res$gene, p=res$p_field))
field$q_field <- qvalue(field$p)$qvalues
field <- field[,-2]
pop <- na.omit(data.frame(gene=res$gene, p=res$p_pop))
pop$q_pop <- qvalue(pop$p)$qvalues
chisq <- qchisq(pop$p,1,lower.tail=F)
lambda <- median(chisq)/qchisq(0.5,1)
pop$p_pop_gif <- pchisq(chisq/lambda,1,lower.tail=F)
pop$q_pop_gif <- qvalue(pop$p_pop_gif)$qvalues
pop <- pop[,-2]
res <- merge(res, inter, by="gene", all=T)
res <- merge(res, field, by="gene", all=T)
res <- merge(res, pop, by="gene", all=T)

write.table(res, file="CG_glm_res.txt", row.names=F, quote=F, sep="\t")
