library(nortest)
library(ggplot2)
library(forcats)

#Read expression data
intragenic_data = read.table("GTEx_v8_RNASeQCv1.1.9_rtc_intragenic_median_tpm.txt", header=T, row.names=1)
intergenic_data = read.table("GTEx_v8_RNASeQCv1.1.9_rtc_intergenic_median_tpm.txt", header=T, row.names=1)
 
#Read list of tissues
tissues = read.table("Tissues.txt", header=F)
colnames(tissues) = c("tissue")

intragenic_data = intragenic_data[match(tissues$tissue,names(intragenic_data))]
intergenic_data = intergenic_data[match(tissues$tissue,names(intergenic_data))]

intragenic_data[intragenic_data == 0] <- NA
intergenic_data[intergenic_data == 0] <- NA

log_intragenic_data = log(intragenic_data+1,2)
log_intergenic_data = log(intergenic_data+1,2)

#Create empty dataframes to store results
df_normality = data.frame(tissue=character(), 
	type=character(), 
	shapiro.p.value=numeric(),
	logshapiro.p.value=numeric(),
	ad.p.value=numeric(),
	logad.p.value=numeric(),
	length=numeric())

# normality tests for every tissue distribution
for(tissue in tissues$tissue){
	intra_tmp = na.omit(intragenic_data[[tissue]])
	inter_tmp = na.omit(intergenic_data[[tissue]])

	logintra_tmp = na.omit(log_intragenic_data[[tissue]])
	loginter_tmp = na.omit(log_intergenic_data[[tissue]])

	#Add result to dataframe
	df_normality = rbind(df_normality, data.frame(tissue=tissue, 
		type=as.character("intragenic"), 
		shapiro.p.value=shapiro.test(intra_tmp)$p.value,
		logshapiro.p.value=shapiro.test(logintra_tmp)$p.value,
		ad.p.value=ad.test(intra_tmp)$p.value,
		logad.p.value=ad.test(logintra_tmp)$p.value,
		length=length(intra_tmp)))

	df_normality = rbind(df_normality, data.frame(tissue=tissue, 
		type=as.character("intergenic"), 
		shapiro.p.value=shapiro.test(inter_tmp)$p.value,
		logshapiro.p.value=shapiro.test(loginter_tmp)$p.value,
		ad.p.value=ad.test(inter_tmp)$p.value,
		logad.p.value=ad.test(loginter_tmp)$p.value,
		length=length(inter_tmp)))

}

df_normality$shapiro.p.adjust = p.adjust(df_normality$shapiro.p.value, method="bonferroni", n=length(df_normality$shapiro.p.value))
df_normality$logshapiro.p.adjust = p.adjust(df_normality$logshapiro.p.value, method="bonferroni", n=length(df_normality$logshapiro.p.value))
df_normality$ad.p.adjust = p.adjust(df_normality$ad.p.value, method="bonferroni", n=length(df_normality$ad.p.value))
df_normality$logad.p.adjust = p.adjust(df_normality$logad.p.value, method="bonferroni", n=length(df_normality$logad.p.value))

# write to output
write.table(df_normality, file = "normality_test.txt", sep = "\t",
            row.names = FALSE)

print("intra kruskal p value")
print(kruskal.test(intragenic_data)$p.value)

print("logintra kruskal p value")
print(kruskal.test(log_intragenic_data)$p.value)

print("inter kruskal p value")
print(kruskal.test(intergenic_data)$p.value)

print("log inter kruskal p value")
print(kruskal.test(log_intergenic_data)$p.value)

wilcox <- data.frame(tissue=character(), wilcox.p.value=numeric())

for(tissue in tissues$tissue){
	intra_tmp = na.omit(intragenic_data[[tissue]])
	inter_tmp = na.omit(intergenic_data[[tissue]])
  
	wilcox = rbind(wilcox, data.frame(tissue=tissue,
		wilcox.p.value=wilcox.test(inter_tmp,intra_tmp)$p.value))
}

wilcox$p.value.adjust = p.adjust(wilcox$wilcox.p.value, method = 'bonferroni')
wilcox = wilcox[order(wilcox$p.value.adjust),]

# write to output
write.table(wilcox, file = "wilcox.txt", sep = "\t",
            row.names = FALSE)

wilcox_sign = wilcox[wilcox$p.value.adjust <= 0.05,]
print(wilcox_sign)

df_signif <- data.frame(expression= numeric(), type=character(), tissue=character())
for(tissue in wilcox_sign$tissue){
	df_intra = data.frame(na.omit(log_intragenic_data[[tissue]]), "intragenic", tissue)
	colnames(df_intra) = c("expression", "type", "tissue")
	df_inter = data.frame(na.omit(log_intergenic_data[[tissue]]), "intergenic", tissue)
	colnames(df_inter) = c("expression", "type", "tissue")

	df_signif = rbind(df_signif, df_intra)
	df_signif = rbind(df_signif, df_inter)
}

signif_plot = ggplot(df_signif, aes(x=tissue, y=expression, fill=type)) + 
	geom_boxplot(outlier.shape = NA) +
	scale_y_continuous(limits = quantile(df_signif$expression, c(0.0, 0.9))) +
	labs(title="", y="log2(TPM+1)", x="") +
	theme_classic() +
	theme(text = element_text(size=15), legend.position="top") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
	scale_x_discrete(limits=wilcox_sign$tissue)

svg(filename = "signif.svg",
	width=8.3,
	height=5.8)
signif_plot
dev.off()

warnings()


