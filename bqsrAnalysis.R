#!/usr/bin/env Rscript

#load modules
require(gsalib)

#check BQSR report was passed
if (length(args)!=1) {
  stop("Usage: bqsrAnalysis.R <BQSR_GatkReport>.\n", call.=FALSE)
}

#Load BQSR GatkReport
data = gsa.read.gatkreport(args[1])

#plot expected vs observed quality scores
plot(data$RecalTable1$QualityScore, data$RecalTable1$EmpericalQuality)

#calculate Pearson correlation
pearson = cor.test(data$RecalTable1$QualityScore, data$RecalTable1$EmpericalQuality, method = "pearson", conf.level = 0.95)

#print result
print(pearson$statistic)
print(pearson$p.value)