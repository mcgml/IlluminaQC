#!/usr/bin/env Rscript

#load modules
require(gsalib)
require(optparse)

#parse command-line args
option_list = list(
  make_option(c("-r", "--report"), action="store", default='', type="character", help="Path to GATK BQSR report")
)
opt=parse_args(OptionParser(option_list=option_list))

#Load BQSR GatkReport
data = gsa.read.gatkreport(opt$report)

#plot expected vs observed quality scores
plot(data$RecalTable1$QualityScore, data$RecalTable1$EmpiricalQuality)

#calculate Pearson correlation
pearson = cor.test(data$RecalTable1$QualityScore, data$RecalTable1$EmpiricalQuality, method = "pearson", conf.level = 0.95)

#print result
print(pearson$estimate, digits = 3)
print(pearson$p.value, digits = 3)