
source('../functions/man.plot.R')

height = read.csv('../data/height_gwas.txt', sep = '\t', header = T)
furlength = read.csv('../data/furlength_gwas.txt', sep = '\t', header = T)
furnish = read.csv('../data/furnish_gwas.txt', sep = '\t', header = T)

my_genome_gwas = furnish

man_data = man.data.frame(my_genome_gwas)

# cols = c('orange', 'gold', 'lightblue')

png('../plots/furnish.png', height = 600, width = 2000, units = 'px', pointsize = 30)

par(mar = c(3, 1.5, 0, 1))

man_plot = man.plot(man_data, threshold = 8.46)
man.highlight(man_plot, threshold = 8.46)
# man.label(man_plot)

dev.off()





-log10(5e-8)


qqnorm(man_data$transformed_pvalue, frame = F)
qqline(man_data$transformed_pvalue, col = 'red', lwd = 2)

man_data$transformed_pvalue[man_data$transformed_pvalue > 0.8 * max(man_data$transformed_pvalue)]

quantile(man_data$transformed_pvalue, probs = seq(0, 1, 0.1))['90%']

quantile(man_data$transformed_pvalue, probs = seq(0.9, 1))

