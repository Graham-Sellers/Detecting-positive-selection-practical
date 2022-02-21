
source('../functions/man.plot2.R')

height = read.csv('../data/height_gwas.txt', sep = '\t', header = T)
furlength = read.csv('../data/furlength_gwas.txt', sep = '\t', header = T)
furnish = read.csv('../data/furnish_gwas.txt', sep = '\t', header = T)

test = read.csv('../data/test_gwas.txt', sep = '\t', header = T)

png('../plots/test.png', height = 600, width = 2000, units = 'px', pointsize = 30)
par(mar = c(3, 1.5, 0, 1))
test_plot = man.plot(test)
dev.off()

png('../plots/furlength.png', height = 600, width = 2000, units = 'px', pointsize = 30)
par(mar = c(3, 1.5, 0, 1))
furlength_plot = man.plot(furlength, threshold = 8.46)
dev.off()

png('../plots/furnish2.png', height = 600, width = 2000, units = 'px', pointsize = 30)
par(mar = c(3, 1.5, 0, 1))
furnish_plot = man.plot(furnish, threshold = 8.46)
dev.off()

png('../plots/height.png', height = 600, width = 2000, units = 'px', pointsize = 30)
par(mar = c(3, 1.5, 0, 1))
height_plot = man.plot(height, threshold = 8.46)
dev.off()


unique(height_plot$gene)
