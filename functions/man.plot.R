
# Here we create the "man.data.frame()" function.
# It converts the input file to a meaningful by-chromosome base pair position.
# and also adds in the -log10 transformed pvalue (nice!)

# manhattan dataframe maker:
man.data.frame = function(df) {
  
  full_df = as.data.frame(df[FALSE,])
  
  for(chromosome in unique(df$chrom)){
    chromosome_df = df[df$chrom == chromosome,]
    chromosome_df = chromosome_df[order(chromosome_df$bp),]
    chromosome_df$relative_bp = c(1, diff(chromosome_df$bp))
    chromosome_df$transformed_pvalue = -log10(chromosome_df$pvalue)
    full_df = rbind(full_df, chromosome_df)
  }
  
  full_df = full_df[,c('chrom', 'relative_bp', 'transformed_pvalue', 'gene')]
  return(full_df)
}

# manhattan plotter:
man.plot = function(df, chroms, threshold, point.cols, line.col) {
  
  # deal with missing inputs and set defaults:
  if(missing(chroms)){
    df = df
  } else {
      df = df[df$chrom %in% chroms,]}
  
  if(missing(threshold)) {
    significance = -log10(5e-8)
  } else {
    significance = threshold}
  
  if(missing(point.cols)) {
    cols = rep(c('grey75', 'grey50'), length(unique(df$chrom)))[1:length(unique(df$chrom))]
    names(cols) = unique(df$chrom)
  } else {
      cols = rep(point.cols, length(unique(df$chrom)))[1:length(unique(df$chrom))]
      names(cols) = unique(df$chrom)}
  
  if(missing(line.col)) {
    line.col = 'orange'
  } else {
    line.col = line.col[1]}
  
  # gap between chromosome:
  space = max(cumsum(df$relative_bp))/length(unique(df$chrom)) * 0.25
  
  # set the colours for the plot:
  df$cols =  cols[as.character(df$chrom)]
  #df$cols = ifelse(df$transformed_pvalue >=8, 'red', cols[as.character(df$chrom)])
  
  # set plot points for each SNP:
  for(c in unique(df$chrom)){
    df[df$chrom == c,][1,]$relative_bp = space}
  
  df$points = cumsum(df$relative_bp)
  
  # set x axis label points for chromosomes:
  chrom_labs = c()
  for(c in unique(df$chrom)){
    daf = df[df$chrom == c,]
    chrom_labs = c(chrom_labs, mean(range(daf$points)))}
  
  # set the upper limit for y:
  y_max = max(df$transformed_pvalue) * 1.2
  x_max = max(cumsum(df$relative_bp))
  
  # make the basic plot with no axes:
  plot(1, 1,
       ylim = c(0, y_max), xlim = c(0, x_max),
       cex = 0, axes = F,
       col = df$cols,
       ylab = NA,
       xlab = NA)
  
  #add the points:
  points(df$points, df$transformed_pvalue,
    col = df$cols,
    pch = 16, cex = 0.4)
  
  # add significance threshold dotted line
  segments(0, significance, max(df$points) + space, significance, lty = 5, lwd = 2, col = line.col)

  # add axes:
  axis(1, at = chrom_labs, labels = unique(df$chrom), lwd = 0, lwd.ticks = 1.5,
       las = 1, pos = -1, tck = -0.01, cex.axis = 0.6)
  axis(2, at = c(seq(0, y_max, 5)), labels = c(seq(0, y_max, 5)), lwd = 0, lwd.tick = 1.5,
       las = 1, pos = 0, tck = -0.01, cex.axis = 0.6, xpd = T)
  segments(0, -1, max(df$points) + space, -1, xpd = T)
  segments(0, -1, 0, y_max, xpd = T)
  
  mtext(expression(paste('-log10(', italic('P'), '_wald)')), side = 2, cex = 1.2)
  
  return(df)
}

man.highlight = function(df, highlight, threshold){
  
  if(missing(threshold)) {
    significance = -log10(5e-8)
  } else {
    significance = threshold}
  
  if(missing(highlight)){
    cols = 'red'
  } else {
    cols = highlight}
  
  high_df = df[df$transformed_pvalue >= significance,]
  points(high_df$points, high_df$transformed_pvalue,
         col = cols,
         pch = 16, cex = 0.4)
}


man.label = function(df, threshold) {
  
  if(missing(threshold)) {
    significance = -log10(5e-8)
  } else {
    significance = threshold}
  
  high_df = df[df$transformed_pvalue >= significance,]
  maxpoints = high_df[FALSE,]
  
  for(gene in unique(high_df$gene)){
    hip = high_df[high_df$gene == gene,]
    max_point = hip[hip$transformed_pvalue == max(hip$transformed_pvalue),]
    maxpoints = rbind(maxpoints, max_point[1,])
    text(max_point[1,]$points, max_point[1,]$transformed_pvalue, gene, pos=4, cex = 0.6, xpd = T)
  }
  
  return(maxpoints)
}

#df = man_data

#high_df = df[df$transformed_pvalue >= 8.64]
#unique(high_df$gene,)

