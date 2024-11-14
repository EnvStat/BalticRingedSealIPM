
error_bars <- function(quantiles, x_data=NULL, y_data=NULL, x_ticks=NULL, alpha=0.5, shape=1, cex=1, main='', xlab='Time', ylab='Value') {
  if(is.null(x_ticks)) {x_ticks = 1:nrow(quantiles)}
  if(is.null(x_data)) {x_data = x_ticks}
  if(is.null(y_data)) {y_data = rep(NaN, length(x_ticks))}
  df.quantiles <- data.frame(time=x_ticks, lb=quantiles[,1], mid=quantiles[,2], ub=quantiles[,3])
  df.data <- data.frame(time=x_data, values=y_data)

  ggplot() +
    geom_errorbar(data=df.quantiles, aes(x=x_ticks, ymin = lb, ymax = ub), width = 0.5, color='black', alpha=alpha) +  # Add error bars (quantiles 25% to 75%)
    geom_point(data=df.quantiles, aes(x = x_ticks, y = mid), color='black', size=0.8, alpha=alpha) +  
    geom_point(data=df.data, aes(x=time, y=values), size=cex, shape=shape) + # actual data
    labs(title = main, x = xlab, y = ylab) +
    theme_classic()
}

ribbon_plot <- function(quantiles, x_data=NULL, y_data=NULL, x_ticks=NULL, ylim=NULL, reference = NULL, asymptotic = NULL,
                        place_legend = 'none', legend_title = 'Group', color_pal='Set1', lwd=1, main='', xlab='Time', ylab='Value') {
  if(is.null(x_ticks)) {x_ticks = 1:nrow(quantiles)}
  if(is.null(x_data)) {x_data = x_ticks}
  if(is.null(y_data)) {y_data = rep(NaN, length(x_ticks))}
  if(!is.list(quantiles)) {
    if(is.null(ylim)) {ylim = c(0, max(quantiles, na.rm = T))}
    df.quantiles <- data.frame(time=x_ticks, lb=quantiles[,1], mid=quantiles[,2], ub=quantiles[,3], group=rep(1, length(x_ticks)))
  } else {
    if(is.null(ylim)) {ylim = c(0, max(sapply(quantiles, max, na.rm = T)))}
    df.quantiles <- data.frame(time=x_ticks, lb=quantiles[[1]][,1], mid=quantiles[[1]][,2], ub=quantiles[[1]][,3], group=1)
    for(i in 2:length(quantiles)) {
      df.quantiles.add <- data.frame(time=x_ticks, lb=quantiles[[i]][,1], mid=quantiles[[i]][,2], ub=quantiles[[i]][,3], group=i)
      df.quantiles <- rbind(df.quantiles, df.quantiles.add)
    }
  }
  df.data <- data.frame(time=x_data, values=y_data)
  
  plot.object <- ggplot() +
    geom_line(data=df.quantiles, aes(x = time, y = mid, group=as.factor(group), color=as.factor(group)), linewidth=lwd) +  
    geom_ribbon(data=df.quantiles, aes(x = time, ymin = lb, ymax = ub, group=as.factor(group), fill=as.factor(group)), alpha=0.2)
  
  if(!is.null(asymptotic)) {
    df.asymptotic <- data.frame(time=c(max(x_ticks)+1, max(x_ticks)+3), 
                                lb=rep(asymptotic[1],2), mid=rep(asymptotic[2],2), ub=rep(asymptotic[3],2), group=rep(1, 2))
    plot.object <- plot.object + 
      geom_line(data=df.asymptotic, aes(x = time, y = mid, group=as.factor(group), color=as.factor(group)), linewidth=lwd) +  
      geom_ribbon(data=df.asymptotic, aes(x = time, ymin = lb, ymax = ub, group=as.factor(group), fill=as.factor(group)), alpha=0.2)
  }
  
  plot.object <- plot.object +
    geom_point(data=df.data, aes(x=time, y=values)) + # actual data
    geom_hline(yintercept = reference, linetype = "dashed") + 
    scale_y_continuous(limits = c(ylim[1], ylim[2])) +
    labs(title = main, x = xlab, y = ylab, color = legend_title, fill = legend_title) +
    theme_classic() +
    theme(legend.position = place_legend, text = element_text(size=25))
  
  if(is.list(quantiles)) {
    plot.object <- plot.object + scale_fill_brewer(palette=color_pal)
  } else {
    plot.object <- plot.object + scale_fill_grey() + scale_color_grey()
  }
  
  plot.object
}


pp_raster <- function(pp_samples, pp_x, bin_width, x_data=NULL, y_data=NULL, exclusions=NULL, x_lab='X', y_lab='Y', cex=1, interpolate=T) {
  breaks <- seq(min(0, min(pp_samples, na.rm=T)), max(pp_samples, na.rm=T), bin_width)
  bins <- matrix(0, length(breaks), length(pp_x))
  rownames(bins) <- 1:nrow(bins)
  tb <- table(findInterval(pp_samples, breaks), rep(pp_x, each=nrow(pp_samples)))
  bins[rownames(bins) %in% rownames(tb),] <- tb
  bins <- t(t(bins) / apply(bins, 2, max, na.rm=T))
  bins[is.na(bins)] <- 0
  
  df <- data.frame(x.vals=rep(pp_x, each=nrow(bins)), y.vals=rep(breaks, ncol(bins)), 
                   value=c(bins))
  ggplot() +
    geom_raster(aes(x=df$x.vals, y=df$y.vals, fill = df$value), interpolate = interpolate) +
    geom_point(aes(y=y_data, x=x_data, color='black', group=exclusions, shape=exclusions), size=cex) +
    scale_fill_continuous(na.value = "transparent", high = "#8EB0D2", low = if(bin_width != 1){"white"}else{"white"},
                          name = 'Relative posterior density', guide = 'colorbar', breaks=c(0,1),
                          labels=c('0', '1')) +
    scale_colour_identity(name='', guide='legend', labels=c('Observed')) +
    scale_shape_manual(values=c(16, 1), guide='none') +
    labs(x = x_lab, y = y_lab) +
    theme_classic()
}


