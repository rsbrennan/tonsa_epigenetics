# mbias plotting

library(scmeth)

mbiasList <- mbiasplot(dir= "~/tonsa_epigenetics/analysis/methylation_extract/bias_check/")

mbiasDf <- do.call(rbind.data.frame, mbiasList)
meanTable <- stats::aggregate(methylation ~ position + read, data=mbiasDf, FUN=mean)
sdTable <- stats::aggregate(methylation ~ position + read, data=mbiasDf, FUN=sd)
seTable <- stats::aggregate(methylation ~ position + read, data=mbiasDf, FUN=function(x){sd(x)/sqrt(length(x))})
sum_mt<-data.frame('position'=meanTable$position,'read'=meanTable$read,
                       'meth'=meanTable$methylation, 'sdMeth'=sdTable$methylation,
                       'seMeth'=seTable$methylation)

sum_mt$upperCI <- sum_mt$meth + (1.96*sum_mt$seMeth)
sum_mt$lowerCI <- sum_mt$meth - (1.96*sum_mt$seMeth)
sum_mt$read_rep <- paste(sum_mt$read, sum_mt$position, sep="_")

g <- ggplot2::ggplot(sum_mt)
g <- g + ggplot2::geom_line(ggplot2::aes_string(x='position', y='meth',
                                                colour='read'))
g <- g + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = 'lowerCI',
                        ymax = 'upperCI', x='position', fill = 'read'),
                        alpha=0.4)
g <- g  + ggplot2::ggtitle('Mbias Plot')
g <- g + ggplot2::ylab('methylation') + theme_bw()

ggsave("~/tonsa_epigenetics/figures/mbias_all.png", g, h=3, w=5, units="in")


