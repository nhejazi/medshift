library(tidyverse)
library(data.table)
library(ggplot2)
library(ggthemes)
library(grid)

res <- list()
for(i in 1:1200){
    ss <- try(load(paste0('output/out', i, '.rda')))
    if(inherits(ss, 'try-error')) {
        ## cat('error ', i)
    } else {
        rtemp <- get(paste0('r', i))
        if(any(sapply(rtemp, length) == 1)) print(i)
        res <- c(res, rtemp)
    }
}

out <- rbindlist(res[sapply(res, length) > 1])

## out %>% filter(estimator == 'TMLE', n == 16200, parameter == 'DE', type == 1) %>%
##     pull(estimate) %>% hist()

## out %>% filter(estimator == 'TMLE', n == 16200, parameter == 'DE', type == 4) %>%
##     pull(estimate) %>% mean()

load('true.rda')
## source('utils.r')
## true <- truth(1e7)
## save(true, file = 'true.rda')
alpha <- 0.05
dplot <- out %>% left_join(true, by = 'parameter') %>%
    group_by(type, parameter, estimator, n) %>%
    summarise(## m = n(),
        rbias = abs(mean(estimate - truth)),
        rnbias = abs(mean(sqrt(n) * (estimate - truth))),
        relse = mean(ses) / sd(estimate),
        relsd = sd(sqrt(n) * estimate / sqrt(eff_bound)),
        relrmse = sqrt(mean(n * (estimate - truth)^2  / eff_bound)),
        coverage = mean(qnorm(alpha / 2) < (estimate - truth) / ses &
                        (estimate - truth) / ses < qnorm(1 - alpha / 2))) %>%
    ## filter(type ==  0, estimator == 'os') %>%
    gather(statistic, value, rbias, rnbias, relse, relsd, relrmse, coverage) %>%
    ungroup() %>%
    mutate(parameter = factor(parameter),
           type = factor(type),
           estimator = factor(estimator),
           statistic = factor(statistic,
                              levels = c('rbias', 'rnbias',
                                         'relsd', 'relrmse',
                                         'coverage', 'relse'),
                              labels = c("group('|',Bias,'|')",
                                         "n^{1/2}~group('|',Bias,'|')",
                                         'rel~sd(hat(theta))',
                                         '(n%*%rel~MSE)^{1/2}',
                                         'Cov(0.95)',
                                         'hat(sd)(hat(theta))/sd(hat(theta))')),
           type = factor(type,
                         levels = 0:5,
                         labels = c('none', 'g', 'e', 'm', 'b', 'd')))

dummy <- data.frame(statistic = c("n^{1/2}~group('|',Bias,'|')",
                                  "group('|',Bias,'|')",
                                  'hat(sd)(hat(theta))/sd(hat(theta))',
                                  'rel~sd(hat(theta))',
                                  '(n%*%rel~MSE)^{1/2}',
                                  'Cov(0.95)'),
                    lim = c(0, 0, 1, 1, 1, 1 - alpha))

ns <- cumsum(rep(sqrt(200), 9))^2
split1 <- c('none')
split2 <- c('g','b')
split3 <- c('m','e','d')

plot1 <- dplot %>% rename(Parameter = parameter, Estimator = estimator) %>%
    filter(type %in% split1) %>%
    ggplot(aes(n, value, colour = Estimator, shape = Parameter)) +
    ylab('') +
    geom_hline(data = dummy, aes(yintercept = lim), colour = 'darkgray', size = 1.5) +
    geom_point(size = 2) +
    geom_line(linetype = "dotted", size = 0.5) +
    theme_igray() +
    scale_x_continuous(breaks = ns, trans = 'sqrt') +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 50, hjust = 1),
          text = element_text(size = 18), panel.spacing = unit(1, "lines"),
          strip.text.y = element_blank()) +
    ## scale_shape_manual(name = 'Estimator:',
    ##                    values = c(3, 4, 1, 2)) +
    facet_grid(statistic ~ type, scales = 'free_y', labeller = label_parsed)
plot2 <- dplot %>% rename(Parameter = parameter, Estimator = estimator) %>%
    filter(type %in% split3) %>%
    ggplot(aes(n, value, colour = Estimator, shape = Parameter)) +
    ylab('') +
    geom_hline(data = dummy, aes(yintercept = lim), colour = 'darkgray', size = 1.5) +
    geom_point(size = 2) +
    geom_line(linetype = "dotted", size = 0.5) +
    theme_igray() +
    scale_x_continuous(breaks = ns, trans = 'sqrt') +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 50, hjust = 1),
          text = element_text(size = 18), panel.spacing = unit(1, "lines"),
          strip.text.y = element_blank()) +
    ## scale_shape_manual(name = 'Estimator:',
    ##                    values = c(3, 4, 1, 2)) +
    facet_grid(statistic ~ type, scales = 'free_y', labeller = label_parsed)
plot3 <- dplot %>% rename(Parameter = parameter, Estimator = estimator) %>%
    filter(type %in% split2) %>%
    ggplot(aes(n, value, colour = Estimator, shape = Parameter)) +
    ylab('') +
    geom_hline(data = dummy, aes(yintercept = lim), colour = 'darkgray', size = 1.5) +
    geom_point(size = 2) +
    geom_line(linetype = "dotted", size = 0.5) +
    theme_igray() +
    scale_x_continuous(breaks = ns, trans = 'sqrt') +
    theme(legend.position = 'right', axis.text.x = element_text(angle = 50, hjust = 1),
          text = element_text(size = 18), panel.spacing = unit(1, "lines")) +
    ## scale_shape_manual(name = 'Estimator:',
    ##                    values = c(3, 4, 1, 2)) +
    facet_grid(statistic ~ type, scales = 'free_y', labeller = label_parsed)

pdf('plot.pdf', width = 13, height = 11)
grid.newpage()
grid.draw(cbind(ggplotGrob(plot1), ggplotGrob(plot2), ggplotGrob(plot3), size = 'last'))
dev.off()

dplot %>%
    filter(type == 'All-consistent', parameter == 'indirect', n == 9800,
           statistic == 'hat(sd)(hat(theta))/sd(hat(theta))')
