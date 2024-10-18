require(tidyverse)
require(grid)
require(gridExtra)

f_gc4 <- data.frame(read.csv('Databases/Forams_GC3_ENc_TPM.csv')) %>%
  mutate(Taxon = paste0(substr(Seq, 1, 4), substr(Seq, 6, 10))) %>%
  group_by(Taxon) %>%
  summarize(Mean.GC4 = mean(GC4))

f_gc12 <- data.frame(read.csv('Databases/Forams_GC12.csv')) %>%
  mutate(Taxon = paste0(substr(Seq, 1, 4), substr(Seq, 6, 10))) %>%
  group_by(Taxon) %>%
  summarize(Mean.GC12 = mean(GC12)*100)

f_utrs <- data.frame(read.csv('Databases/Forams_UTRs.csv')) %>%
  mutate(Taxon = paste0(substr(Seq, 1, 4), substr(Seq, 6, 10))) %>%
  group_by(Taxon) %>%
  summarize(A = sum(A), T = sum(T), G = sum(G), C = sum(C)) %>%
  mutate(Total = A + T + G + C) %>%
  mutate(Mean.UTR.GC = 100*(G + C)/(A + T + G + C)) %>%
  select(Taxon | Mean.UTR.GC)

f_compdb <- f_gc4 %>%
  merge(f_gc12, by = 'Taxon') %>%
  merge(f_utrs, by = 'Taxon')

compdb <- f_combdb

add_style <- function(x){
  x +
    geom_point(size = 3, pch = 21, stroke = 0.75) +
    geom_smooth(method = 'lm', color = 'black', se = F, size = 0.5, fullrange = F) +
    geom_abline(slope = 1, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 80), expand = c(0, 0)) +
    scale_color_manual(values = c('black', 'blue')) +
    theme_classic() +
    theme(
      legend.position = 'none',
      axis.text = element_text(size = 12, color = 'black'),
      axis.title = element_blank()
    )
}

fig_4a <- ggplot(compdb, aes(Mean.UTR.GC, Mean.GC4)) %>%
  add_style

fig_4b <- ggplot(compdb, aes(Mean.UTR.GC, Mean.GC12)) %>%
  add_style()

fig_4c <- ggplot(compdb, aes(Mean.GC4, Mean.GC12)) %>%
  add_style()

grid.arrange(fig_4a, fig_4b, fig_4c, nrow = 1)


















