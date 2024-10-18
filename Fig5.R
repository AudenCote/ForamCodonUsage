require(tidyverse)
require(grid)
require(gridExtra)

f_gc4 <- data.frame(read.csv('Databases/Forams_GC3_ENc_TPM.csv'))

f_gc12 <- data.frame(read.csv('Databases/Forams_GC12.csv')) %>%
  mutate(Taxon = paste0(substr(Seq, 1, 4), substr(Seq, 6, 10)))

f_compdb <- f_gc4 %>%
  merge(f_gc12, by = 'Seq')

neutdb <- f_compdb %>%
  mutate(GC12 = GC12 * 100)

tpm_thresholds <- neutdb %>%
  group_by(Taxon) %>%
  summarize(TPM.Lower = quantile(TPM, 0.1), TPM.Higher = quantile(TPM, 0.9))

tpmdb <- neutdb %>%
  group_by(Taxon) %>%
  summarize(Mean.GC4 = mean(GC4), Mean.GC12 = mean(GC12), NSeqs = n())

low_tpm_neutdb <- neutdb %>%
  merge(tpm_thresholds, by = 'Taxon') %>%
  filter(TPM < TPM.Lower)

hi_tpm_neutdb <- neutdb %>%
  merge(tpm_thresholds, by = 'Taxon') %>%
  filter(TPM > TPM.Higher)

slopes <- matrix(nrow = length(unique(neutdb$Taxon)), ncol = 7)
for(i in 1:length(unique(neutdb$Taxon))){
  spec = unique(neutdb$Taxon)[i]
  
  lo.specdb <- filter(low_tpm_neutdb, Taxon == spec)
  hi.specdb <- filter(hi_tpm_neutdb, Taxon == spec)
  
  low.spec.fit <- lm(GC12 ~ GC4, data = lo.specdb)
  low.spec.slope <- summary(low.spec.fit)$coefficients[2,1]
  low.spec.p.pearson <- summary(low.spec.fit)$coefficients[2,4]
  low.spec.p.spearman <- cor.test(lo.specdb$GC4, lo.specdb$GC12, method = 'spearman')$p.value
  
  hi.spec.fit <- lm(GC12 ~ GC4, data = hi.specdb)
  hi.spec.slope <- summary(hi.spec.fit)$coefficients[2,1]
  hi.spec.p.pearson <- summary(hi.spec.fit)$coefficients[2,4]
  hi.spec.p.spearman <- cor.test(hi.specdb$GC4, hi.specdb$GC12, method = 'spearman')$p.value
  
  slopes[i, ] <- c(spec, low.spec.slope, low.spec.p.pearson, low.spec.p.spearman, hi.spec.slope, hi.spec.p.pearson, hi.spec.p.spearman)
}

slopes <- data.frame(slopes)
colnames(slopes) <- c('Taxon', 'Slope.Low', 'P.Pearson.Low', 'P.Spearman.Low', 'Slope.Hi', 'P.Pearson.Hi', 'P.Spearman.Hi')

tpmdb <- tpmdb %>%
  merge(slopes, by = 'Taxon') %>%
  mutate(Low.Sig = ifelse(as.numeric(P.Spearman.Low) < 0.05, T, F), Hi.Sig = ifelse(as.numeric(P.Spearman.Hi) < 0.05, T, F))


add_style <- function(x){
  x +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_y_continuous(limits = c(-0.85, 0.85), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 85), expand = c(0, 0)) +
    scale_shape_manual(values = c(21, 19)) +
    scale_color_manual(values = c('black', 'blue')) +
    theme_classic() +
    theme(
      legend.position = 'none',
      axis.text = element_text(size = 12, color = 'black'),
      axis.title = element_blank()
    )
}

low_plot <- ggplot(tpmdb, aes(Mean.GC4, as.numeric(Slope.Low), shape = Low.Sig)) %>%
  add_style()

hi_plot <- ggplot(tpmdb, aes(Mean.GC4, as.numeric(Slope.Hi), shape = Hi.Sig)) %>%
  add_style()

grid.arrange(low_plot, hi_plot, hi_plot, nrow = 1)















