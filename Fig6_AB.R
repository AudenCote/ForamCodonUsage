library(tidyverse)

tpmdb <- data.frame(read.csv('Databases/Forams_GC3_ENc_TPM.csv')) %>%
  mutate(Taxon = paste(substr(Seq, 1, 4), substr(Seq,6,10), sep = '')) %>%
  merge(data.frame(read.csv('Databases/Forams_GC12.csv')), by = 'Seq') %>%
  mutate(GC12 = GC12*100) %>%
  mutate(LogTPM = log10(TPM)) %>%
  na.omit() %>%
  filter(LogTPM != -Inf & LogTPM != Inf) %>%
  select(Taxon | Seq | GC12 | GC4 | LogTPM | ObsWrightENc_6Fold)

sigs <- matrix(nrow = length(unique(tpmdb$Taxon)), ncol = 7)
for(i in 1:length(unique(tpmdb$Taxon))){
  spec = unique(tpmdb$Taxon)[i]
  specdb <- filter(tpmdb, Taxon == spec)
  
  gc4.spec.fit <- lm(GC4 ~ LogTPM, data = specdb)
  gc4.spec.p <- cor.test(specdb$GC4, specdb$LogTPM, method = 'spearman')$p.value
  gc4.spec.slope <- summary(gc4.spec.fit)$coefficients[2,1]
  
  gc12.spec.fit <- lm(GC12 ~ LogTPM, data = specdb)
  gc12.spec.p <- cor.test(specdb$GC12, specdb$LogTPM, method = 'spearman')$p.value
  gc12.spec.slope <- summary(gc12.spec.fit)$coefficients[2,1]
  
  enc.spec.fit <- lm(ObsWrightENc_6Fold ~ LogTPM, data = specdb)
  enc.spec.p <- cor.test(specdb$ObsWrightENc_6Fold, specdb$LogTPM, method = 'spearman')$p.value
  enc.spec.slope <- summary(enc.spec.fit)$coefficients[2,1]
  
  sigs[i, ] <- c(spec, gc4.spec.p, gc4.spec.slope, gc12.spec.p, gc12.spec.slope, enc.spec.p, enc.spec.slope)
}

sigs <- data.frame(sigs)
colnames(sigs) <- c('Taxon', 'GC4.P', 'GC4.Slope', 'GC12.P', 'GC12.Slope', 'ENc.P', 'ENc.Slope')

tpmdb <- tpmdb %>%
  group_by(Taxon) %>%
  summarize(Mean.GC4 = mean(GC4), Mean.GC12 = mean(GC12)) %>%
  merge(sigs, by = 'Taxon') %>%
  mutate(GC4.Sig = ifelse(as.numeric(GC4.P) < 0.05, 'P < 0.05', 'P > 0.05')) %>%
  mutate(GC12.Sig = ifelse(as.numeric(GC12.P) < 0.05, 'P < 0.05', 'P > 0.05')) %>%
  mutate(ENc.Sig = ifelse(as.numeric(ENc.P) < 0.05, 'P < 0.05', 'P > 0.05')) %>%
  mutate(NonForam = ifelse(Taxon %in% nonforams$Taxon, T, F)) %>%
  filter(NonForam == F)

gc12_min_expected = 43.05557
gc12_max_expected = 47.16760

ggplot(tpmdb, aes(Mean.GC12, as.numeric(GC12.Slope), shape = GC12.Sig, color = NonForam)) +
  geom_point(size = 3, stroke = 0.75) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = gc12_min_expected, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = gc12_max_expected, linetype = 'dashed', color = 'red') +
  scale_shape_manual(values = c(19,21)) +
  scale_color_manual(values = c('black', 'blue')) +
  scale_y_continuous(lim = c(-10.5, 10.5)) +
  scale_x_continuous(lim = c(35, 55)) +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_blank()
  ) +
  labs(x = '', y = '')

