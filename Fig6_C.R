library(tidyverse)

aadb <- data.frame(read.csv('Databases/Forams_AA.csv')) %>%
  select(!Non.polar & !Polar & !Neg.charge & !Pos.charge) %>%
  pivot_longer(!Taxon & !Seq & !TPM, values_to = 'Count', names_to = 'AA') %>%
  mutate(TPM = log10(TPM)) %>%
  filter(TPM != Inf & TPM != -Inf) %>%
  mutate(AA = ifelse(AA %in% c('F', 'Y', 'M', 'I', 'N', 'K'), 'FYMINK',
                      ifelse(AA %in% c('G', 'A', 'R', 'P'), 'GARP', 'Other'))) %>%
  group_by(Seq, TPM, Taxon, AA) %>%
  summarize(Count = sum(Count)) %>%
  na.omit()

sigs <- matrix(nrow = length(unique(aadb$Taxon))*20, ncol = 4)
for(i in 1:length(unique(aadb$Taxon))){
  for(j in 1:3){
    spec = unique(aadb$Taxon)[i]
    aa = unique(aadb$AA)[j]
    spec_aa_db <- filter(aadb, Taxon == spec & AA == aa)
    
    spec.fit <- lm(Count ~ TPM, data = spec_aa_db)
    spec.p <- summary(spec.fit)$coefficients[2,4]
    spec.slope <- summary(spec.fit)$coefficients[2,1]
    
    sigs[(i-1)*20+j, ] <- c(spec, aa, spec.p, spec.slope)
  }
}

sigs <- data.frame(sigs)
colnames(sigs) <- c('Taxon', 'AA', 'P', 'Slope')

sigs <- sigs %>%
  mutate(Sig = ifelse(P < 0.05, T, F))

gcdb <- data.frame(read.csv('Databases/Forams_GC3_ENc_TPM.csv')) %>%
  mutate(Taxon = paste(substr(Seq, 1, 4), substr(Seq,6,10), sep = '')) %>%
  group_by(Taxon) %>%
  summarize(Mean.GC4 = mean(GC4)) %>%
  na.omit()

sigs <- sigs %>%
  merge(gcdb, by = 'Taxon', all = T)

sigs <- na.omit(sigs)

ggplot(sigs, aes(Mean.GC4, as.numeric(Slope)*100, fill = AA)) +
  geom_point(size = 3, pch = 21, stroke = 1.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 50, linetype = 'dashed') +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_fill_manual(values = c('darkgreen', 'purple', 'gray')) +
  scale_color_manual(values = c('NA', 'blue')) +
  scale_x_continuous(limits = c(0, 80)) +
  #scale_shape_manual(values = c(19, 21)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_blank(),
    #axis.title = element_text(size = 20),
    legend.position = 'none'
  ) +
  labs(x = '', y = '')







