require(tidyverse)

aadb <- data.frame(read.csv('Databases/Forams_AA.csv')) %>%
  select(!Neg.charge & !Non.polar & !Polar & !Pos.charge)

tpmdb <- aadb %>%
  group_by(Taxon) %>%
  summarize(Low.TPM = quantile(TPM, 0.2), High.TPM = quantile(TPM, 0.8))

gcdb <- data.frame(read.csv('Databases/Forams_GC3_ENc_TPM.csv')) %>%
  select(Seq | GC4 | TPM) %>%
  mutate(Taxon = paste0(substr(Seq, 1, 4), substr(Seq, 6, 10))) %>%
  merge(tpmdb, by = 'Taxon') %>%
  mutate(Expression = ifelse(TPM < Low.TPM, 'Low',
                                    ifelse(TPM > High.TPM, 'High', 'Middle')))

gcdb_by_exp <- gcdb %>%
  group_by(Taxon, Expression) %>%
  summarize(Mean.GC4 = mean(GC4))

gcdb <- gcdb %>%
  group_by(Taxon) %>%
  summarize(Mean.GC4 = mean(GC4))

utrdb <- data.frame(read.csv('Databases/Forams_UTRs.csv')) %>%
  mutate(Taxon = paste0(substr(Seq, 1, 4), substr(Seq, 6, 10))) %>%
  filter(Len > 50 & Len < 500) %>%
  mutate(UTR.GC = (G + C)/(A + T + G + C)) %>%
  group_by(Taxon) %>%
  summarize(UTR.GC = mean(UTR.GC) * 100)

aa_names = data.frame(AA = c('F', 'Y', 'M', 'I', 'N', 'K', 'G', 'A', 'R', 'P', 'C', 'D', 'E', 'H', 'L', 'Q', 'S', 'T', 'V', 'W'),
                      AA.name = c('Phenylalanine', 'Tyrosine', 'Methionine', 'Isoleucine', 'Asparagine', 'Lysine', 'Glycine', 'Alanine', 'Arginine', 'Proline', 'Cysteine', 'Aspartic acid', 'Glutamic acid', 'Histidine', 'Leucine', 'Glutamine', 'Serine', 'Threonine', 'Valine', 'Tryptophan'))

aadb <- aadb %>%
  pivot_longer(!Seq & !Taxon & !TPM, values_to = 'Freq', names_to = 'AA') %>%
  merge(tpmdb, by = 'Taxon') %>%
  mutate(Expression = ifelse(TPM < Low.TPM, 'Low', 
                                    ifelse(TPM > High.TPM, 'High', 'Middle'))) 

aadb_by_exp <- aadb %>%
  group_by(Taxon, AA, Expression) %>%
  summarize(Mean.Freq = mean(Freq) * 100) %>%
  merge(gcdb_by_exp, by = c('Taxon', 'Expression')) %>%
  merge(aa_names, by = 'AA') %>%
  mutate(AA.Group = ifelse(AA %in% c('G', 'A', 'R', 'P'), 'GARP',
                           ifelse(AA %in% c('F', 'Y', 'M', 'I', 'N', 'K'), 'FYMINK', 'Other')))
  
aadb <- aadb %>%
  group_by(Taxon, AA) %>%
  summarize(Mean.Freq = mean(Freq) * 100) %>%
  merge(gcdb, by = 'Taxon') %>%
  merge(utrdb, by = 'Taxon') %>%
  merge(aa_names, by = 'AA') %>%
  mutate(AA.Group = ifelse(AA %in% c('G', 'A', 'R', 'P'), 'GARP',
                           ifelse(AA %in% c('F', 'Y', 'M', 'I', 'N', 'K'), 'FYMINK', 'Other')))

aadb$AA.name = factor(aadb$AA.name, levels = aa_names$AA.name)

panel_plot <- ggplot(aadb) +
  geom_point(aes(Mean.GC4, Mean.Freq, color = AA.Group)) +
  geom_smooth(mapping = aes(Mean.GC4, Mean.Freq), method = 'lm', color = 'black', se = F) +
  geom_smooth(mapping = aes(UTR.GC, Mean.Freq), method = 'lm', color = 'darkred', se = F, linetype = 'dashed') +
  scale_color_manual(values = c('purple', 'darkgreen', 'gray')) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = 'bold'),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.position = 'none'
  ) +
  facet_wrap(~AA.name, scales = 'free_y') +
  labs(x = 'Average silent-site GC content (%)', y = 'Average amino acid frequency (%)')


sigs <- matrix(nrow = 60, ncol = 4)
for(i in 1:3){
  for(j in 1:20){
    exp = unique(aadb_by_exp$Expression)[i]
    aa = unique(aadb_by_exp$AA)[j]
    exp_aa_db <- filter(aadb_by_exp, Expression == exp & AA == aa)
    
    spec.fit <- lm(Mean.Freq ~ Mean.GC4, data = exp_aa_db)
    spec.p <- summary(spec.fit)$coefficients[2,4]
    spec.slope <- summary(spec.fit)$coefficients[2,1]
      
    sigs[(i-1)*20+j, ] <- c(exp, aa, spec.p, spec.slope)
  }
}
  
sigs <- data.frame(sigs)
colnames(sigs) <- c('Expression', 'AA', 'P', 'Slope')

sigs_wide <- sigs %>%
  select(!P) %>%
  pivot_wider(names_from = 'Expression', values_from = 'Slope') %>%
  mutate(AA.Group = ifelse(AA %in% c('G', 'A', 'R', 'P'), 'GARP',
                           ifelse(AA %in% c('F', 'Y', 'M', 'I', 'N', 'K'), 'FYMINK', 'Other')))

xy <- data.frame(x = c(-0.05, 0.05),
                 y = c(-0.05, 0.05))

slopes_plot <- ggplot() +
  geom_point(data = sigs_wide, mapping = aes(as.numeric(Low), as.numeric(High), fill = AA.Group), size = 3.5, pch = 21, color = 'NA') +
  geom_line(data = xy, mapping = aes(as.numeric(x), as.numeric(y)), size = 1, linetype = 'dashed') +
  scale_fill_manual(values = c('purple', 'darkgreen', 'gray')) +
  scale_x_continuous(limits = c(-0.055, 0.055), expand = c(0, 0), breaks = c(-0.05, 0, 0.05)) +
  scale_y_continuous(limits = c(-0.055, 0.055), expand = c(0, 0), breaks = c(-0.05, 0, 0.05)) +
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_blank()
  ) +
  labs(x = '', y = '', color = '')

slopes_plot


