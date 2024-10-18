library(tidyverse)

recids <- data.frame(read.csv('Databases/UTR_RecID_Corr.csv'))
tpmdb <- data.frame(read.csv('Databases/Forams_GC3_ENc_TPM.csv')) %>%
  select(Seq | TPM | GC4)

meangc4s <- tpmdb %>%
  mutate(Taxon = paste(substr(Seq, 1, 4), substr(Seq, 6, 10), sep = '')) %>%
  group_by(Taxon) %>%
  summarize(Mean.GC4 = mean(GC4))

convar <- data.frame(read.csv('Databases/Forams_Akashi.csv')) %>%
  select(Seq | GC4 | Type) %>%
  pivot_wider(values_from = GC4, names_from = Type) %>%
  mutate(Taxon = paste(substr(Seq, 1, 4), substr(Seq, 6, 10), sep = '')) %>%
  merge(recids, by.x = 'Seq', by.y = 'UTR') %>%
  merge(tpmdb, by.x = 'ORF', by.y = 'Seq')

tpm_perc <- convar %>%
  group_by(Taxon) %>%
  summarize(TPM_P10 = quantile(TPM, 0.1), TPM_P90 = quantile(TPM, 0.90), TPM_P50 = quantile(TPM, 0.5), N = n()) %>%
  filter(N > 300)

convar <- convar %>%
  merge(tpm_perc, by = 'Taxon') %>%
  mutate(Perc = ifelse(TPM < TPM_P10, 'TPM < 10%', 
                       ifelse(TPM > TPM_P90, 'TPM > 90%', '10% < TPM > 90%'))) %>%
  na.omit() %>%
  filter(TPM < TPM_P10) %>%
  group_by(Taxon) %>%
  summarize(Variable = mean(Variable), Conserved = mean(Conserved)) %>%
  mutate(diff = log(Variable/Conserved)) %>%
  merge(meangc4s, by = 'Taxon')

ggplot() +
  geom_point(data = convar, mapping = aes(Mean.GC4, diff), size = 4, pch = 21, stroke = 1) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 50, linetype = 'dashed') +
  scale_y_continuous(limits = c(-0.1, 0.6)) +
  scale_x_continuous(limits = c(0, 75), breaks = c(0, 25, 50, 75)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15, color = 'black'),
    axis.title = element_text(size = 15)
  ) +
  labs(x = '', y = '')


