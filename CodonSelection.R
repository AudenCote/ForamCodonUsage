library(tidyverse)
library(FactoMineR)
library(factoextra)
library(data.table)
library(ggrepel)
library(grid)
library(gridExtra)
library(gt)
library(dplyr)

exp <- data.frame(read.csv('Databases/Forams_GC3_ENc_TPM.csv')) %>%
  select(Seq | TPM) %>%
  mutate(LogTPM = log10(TPM))

codon_by_aa <- data.frame(read_tsv('Databases/CodonDB.tsv'))

taxdb <- data.frame(read_tsv('Databases/TaxDB.tsv')) %>%
  rename(Taxon = Tax)

freqs <- data.frame(read.csv('Databases/Forams_WithinFamilyCodonFreqs_PerSeq.csv')) %>%
  mutate(Taxon = paste(substr(Seq, 1, 4), substr(Seq,6,10), sep = '')) %>%
  merge(exp, by = 'Seq')

freqs <- freqs %>% filter(TPM > 0)

sigs <- matrix(nrow = length(unique(freqs$Taxon)) * length(unique(freqs$Codon)), ncol = 4)
for(i in 1:length(unique(freqs$Taxon))){
  spec = unique(freqs$Taxon)[i]
  specdb <- filter(freqs, Taxon == spec)
  
  for(j in 1:length(unique(specdb$Codon))){
    codon = unique(specdb$Codon)[j]
    codondb = filter(specdb, Codon == codon) %>% na.omit()
    
    cod.fit <- lm(Freq ~ LogTPM, data = codondb)
    cod.p <- cor.test(codondb$Freq, codondb$LogTPM, method = 'spearman', exact = F)$p.value
    cod.slope <- summary(cod.fit)$coefficients[2,1]
    
    sigs[(i-1)*length(unique(freqs$Codon)) + j, ] <- c(spec, codon, cod.slope, cod.p)
  }
  
}

sigs <- data.frame(sigs)
colnames(sigs) <- c('Taxon', 'Codon', 'Slope', 'P')

corrdb <- sigs %>%
  na.omit() %>%
  merge(codon_by_aa, by = 'Codon')

aa_order = c('K', 'I', 'T', 'R', 'N', 'S', 'M', 'L', 'Y', 'F', 'C', 'W', 'Q', 'P', 'H', 'E', 'V', 'A', 'G', 'D')
codon_order = c('AAA', 'ATA', 'ACA', 'AGA', 'AAT', 'ATT', 'ACT', 'AGT', 'AAC', 'ATC', 'ACC', 'AGC', 'AAG', 'ATG', 'ACG', 'AGG', 'TTA', 'TCA', 'TAT', 'TTT', 'TCT', 'TGT', 'TAC', 'TTC', 'TCC', 'TGC', 'TTG', 'TCG', 'TGG', 'CAA', 'CTA', 'CCA', 'CGA', 'CAT', 'CTT', 'CCT', 'CGT', 'CAC', 'CTC', 'CCC', 'CGC', 'CAG', 'CTG', 'CCG', 'CGG', 'GAA', 'GTA', 'GCA', 'GGA', 'GAT', 'GTT', 'GCT', 'GGT', 'GAC', 'GTC', 'GCC', 'GGC', 'GAG', 'GTG', 'GCG', 'GGG')

corrdb$Codon <- factor(corrdb$Codon, levels = codon_order)
corrdb$AA <- factor(corrdb$AA, levels = aa_order)

gcdb <- data.frame(read.csv('../Databases/Forams_GC3_ENc_TPM.csv')) %>%
  select(!TPM)

gc4_order <- gcdb %>%
  mutate(Taxon = paste(substr(Seq, 1, 4), substr(Seq, 6, 10), sep = '')) %>%
  group_by(Taxon) %>%
  summarize(Mean.GC4 = mean(GC4))

corrdb$P.Adj <-p.adjust(corrdb$P, method = 'BH')

corrdb <- corrdb %>%
  merge(gc4_order, by = 'Taxon') %>%
  mutate(Sig = ifelse(P.Adj > 0.05, F, T)) %>%
  merge(taxdb, by = 'Taxon') 

aa_cats <- corrdb %>%
  group_by(AA) %>%
  summarize(Family = n()/length(unique(codondb$Tax)))

corrdb <- corrdb %>%
  merge(aa_cats, by = 'AA')

corrdb$Slope = as.numeric(corrdb$Slope)

corrdb <- corrdb %>%
  mutate(bin = ifelse(Slope < -0.15, 1,
                      ifelse(Slope < -0.1, 2,
                             ifelse(Slope < -0.05, 3,
                                    ifelse(Slope < 0, 4,
                                           ifelse(Slope == 0, 5,
                                                  ifelse(Slope < 0.05, 6,
                                                         ifelse(Slope < 0.1, 7,
                                                                ifelse(Slope < 0.15, 8, 9)))))))))

corrdb <- corrdb %>%
  filter(AA != 'M' & AA != 'W')

fig3b <- ggplot(corrdb, aes(Codon, reorder(Taxon, Mean.GC4), fill = as.factor(bin))) +
  geom_tile(color = 'white') +
  theme_classic() +
  scale_fill_manual(values = c("#0b389d", "#6581bf", "#a1b3d7", '#d7dcea', "white", '#F6CDCD', '#F29F9F', '#EE7272', '#E61717')) +
  theme(
    axis.text.x = element_text(color = 'black', size = 15, angle = 90, vjust = 0.5, family = 'Monaco'),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 15),
    legend.position = 'none',
    axis.title = element_blank(),
    strip.background = element_blank()#,
  ) +
  facet_grid(~ Family + AA, scales = "free", space = "free")


codonfreqs <- data.frame(read.csv('../Databases/Forams_Codons.csv')) %>%
  mutate(Freq = Count/TotalCodons) %>%
  select(Tax | Codon | Freq) %>%
  rename(Taxon = Tax)

cadb <- corrdb %>%
  merge(codonfreqs, by = c('Taxon', 'Codon')) %>%
  mutate(Slope = Slope + 1) %>%
  select(Taxon | Codon | Slope) %>%
  pivot_wider(names_from = 'Codon', values_from = 'Slope') %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Taxon')

res.ca <- CA(cadb, ncp = 5, graph = TRUE)

ca_labels.x = fviz_ca_row(res.ca, geom = "point")$labels$x
ca_labels.y = fviz_ca_row(res.ca, geom = "point")$labels$y

ca_plot_data <- data.table(fviz_ca_row(res.ca, geom = "point")$data) %>%
  mutate(Taxon = name) %>%
  select(Taxon | x | y) %>%
  merge(select(corrdb, Taxon | Mean.GC4 | Clade), by = 'Taxon') %>%
  mutate(Clade = ifelse(Clade == 'Xenophyophorea' | Clade == 'Monothalamid', 'Monothalamids',
                        ifelse(Clade == 'Tubothalamid', 'Tubothalamea', Clade)))

figs14 <- rscu_plot <- ggplot() +
  geom_point(data = ca_plot_data, aes(x, y, color = Clade), size = 3) +
  scale_color_manual(values = c('red', 'darkblue', 'orange', 'lightblue')) +
  theme_classic() +
  labs(x = ca_labels.x, y = ca_labels.y) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  labs(color = 'Clade')

rscu_plot









