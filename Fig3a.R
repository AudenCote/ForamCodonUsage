library(tidyverse)

codondb <- data.frame(read.csv('Databases/Forams_Codons.csv'))

gcdb <- data.frame(read.csv('Databases/Forams_GC3_ENc_TPM.csv')) %>%
  select(!TPM)

taxdb <- data.frame(read_tsv('Databases/TaxDB.tsv'))

aa_sums <- codondb %>%
  group_by(Tax, AA) %>%
  summarize(TotalAA = sum(Count))

aa_cats <- codondb %>%
  group_by(AA) %>%
  summarize(Family = n()/length(unique(codondb$Tax)))

cons_aa_order <- codondb %>%
  group_by(AA) %>%
  summarize(TotalAA = sum(Count))

cons_aa_order <- aa_sums[order(-cons_aa_order$TotalAA),]

gc4_order <- gcdb %>%
  mutate(Tax = paste(substr(Seq, 1, 4), substr(Seq, 6, 10), sep = '')) %>%
  group_by(Tax) %>%
  summarize(Mean.GC4 = mean(GC4))

gc4_order <- gc4_order[order(gc4_order$Mean.GC4),]

codondb <- codondb %>%
  merge(aa_sums, by = c('Tax', 'AA')) %>%
  merge(gc4_order, by = 'Tax') %>%
  merge(taxdb, by = 'Tax') 

aa_order = c('K', 'I', 'T', 'R', 'N', 'S', 'M', 'L', 'Y', 'F', 'C', 'W', 'Q', 'P', 'H', 'E', 'V', 'A', 'G', 'D')
codon_order = c('AAA', 'ATA', 'ACA', 'AGA', 'AAT', 'ATT', 'ACT', 'AGT', 'AAC', 'ATC', 'ACC', 'AGC', 'AAG', 'ATG', 'ACG', 'AGG', 'TTA', 'TCA', 'TAT', 'TTT', 'TCT', 'TGT', 'TAC', 'TTC', 'TCC', 'TGC', 'TTG', 'TCG', 'TGG', 'CAA', 'CTA', 'CCA', 'CGA', 'CAT', 'CTT', 'CCT', 'CGT', 'CAC', 'CTC', 'CCC', 'CGC', 'CAG', 'CTG', 'CCG', 'CGG', 'GAA', 'GTA', 'GCA', 'GGA', 'GAT', 'GTT', 'GCT', 'GGT', 'GAC', 'GTC', 'GCC', 'GGC', 'GAG', 'GTG', 'GCG', 'GGG')

codondb$Codon <- factor(codondb$Codon, levels = codon_order)
codondb$AA <- factor(codondb$AA, levels = aa_order)

aadb <- codondb %>%
  group_by(AA) %>%
  summarize(Codons = n_distinct(Codon)) %>%
  filter(Codons > 1) %>%
  mutate(bin1 = 1/(Codons*5), bin2 = 2/(Codons*5), bin3 = 3/(Codons*5), bin4 = 4/(Codons*5), bin5 = 5/(Codons*5), bin6 = 1/Codons + (1-(1/Codons))/5, bin7 = 1/Codons + 2*(1-(1/Codons))/5, bin8 = 1/Codons + 3*(1-(1/Codons))/5, bin9 = 1/Codons +4*(1-(1/Codons))/5) %>%
  select(!Codons)

codondb <- codondb %>%
  merge(aadb, by = 'AA') %>%
  mutate(CodonProp = (Count/TotalAA)) %>%
  mutate(Bin = ifelse(CodonProp < bin1, 1,
                      ifelse(CodonProp < bin2, 2,
                             ifelse(CodonProp < bin3, 3,
                                    ifelse(CodonProp < bin4, 4,
                                           ifelse(CodonProp < bin5, 5,
                                                  ifelse(CodonProp < bin6, 6,
                                                         ifelse(CodonProp < bin7, 7,
                                                                ifelse(CodonProp < bin8, 8,
                                                                       ifelse(CodonProp < bin9, 9, 10)))))))))) %>%
  mutate(Whole_Name = paste0(Whole_Name, ' (', str_sub(Tax, 6, 10), ')')) %>%
  merge(aa_cats, by = 'AA')

ggplot(codondb, aes(Codon, reorder(Whole_Name, Mean.GC4), fill = as.factor(Bin))) +
  geom_tile(color = 'white') +
  theme_classic() +
  scale_fill_manual(values = c("#0b389d", "#6581bf", "#a1b3d7", '#d7dcea', "white", "white", '#F6CDCD', '#F29F9F', '#EE7272', '#E61717')) +
  theme(
    axis.text.x = element_text(color = 'black', size = 15, angle = 90, vjust = 0.5, family = 'Monaco'),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 15),
    legend.position = 'none',
    axis.title = element_blank(),
    strip.background = element_blank()#,
  ) +
  facet_grid(~Family + AA, scales = "free", space = "free")

gc4_sidebar <- ggplot(codondb, aes('1', reorder(Whole_Name, as.numeric(Mean.GC4)), fill = as.numeric(Mean.GC4))) +
  geom_tile(color = 'white') +
  theme_classic() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient(low = 'white', high = 'black') +
  theme(legend.position = 'none')










