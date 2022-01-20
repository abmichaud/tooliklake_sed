library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)

# Significant taxa
metadatat<- read_excel(path="tl_metadata.xlsx")
tooll.tax <- read_tsv("tooll.tax") %>% 
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string = taxonomy, pattern = ";$", replacement = "")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")

tooll.otu.single <- read_tsv("tooll.shared", col_types = cols(Group=col_character())) %>%
  select(-label, -numOtus)
sample <- tooll.otu.single$Group
tooll.otu.single <- select(tooll.otu.single, -Group)
tooll.otu.single <- tooll.otu.single[,colSums(tooll.otu.single)>1]
tooll.otu.single <- cbind(tooll.otu.single, sample)

rel_abun_rowSums <- function(x, y) {
  mutate(y, n=rowSums(x))
}

total_reads_single <- tooll.otu.single %>% 
  group_by(sample) %>% 
  group_map(rel_abun_rowSums) %>%
  bind_rows()

tooll.otu.single <- tooll.otu.single %>% mutate(reads=total_reads_single$n)

tooll.otu.single <- tooll.otu.single %>%
  pivot_longer(cols = c(-sample, -reads), names_to = "otu", values_to = "count") %>%
  mutate(rel_abund=(count/reads*100))

# @input reformatted OTU table and OTU taxonomy table
# @return % rel abund and taxonomy for every OTU from every sample
tl.otu.tax <- inner_join(tooll.otu.single, tooll.tax)

# @input metadata for all samples and % rel abund and taxonomy for every OTU from every sample
# @return complete dataframe of all OTUs from every sample with everything we want to know about each OTU
tl.otu.tax.meta <- inner_join(tl.otu.tax, metadatat, by=c("sample"="group"))

#Functional groups for Figure
feox<-filter(tl.otu.tax.meta, genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella")) %>%
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarize(rel_abund_sum=sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

ferb<-filter(tl.otu.tax.meta, genus%in%c("Geobacter", "Geobacteraceae_uncultured", "Geotalea", "Citrifermentans", "Shewanella", 
                                      "Anaeromyxobacter", "Pelobacter", "Desulfuromonas", "Carboxydothermus", "Geothrix")) %>%
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarise(rel_abund_sum = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

# used Dedysh, S. N., & Knief, C. (2018). Diversity and Phylogeny of Described Aerobic Methanotrophs. 
# Methane Biocatalysis: Paving the Way to Sustainability, 17–42. and
# ISME Journal (2019) 13:1209–1225 that shows Methylomonadaceae is a methanotroph
methyl<-filter(tl.otu.tax.meta, family%in%c("Methylobacteriaceae", "Methylophilaceae", "Methylomonadaceae", "Methylococcaceae", "Methylothermaceae", 
                                         "Methylocystaceae", "Beijerinckiaceae", "Methylacidophilaceae")) %>% 
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarise(rel_abund_sum = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

# using Holmes and Smith 2016 Advances in Applied Microbiology section 2.2 Phylogeny of Methanogens to determine methanogens families
methano<-filter(tl.otu.tax.meta, order%in%c("Methanobacteriales", "Methanococcales", "Methanomicrobiales", "Methanosarcinales", "Methanopyrales")) %>% 
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarise(rel_abund_sum = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

#Summaries of other relevant taxa discussed in the text,
#not called out specifically in the figures
tl.gal<- filter(tl.otu.tax.meta, family%in%c("Gallionellaceae")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

tl.geo<-filter(tl.otu.tax.meta, family%in%c("Geobacteraceae")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

tl.lepto<-filter(tl.otu.tax.meta, genus%in%c("Leptothrix")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

tl.methylpren<-filter(tl.otu.tax.meta, family%in%c("Methanoperedenaceae")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

tl.geothrix<-filter(tl.otu.tax.meta, genus%in%c("Geothrix")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

geochem <- read_excel(path="TLC.xlsx")
fe2 <- geochem$Fe2
o2 <- geochem$o2
geochem_depth <- geochem$depth
depth <- c(0, -0.25, -1.25, -2.25, -3.25, -4.25, -5.25, -6.25)
sfe_depth <- c(-0.25, -1.25, -2.25, -3.25, -4.25, -5.25, -6.25)
sfe <- c(63.6, 14.2, 4.3, 3.4, 0.4, 7.1, 5.7)

geochemistry <- data.frame(geochem_depth, o2, fe2)

#png("TL_stnC.eps", height = 1600, width = 2400, units = 'px', res=300)
#Figure 3 for Hudson et al ESPI
pdf("TL_stnC.pdf", height = 6, width = 6.5)
par(mfrow=c(1,2), pin=c(5,3), mar = c(4,4,4,0.25))
plot(ferb$rel_abund_sum, ferb$depth, col="darkred", pch=15, ylim = c(-8,1.01), xlim = c(0,2), axes=FALSE, ann=FALSE)
lines(ferb$rel_abund_sum, ferb$depth, col="darkred", lty=3, lwd=2)

points(feox$rel_abund_sum , feox$depth, col="darkorange3", pch=19)
lines(feox$rel_abund_sum , feox$depth, col="darkorange3", lty=3, lwd=2)

points(methyl$rel_abund_sum, methyl$depth, col="blue2", pch=17)
lines(methyl$rel_abund_sum, methyl$depth, col="blue2", lty=3, lwd=2)

points(methano$rel_abund_sum, methano$depth, col="black", pch=18, cex=1.25)
lines(methano$rel_abund_sum, methano$depth, col="black", lty=3, lwd=2)

#points(methanopren, depth, col="darkblue", pch=15)
#lines(methanopren, depth, col="darkblue", lty=3, lwd=2)

axis(3)
axis(2, las=2)

mtext("Relative Abundance (%)", side=3, line = 2.5)
mtext("Depth (cm)", side = 2, line = 2.5, las=0)

legend(0.25, -6.3, legend = c("Iron reducing", "Iron oxidizing", "Methane oxidizing", "Methanogenic"), 
       col = c("darkred", "darkorange3", "blue2", "black"), pch=c(15,19,17,18), lty=c(3,3,3,3), 
       lwd=c(2,2,2,2), ncol=1, bty="n")

plot(o2, geochem_depth, col="blue2", pch=15, ylim = c(-8,1.01), xlim = c(0, 310), axes=FALSE, ann=FALSE)
lines(o2, geochem_depth, col="blue2", lty=3, lwd=2)
xtick<-seq(0,300, by=100)
axis(3, at=xtick)
axis(2, las=2)
mtext(expression(paste("Concentration (",mu,"M)")), side=3, line = 2.5)

points(fe2, geochem_depth, col="darkgreen", pch=17)
lines(fe2, geochem_depth, col="darkgreen", lty=3, lwd=2)
par(new = TRUE)
plot(sfe, sfe_depth, col="darkorange3", pch=19, xlim = c(0, 70), ylim = c(-8, 1.01), xaxt = "n", yaxt = "n", 
     ylab = "", xlab = "", axes=FALSE, ann=FALSE)
lines(sfe, sfe_depth, col="darkorange3", lty=3, lwd=2)
x2tick<- seq(0,80, by=20)
axis(side = 1, at=x2tick)
fe_xlab=expression(paste("Fe(III)-Oxides (",mu,"mole/gdw)"))
mtext(fe_xlab, side = 1, line = 2.5)
legend_fe = expression(paste("Fe(II)"))
legend(10, -6.3, legend = c("Oxygen", legend_fe, "Fe(III)-Oxides"), col = c("blue2", "darkgreen", "darkorange3"), 
       pch = c(15, 17, 19), lty = c(3,3), lwd = c(2,2), ncol = 1, bty = "n")

dev.off()

## Seep site processing
metadatas<- read_excel(path="TLS_metadata.xlsx")
tools.tax <- read_tsv("tools.tax") %>% 
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string = taxonomy, pattern = ";$", replacement = "")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")

# @input mothur formatted and name shortened shared file (# of OTUS in each sample)
# @return tidy dataframe in long format with OTUs in % relative abundance
tools.otu.single <- read_tsv("tools.shared", col_types = cols(Group=col_character())) %>%
  select(-label, -numOtus)
sample <- tools.otu.single$Group
tools.otu.single <- select(tools.otu.single, -Group)
tools.otu.single <- tools.otu.single[,colSums(tools.otu.single)>1]
tools.otu.single <- cbind(tools.otu.single, sample)

rel_abun_rowSums <- function(x, y) {
  mutate(y, n=rowSums(x))
}

total_reads_single <- tools.otu.single %>% 
  group_by(sample) %>% 
  group_map(rel_abun_rowSums) %>%
  bind_rows()

tools.otu.single <- tools.otu.single %>% mutate(reads=total_reads_single$n)

tools.otu.single <- tools.otu.single %>%
  pivot_longer(cols = c(-sample, -reads), names_to = "otu", values_to = "count") %>%
  mutate(rel_abund=(count/reads*100))

tls.otu.tax <- inner_join(tools.otu.single, tools.tax)

# @input metadata for all samples and % rel abund and taxonomy for every OTU from every sample
# @return complete dataframe of all OTUs from every sample with everything we want to know about each OTU
tls.otu.tax.meta <- inner_join(tls.otu.tax, metadatas, by=c("sample"="group"))

#Functional groups for Figure
feox.s<-filter(tls.otu.tax.meta, genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella")) %>%
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarize(rel_abund_sum=sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

ferb.s<-filter(tls.otu.tax.meta, genus%in%c("Geobacter", "Geobacteraceae_uncultured", "Geotalea", "Citrifermentans", "Shewanella", 
                                         "Anaeromyxobacter", "Pelobacter", "Desulfuromonas", "Carboxydothermus", "Geothrix")) %>%
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarise(rel_abund_sum = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

# used Dedysh, S. N., & Knief, C. (2018). Diversity and Phylogeny of Described Aerobic Methanotrophs. 
# Methane Biocatalysis: Paving the Way to Sustainability, 17–42. and
# ISME Journal (2019) 13:1209–1225 that shows Methylomonadaceae is a methanotroph
methyl.s<-filter(tls.otu.tax.meta, family%in%c("Methylobacteriaceae", "Methylophilaceae", "Methylomonadaceae", "Methylococcaceae", "Methylothermaceae", 
                                            "Methylocystaceae", "Beijerinckiaceae", "Methylacidophilaceae")) %>% 
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarise(rel_abund_sum = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

# using Holmes and Smith 2016 Advances in Applied Microbiology section 2.2 Phylogeny of Methanogens to determine methanogens families
methano.s<-filter(tls.otu.tax.meta, order%in%c("Methanobacteriales", "Methanococcales", "Methanomicrobiales", "Methanosarcinales", "Methanopyrales")) %>% 
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarise(rel_abund_sum = sum(rel_abund)) %>%
  inner_join(., metadatat, by=c("sample" = "group"))

#Other taxa discussed in the text, but not detailed in the Figure
tls.gal<- filter(tls.otu.tax.meta, family%in%c("Gallionellaceae")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

tls.geo<-filter(tls.otu.tax.meta, family%in%c("Geobacteraceae")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

tls.lepto<-filter(tls.otu.tax.meta, genus%in%c("Leptothrix")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

tls.methylpren<-filter(tls.otu.tax.meta, family%in%c("Methanoperedenaceae")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

# Using the 5 established orders of methanogens as stated in Liu 2010 Taxonomy of Methanogens in Handbook of Hydrocarbon and Lipid microbiology
tls.methan<-filter(tls.otu.tax.meta, order%in%c("Methanobacteriales", "Methanococcales", "Methanomicrobiales", "Methanosarcinales", "Methanopyrales")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

tls.geothrix<-filter(tls.otu.tax.meta, genus%in%c("Geothrix")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

tls.rhodo<-filter(tls.otu.tax.meta, genus%in%c("Rhodoferax")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadatas, by=c("sample" = "group"))

ferb_seep <- tls.geo$rel_abund + tls.geothrix$rel_abund
feob_seep <- tls.gal$rel_abund + tls.lepto$rel_abund
methyl_seep <- tls.methyl$rel_abund
methan_seep <- tls.methan$rel_abund
methylpren_seep <- tls.methylpren$rel_abund
rhodo_seep <- tls.rhodo$rel_abund

microbe<-c('Iron Oxidizing', 'Iron Reducing', 'Methane Oxidizing', 'Methanogenic')
seep_relabund<-c(feox.s$rel_abund_sum, ferb.s$rel_abund_sum, methyl.s$rel_abund_sum, methano.s$rel_abund_sum)

seep_relabund<-data.frame(microbe, seep_relabund)

seep <- ggplot(seep_relabund, aes(x=microbe, y=seep_relabund)) + 
  geom_bar(stat = "identity") +
  ylab("Relative Abundance (%)") +
  xlab("Functional Groups of Microorganisms") +
  coord_flip() +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("Figure4B.pdf")
print(seep)
dev.off()