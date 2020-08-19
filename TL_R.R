library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)

# Significant taxa
metadata<- read_excel(path="TL_metadata.xlsx")
tooll.tax <- read_tsv("tooll.otu.taxonomy") %>% 
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string = taxonomy, pattern = ";$", replacement = "")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")

# @input mothur formatted and name shortened shared file (# of OTUS in each sample)
# @return tidy dataframe in long format with OTUs in % relative abundance
tooll.otu <- read_tsv("tooll.otu.table", col_types = cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(sample=Group) %>%
  pivot_longer(cols = -sample, names_to = "otu", values_to = "count") %>%
  mutate(rel_abund=(count/18010)*100)

# @input reformatted OTU table and OTU taxonomy table
# @return % rel abund and taxonomy for every OTU from every sample
tl.otu.tax <- inner_join(tooll.otu, tooll.tax)

# @input metadata for all samples and % rel abund and taxonomy for every OTU from every sample
# @return complete dataframe of all OTUs from every sample with everything we want to know about each OTU
tl.otu.tax.meta <- inner_join(tl.otu.tax, metadata, by=c("sample"="group")) %>%
  filter(., sample!="TL1")

tl.gal<- filter(tl.otu.tax.meta, family%in%c("Gallionellaceae")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

tl.geo<-filter(tl.otu.tax.meta, family%in%c("Geobacteraceae")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

tl.lepto<-filter(tl.otu.tax.meta, genus%in%c("Leptothrix")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))
#There is hardly any leptothrix (<0.005% rel abund) do not use in final analysis, will not impact FeOB rel abund

# used Dedysh, S. N., & Knief, C. (2018). Diversity and Phylogeny of Described Aerobic Methanotrophs. 
# Methane Biocatalysis: Paving the Way to Sustainability, 17–42. and
# ISME Journal (2019) 13:1209–1225 that shows Methylomonadaceae is a methanotroph
tl.methyl<-filter(tl.otu.tax.meta, family%in%c("Methylomonadaceae", "Methylococcaceae")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

tl.methylpren<-filter(tl.otu.tax.meta, family%in%c("Methanoperedenaceae")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

# using Holmes and Smith 2016 Advances in Applied Microbiology section 2.2 Phylogeny of Methanogens to determine methanogens families
tl.methan<-filter(tl.otu.tax.meta, family==c("Methanoregulaceae", "Methanosaetaceae", "Methanobacteriaceae")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

tl.geothrix<-filter(tl.otu.tax.meta, genus%in%c("Geothrix")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

tl.cyano<-filter(tl.otu.tax.meta, class%in%c("Cyanobacteriia")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))


geochem <- read_excel(path="TLC.xlsx")
fe2 <- geochem$Fe2
o2 <- geochem$o2
geochem_depth <- geochem$depth
depth <- c(0, -0.25, -1.25, -2.25, -3.25, -4.25)
sfe_depth <- c(-0.25, -1.25, -2.25, -3.25, -4.25, -5.25, -6.25)
sfe <- c(63.6, 14.2, 4.3, 3.4, 0.4, 7.1, 5.7)

## plot multiple dots and lines on same plot, depth profile of OTUs of interest in this case
geo <- tl.geo$rel_abund + tl.geothrix$rel_abund
gal <- tl.gal$rel_abund
methyl <- tl.methyl$rel_abund
methan <- tl.methan$rel_abund
methanopren <- tl.methylpren$rel_abund


#png("TL_stnC.eps", height = 1600, width = 2400, units = 'px', res=300)
pdf("TL_stnC.pdf", height = 6, width = 6.5)

par(mfrow=c(1,2), pin=c(5,3), mar = c(4,4,4,0.25))
plot(geo, depth, col="darkred", pch=15, ylim = c(-6,1.05), xlim = c(0,2.5), axes=FALSE, ann=FALSE)
axis(3)
axis(2, las=2)
lines(geo, depth, col="darkred", lty=3, lwd=2)
mtext("Relative Abundance (%)", side=3, line = 2.5)
mtext("Depth (cm)", side = 2, line = 2.5, las=0)

points(gal, depth, col="darkorange2", pch=16)
lines(gal, depth, col="darkorange2", lty=3, lwd=2)

points(methyl, depth, col="blue2", pch=17)
lines(methyl, depth, col="blue2", lty=3, lwd=2)

points(methan, depth, col="black", pch=18, cex=1.25)
lines(methan, depth, col="black", lty=3, lwd=2)

#points(methanopren, depth, col="darkblue", pch=15)
#lines(methanopren, depth, col="darkblue", lty=3, lwd=2)

legend(0.5, -4.5, legend = c("Iron reducing", "Iron oxidizing", "Methane oxidizing", "Methanogenic"), 
       col = c("darkred", "darkorange2", "blue2", "black"), pch=c(15,16,17,18), lty=c(3,3,3,3), 
       lwd=c(2,2,2,2), ncol=1, bty="n")

plot(o2, geochem_depth, col="blue2", pch=17, ylim = c(-6,1.05), xlim = c(0, 310), axes=FALSE, ann=FALSE)
lines(o2, geochem_depth, col="blue2", lty=3, lwd=2)
xtick<-seq(0,300, by=100)
axis(3, at=xtick)
axis(2, las=2)
mtext(expression(paste("Concentration (",mu,"M)")), side=3, line = 2.5)

points(fe2, geochem_depth, col="darkred", pch=15)
lines(fe2, geochem_depth, col="darkred", lty=3, lwd=2)
par(new = TRUE)
plot(sfe, sfe_depth, col="orange", pch=16, xlim = c(0, 70), ylim = c(-6, 1.05), xaxt = "n", yaxt = "n", 
     ylab = "", xlab = "", axes=FALSE, ann=FALSE)
lines(sfe, sfe_depth, col="orange", lty=3, lwd=2)
axis(side = 1)
fe_xlab=expression(paste("Solid phase Fe(III) (",mu,"mole/gdw)"))
mtext(fe_xlab, side = 1, line = 2.5)
legend_fe = expression(paste("Fe"^"2+"))
legend(20, -5, legend = c("Oxygen", legend_fe, "Fe(III)"), col = c("blue2", "darkred", "orange"), 
       pch = c(17, 15, 16), lty = c(3,3), lwd = c(2,2), ncol = 1, bty = "n")

dev.off()
