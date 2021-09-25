# Load necessary libraries
library(readxl)
library(ggplot2)
library(colortools)
library(ggrepel)
library(reshape2)
library(rcompanion)
library(lmPerm)
library(RVAideMemoire)
library(labdsv)
library(metacom)
library(BiodiversityR)
library(tidyverse)

# Load necessary data files, parse data
RichEvenSamp         <- read_excel("RichEvenSamp.xlsx")
Forage_BasalResource <- read_excel("Forage_BasalResource.xlsx")
Forage_Permanova     <- read_csv("Forage_Permanova.csv")
RichEvenSamp$nseines <- Forage_Permanova$nseines
RichEvenSamp$ntrawls <- Forage_Permanova$ntrawls
RichEvenSamp$Bay     <- as.factor(RichEvenSamp$Bay)
RichEvenSamp$Year    <- as.factor(RichEvenSamp$Year)

# Calculate average abundance by estuary, transform, and assign column names
Forage_Permanova_avg           <- recast(Forage_Permanova, formula = Bay ~ variable, measure.var = c(9:58), fun.aggregate = mean)
Forage_Permanova_4rt           <- (Forage_Permanova_avg[, 2:51])^0.25
rownames(Forage_Permanova_4rt) <- Forage_Permanova_avg$Bay
colnames(Forage_Permanova_4rt) <- c("Alabama Shad",
                                    "Skipjack Herring",
                                    "Cuban Anchovy",
                                    "Striped Anchovy",
                                    "Dusky Anchovy",
                                    "Bay Anchovy",
                                    "Anchovies",
                                    "Silver Perch",
                                    "Menhadens",
                                    "Blue Runner",
                                    "Crevalle Jack",
                                    "Horse-eye Jack",
                                    "Atlantic Bumper",
                                    "Herrings",
                                    "Sheepshead Minnow",
                                    "Round Scad",
                                    "Irish Pompano",
                                    "Gizzard Shad",
                                    "Threadfin Shad",
                                    "Eucinostomus Mojarras",
                                    "Marsh Killifish",
                                    "Gulf Killifish",
                                    "Seminole Killifish",
                                    "Longnose Killifish",
                                    "Eastern Mosquitofish",
                                    "Tomtate",
                                    "Scaled Sardine",
                                    "Reef Silverside",
                                    "False Silver Halfbeak",
                                    "Halfbeak",
                                    "Halfbeaks",
                                    "Pinfish",
                                    "Spot",
                                    "Bluefin Killifish",
                                    "Rainwater Killifish",
                                    "Rough Silverside",
                                    "Menidia Silversides",
                                    "Striped Mullet",
                                    "White Mullet",
                                    "Fantail Mullet",
                                    "Mullets",
                                    "Leatherjacket",
                                    "Atlantic Thread Herring",
                                    "Pigfish",
                                    "Gulf Butterfish",
                                    "Harvestfish",
                                    "Sailfin Molly",
                                    "Spanish Sardine",
                                    "Atlantic Moonfish",
                                    "Longspine Porgy")

# Order matrix by reciprocal averaging for shadeplot
Forage_order   <- OrderMatrix(Forage_Permanova_4rt, binary = FALSE, outputScores = TRUE)
taxID          <- sort(Forage_order$speciesscores, decreasing = TRUE, index.return = TRUE)$ix
Forage_ordered <- Forage_Permanova_4rt[order(Forage_order$sitescores, decreasing = TRUE), taxID]

# Assign numerical labels to taxa
for (i in 1:ncol(Forage_ordered)) {
    if (i < 10) {
        colnames(Forage_ordered)[i] <- paste(colnames(Forage_ordered)[i], "   (", i, ")", sep = "")
    } else {
        colnames(Forage_ordered)[i] <- paste(colnames(Forage_ordered)[i], " (", i, ")", sep = "")
    }
}

# Assign grouping factor for estuaries
Forage_ordered$bay <- factor(rownames(Forage_ordered), levels = c("AP", "CK", "TB", "CH"))


# Wide format -> long format
melt_forage                                  <- melt(Forage_ordered)
melt_forage[melt_forage$value == 0, "value"] <- NA

# Produce shadeplot
forage_shade <- ggplot(melt_forage, aes(bay, reorder(variable,desc(variable)))) +
  geom_tile(aes(fill = value),show.legend = FALSE) +
  scale_fill_gradient(low = "gray90",high = "black",na.value = "white")+

  theme_void() +
  scale_x_discrete(position = "top") +
  theme(axis.text = element_text(hjust = 1),axis.text.x = element_text(hjust = 0.5,size = 24))
ggsave("Forage_ShadePlot.tiff",forage_shade,width = 6,height = 9,dpi = 300)

# Prepare data for rank abundance and frequency plots
rownames(Forage_BasalResource) <- Forage_BasalResource$Taxa
basal                          <- data.frame(Forage_BasalResource$Basal)
basal$tax                      <- rownames(Forage_BasalResource)
basal                          <- data.frame(basal[taxID, ])
rownames(basal)                <- as.character(1:nrow(basal))
basal                          <- data.frame(basal = basal[, 1])
Forage_tax                     <- Forage_Permanova[, 9:58]
Forage_tax                     <- Forage_tax[, taxID]
colnames(Forage_tax)           <- as.character(1:ncol(Forage_tax))
Forage_fact                    <- Forage_Permanova[, 2:6]
Forage_fact$Bay                <- factor(Forage_fact$Bay, levels = c("AP", "CK", "TB", "CH"))
Forage_fact$Zone               <- as.factor(Forage_fact$Zone)
Forage_fact$Season             <- as.factor(Forage_fact$Season)
Forage_Permanova$Bay           <- as.factor(Forage_Permanova$Bay)


# Calculate rank abundance by estuary
forage_rankabun_AB        <- rankabundance(1 * Forage_tax, y = as.data.frame(Forage_Permanova), factor = "Bay", level = "AP")
forage_rankabun_AB        <- data.frame(bay = "AB", rank = forage_rankabun_AB[, "rank"], pctn = forage_rankabun_AB[, "proportion"])
forage_rankabun_AB        <- merge(forage_rankabun_AB, basal, by = 0)
forage_rankabun_AB        <- forage_rankabun_AB[order(forage_rankabun_AB$rank), ]
forage_rankabun_AB$cumul  <- cumsum(forage_rankabun_AB$pctn)

forage_rankabun_CK        <- rankabundance(1 * Forage_tax, y = as.data.frame(Forage_Permanova), factor = "Bay", level = "CK")
forage_rankabun_CK        <- data.frame(bay = "CK", rank = forage_rankabun_CK[, "rank"], pctn = forage_rankabun_CK[, "proportion"])
forage_rankabun_CK        <- merge(forage_rankabun_CK, basal, by = 0)
forage_rankabun_CK        <- forage_rankabun_CK[order(forage_rankabun_CK$rank), ]
forage_rankabun_CK$cumul  <- cumsum(forage_rankabun_CK$pctn)

forage_rankabun_TB        <- rankabundance(1 * Forage_tax, y = as.data.frame(Forage_Permanova), factor = "Bay", level = "TB")
forage_rankabun_TB        <- data.frame(bay = "TB", rank = forage_rankabun_TB[, "rank"], pctn = forage_rankabun_TB[, "proportion"])
forage_rankabun_TB        <- merge(forage_rankabun_TB, basal, by = 0)
forage_rankabun_TB        <- forage_rankabun_TB[order(forage_rankabun_TB$rank), ]
forage_rankabun_TB$cumul  <- cumsum(forage_rankabun_TB$pctn)

forage_rankabun_CH        <- rankabundance(1 * Forage_tax, y = as.data.frame(Forage_Permanova), factor = "Bay", level = "CH")
forage_rankabun_CH        <- data.frame(bay = "CH", rank = forage_rankabun_CH[, "rank"], pctn = forage_rankabun_CH[, "proportion"])
forage_rankabun_CH        <- merge(forage_rankabun_CH, basal, by = 0)
forage_rankabun_CH        <- forage_rankabun_CH[order(forage_rankabun_CH$rank), ]
forage_rankabun_CH$cumul  <- cumsum(forage_rankabun_CH$pctn)

# Parse all rank abundances into one df for plotting
forage_rankabun           <- rbind(forage_rankabun_AB, forage_rankabun_CH, forage_rankabun_CK, forage_rankabun_TB)
colnames(forage_rankabun) <- c("taxa", "bay", "rank", "pctn", "basal", "cumul")
forage_rankabun$bay       <- as.factor(forage_rankabun$bay)
forage_rankabun           <- forage_rankabun[order(forage_rankabun$rank), ]
rankabun_plots            <- list()

# Construct rank abundance plots by estuary
for (i in levels(forage_rankabun$bay)){
  abun_plot <- ggplot(filter(forage_rankabun, bay == i), aes(x = rank, y = pctn,fill = basal,label = taxa)) +
    geom_col() +
    theme_classic() +
    scale_fill_manual(values = c("darkgreen","darkblue","chocolate4","aquamarine3"),limits = c("Benthic","Planktonic","Detritivore","Mixed")) +
    geom_text(vjust = "bottom",nudge_y = 0.5,size = 3,fontface = "bold") +
    labs(x = "Species rank",y = "Relative abundance",fill = "Basal resource") +
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
    scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0,0.1))) +
    theme()
  plot(abun_plot)
  rankabun_plots[[i]] <- abun_plot
}

# Calculate rank frequency as presence-absence rank abundance by estuary
forage_rankfreq_AB <- rankabundance(1 * (Forage_tax > 0), y = as.data.frame(Forage_Permanova), factor = "Bay", level = "AP")
forage_rankfreq_AB <- data.frame(bay = "AB", rank = forage_rankfreq_AB[, "rank"], f = forage_rankfreq_AB[, "abundance"]/nrow(Forage_tax))
forage_rankfreq_AB <- merge(forage_rankfreq_AB, basal, by = 0)
forage_rankfreq_CK <- rankabundance(1 * (Forage_tax > 0), y = as.data.frame(Forage_Permanova), factor = "Bay", level = "CK")
forage_rankfreq_CK <- data.frame(bay = "CK", rank = forage_rankfreq_CK[, "rank"], f = forage_rankfreq_CK[, "abundance"]/nrow(Forage_tax))
forage_rankfreq_CK <- merge(forage_rankfreq_CK, basal, by = 0)
forage_rankfreq_TB <- rankabundance(1 * (Forage_tax > 0), y = as.data.frame(Forage_Permanova), factor = "Bay", level = "TB")
forage_rankfreq_TB <- data.frame(bay = "TB", rank = forage_rankfreq_TB[, "rank"], f = forage_rankfreq_TB[, "abundance"]/nrow(Forage_tax))
forage_rankfreq_TB <- merge(forage_rankfreq_TB, basal, by = 0)
forage_rankfreq_CH <- rankabundance(1 * (Forage_tax > 0), y = as.data.frame(Forage_Permanova), factor = "Bay", level = "CH")
forage_rankfreq_CH <- data.frame(bay = "CH", rank = forage_rankfreq_CH[, "rank"], f = forage_rankfreq_CH[, "abundance"]/nrow(Forage_tax))
forage_rankfreq_CH <- merge(forage_rankfreq_CH, basal, by = 0)

# Parse all rank frequencies into one df for plotting
forage_rankfreq           <- rbind(forage_rankfreq_AB, forage_rankfreq_CH, forage_rankfreq_CK, forage_rankfreq_TB)
colnames(forage_rankfreq) <- c("taxa", "bay", "rank", "freq", "basal")
forage_rankfreq$bay       <- as.factor(forage_rankfreq$bay)
rankfreq_plots            <- list()

# Construct rank frequency plots by estuary
for (i in levels(forage_rankfreq$bay)){
  freq_plot <- ggplot(filter(forage_rankfreq, bay == i), aes(x = rank, y = freq,fill = basal,label = taxa)) +
    geom_col() +
    scale_fill_manual(values = c("darkgreen","darkblue","chocolate4","aquamarine3"),limits = c("Benthic","Planktonic","Detritivore","Mixed")) +
    geom_text(vjust = "bottom",nudge_y = 0.005,size = 3,fontface = "bold") +
    theme_classic() +
    labs(x = "Species rank",y = "Frequency of occurrence",fill = "Basal resource") +
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
    scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0,.1))) +
    theme()
  plot(freq_plot)
  rankfreq_plots[[i]] <- freq_plot
}

# Detrend richness and evenness for sample size
detrendMod_Rich                                          <- lm(S ~ nseines + ntrawls, data = RichEvenSamp)
anova(detrendMod_Rich)
S_detrend                                                <- detrendMod_Rich$residuals + detrendMod_Rich$coefficients[1]
RichEvenSamp$S_detrend                                   <- S_detrend

detrendMod_Even                                          <- lm(`J'` ~ nseines + ntrawls, data = RichEvenSamp)
summary(detrendMod_Even)
J_detrend                                                <- detrendMod_Even$residuals + detrendMod_Even$coefficients[1]
RichEvenSamp[is.na(RichEvenSamp$`J'`) == 0, "J_detrend"] <- J_detrend

# Perform anovas for richness and evenness by bay and season
detrendRich_perm_bay       <- perm.anova(S_detrend ~ Bay, data = RichEvenSamp, nperm = 999)
View(detrendRich_perm_bay)
detrendRich_pair_bay       <- pairwisePermutationTest(S_detrend ~ Bay, data = RichEvenSamp)
detrendRich_PM_bay         <- pairwisePermutationMatrix(S_detrend ~ Bay, data = RichEvenSamp)

detrendRich_perm_Season    <- perm.anova(S_detrend ~ as.factor(Season), data = RichEvenSamp, nperm = 999)
View(detrendRich_perm_Season)
detrendRich_pair_Season    <- pairwisePermutationTest(S_detrend ~ as.factor(Season), data = RichEvenSamp)
detrendRich_PM_Season      <- pairwisePermutationMatrix(S_detrend ~ as.factor(Season), data = RichEvenSamp)


detrendEven_perm_bay       <- perm.anova(J_detrend ~ Bay, data = RichEvenSamp, nperm = 999) # Not significant
View(detrendEven_perm_bay)

detrendEven_perm_Season    <- perm.anova(J_detrend ~ as.factor(Season), data = RichEvenSamp, nperm = 999)
View(detrendEven_perm_Season)
detrendEven_pair_Season    <- pairwisePermutationTest(J_detrend ~ as.factor(Season), data = RichEvenSamp)
detrendEven_PM_Season      <- pairwisePermutationMatrix(J_detrend ~ as.factor(Season), data = RichEvenSamp)

# Construct plots for richness and evenness by bay and season
RichBay_plot    <- ggplot(RichEvenSamp, aes(x=factor(Bay,levels = c("AP","CK","TB","CH")), y=S_detrend)) + 
  theme_classic() +
  scale_x_discrete(labels = c("AB","CK","TB","CH")) +
  labs(x = NULL,y = "Species richness") +
  geom_boxplot(size = 1,outlier.size = 1)
ggsave("RichBay_plot.tiff",RichBay_plot,width = 4,height = 4,dpi = 300)

RichSeason_plot <- ggplot(RichEvenSamp, aes(x=as.factor(Season), y=S_detrend)) + 
  theme_classic() +
  scale_x_discrete(labels  = c("Wi","Sp","Su","Fa")) +
  labs(x = NULL,y = "Species richness") +
  geom_boxplot(size = 1,outlier.size = 1)
ggsave("RichSeason_plot.tiff",RichSeason_plot,width = 4,height = 4,dpi = 300)

EvenBay_plot    <- ggplot(RichEvenSamp, aes(x=factor(Bay,levels = c("AP","CK","TB","CH")), y=J_detrend)) + 
  theme_classic() +
  scale_x_discrete(labels = c("AB","CK","TB","CH")) +
  labs(x = NULL,y = "Species evenness") +
  geom_boxplot(size = 1,outlier.size = 1)
ggsave("EvenBay_plot.tiff",EvenBay_plot,width = 4,height = 4,dpi = 300)

EvenSeason_plot <- ggplot(RichEvenSamp, aes(x=as.factor(Season), y=J_detrend)) + 
  theme_classic() +
  scale_x_discrete(labels  = c("Wi","Sp","Su","Fa")) +
  labs(x = NULL,y = "Species evenness") +
  geom_boxplot(size = 1,outlier.size = 1)
ggsave("EvenSeason_plot.tiff",EvenSeason_plot,width = 4,height = 4,dpi = 300)

# Conduct modified indval analysis for all samples
clust_all          <- rep(1, nrow(Forage_tax))
forage_indval_all  <- indval(Forage_tax, clust_all)

avg_abund          <- apply(Forage_tax, 2, mean)
rel_abund          <- avg_abund/sum(avg_abund)
indvals_all        <- forage_indval_all$relfrq * rel_abund
indvals_all[, 2]   <- forage_indval_all$relfrq
indvals_all[, 3]   <- rel_abund

# Conduct modified indval analysis by estuary
indvals_bay           <- list()
for (i in levels(Forage_fact$Bay)) {
    clust_bay         <- rep(1, nrow(Forage_tax[Forage_fact$Bay == i, ]))
    forage_indval_bay <- indval(Forage_tax[Forage_fact$Bay == i, colSums(Forage_tax[Forage_fact$Bay == i, ]) > 0], clust_bay)
    
    avg_abund         <- apply(Forage_tax[Forage_fact$Bay == i, colSums(Forage_tax[Forage_fact$Bay == i, ]) > 0], 2, mean)
    rel_abund         <- avg_abund/sum(avg_abund)
    indvals_bay[[i]]  <- data.frame(forage_indval_bay$relfrq * rel_abund, rel_abund, forage_indval_bay$relfrq)
}


# Save plot and numerical id data
save(forage_shade, file = "ShadePlot.RData")

save(RichBay_plot, RichSeason_plot, EvenBay_plot, EvenSeason_plot, file = "SummStatPlots.RData")

save(rankfreq_plots, rankabun_plots, file = "RankPlots.RData")

save(taxID, file = "TaxaID.RData")
