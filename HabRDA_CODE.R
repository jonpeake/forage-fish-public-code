# Load Required Packages
library(readxl)
library(vegan)
library(labdsv)
library(ggplot2)
library(colortools)
library(ggrepel)
library(reshape2)
library(readr)
library(RColorBrewer)
library(egg)
library(tidyverse)

# Import data
ForageAbioticComps <- read_csv("Forage_AbioticComps.csv")

# Parse data
taxa <- ForageAbioticComps[, 9:58]
samp <- ForageAbioticComps[, 7:8]
abio <- ForageAbioticComps[, 59:64]
fact <- ForageAbioticComps[, 2:6]
vars <- data.frame(abio, samp)

# Calculate vif to test for multicollinearity
vif  <- diag(solve(cor(vars)))
1 - (cor(vars)^2)

# Clean workspace
rm(ForageAbioticComps)

# Set up null dbRDA model with conditions for variable selection on square root bray curtis of 4th root transform
modNull <- dbrda(taxa^0.25 ~ Condition(nseines) + Condition(ntrawls), sqrt.dist = T, dist = "bray", data = vars)

# Conduct stepwise selection of dbrda model
modBest <- step(modNull, scope = c(lower = taxa^0.25 ~ Condition(nseines) + Condition(ntrawls), upper = taxa^0.25 ~ bottomveg + temp + pH + sal + DO + vis + Condition(nseines) + Condition(ntrawls)), 
    direction = "both", trace = 2)

# Run permutation test for full optimal dbRDA model, by axis, and by term
anova_best_Full <- anova(modBest, parallel = 2)
anova_best_Full
anova_best_Axis <- anova(modBest, by = "axis", parallel = 2, model = "reduced", cutoff = 0.05)
anova_best_Axis
anova_best_Term <- anova(modBest, by = "terms", parallel = 2, model = "reduced")
anova_best_Term

# Assign column names of species for plotting
load("TaxaID.RData")
taxa           <- taxa[, taxID]
colnames(taxa) <- as.character(1:ncol(taxa))

# Calculate species scores based on hellinger transform
sppscores(modBest) <- decostand(taxa^0.25, "hellinger")

# Isolate top 15 taxa
clust_all          <- rep(1, nrow(taxa))
indval_all         <- indval(taxa, clust_all)
avg_abund          <- apply(taxa, 2, mean)
rel_abund          <- avg_abund/sum(avg_abund)
indvals_all        <- data.frame(indval = indval_all$relfrq * rel_abund, relfrq = indval_all$relfrq, relabn = rel_abund)
names(indvals_all) <- c("indval", "relfrq", "relabn")
indval_top         <- indvals_all[order(indvals_all$indval, decreasing = TRUE), ]
indval_top15       <- indval_top[1:15, ]
taxa_subset        <- rownames(indval_top15)

# Convert factor columns to factors
fact$Bay <- factor(fact$Bay, levels = c("AP", "CK", "TB", "CH"))
fact$Zone <- as.factor(fact$Zone)
fact$Year <- as.factor(fact$Year)
fact$Season <- as.factor(fact$Season)

# Isolate species scores
scores_spp        <- data.frame(scores(modBest, display = "sp"))
scores_spp_subset <- scores_spp[taxa_subset, ]

# Isolate site scores and calculate 95% ellipse of centroid by season and bay
scores_sites           <- scores(modBest, display = "wa")
scores_sites           <- data.frame(scores_sites, Bay = fact$Bay, Season = fact$Season, BaySeason = interaction(fact$Bay, fact$Season))
ord_ellipse_bay        <- ordiellipse(modBest, scores_sites$Bay, display = "wa", kind = "se", conf = 0.95, draw = "none")
ord_ellipse2_bay       <- ordiellipse(modBest, scores_sites$Bay, display = "wa", kind = "sd", conf = 0.95, draw = "none")
ord_ellipse_season     <- ordiellipse(modBest, scores_sites$Season, display = "wa", kind = "se", conf = 0.95, draw = "none")
ord_ellipse2_season    <- ordiellipse(modBest, scores_sites$Season, display = "wa", kind = "sd", conf = 0.95, draw = "none")
coords_ellipse_bay     <- data.frame()
coords_ellipse_season  <- data.frame()
coords_ellipse2_bay    <- data.frame()
coords_ellipse2_season <- data.frame()

# Produce ellipse path coordinates, groups, and color for ggplot
for (g in levels(scores_sites$Bay)) {
    coords_ellipse_bay <- rbind(coords_ellipse_bay, cbind(as.data.frame(with(scores_sites[scores_sites$Bay, ], vegan:::veganCovEllipse(2.5 * ord_ellipse_bay[[g]]$cov, 2.5 * ord_ellipse_bay[[g]]$center, 
        2.5 * ord_ellipse_bay[[g]]$scale))), Bay = g))
}

for (h in levels(scores_sites$Season)) {
    coords_ellipse_season <- rbind(coords_ellipse_season, cbind(as.data.frame(with(scores_sites[scores_sites$Season, ], vegan:::veganCovEllipse(2.5 * ord_ellipse_season[[h]]$cov, 2.5 * 
        ord_ellipse_season[[h]]$center, 2.5 * ord_ellipse_season[[h]]$scale))), Season = h))
}

coords_ellipse_bay$Bay <- factor(coords_ellipse_bay$Bay,levels = c("AP","CK","TB","CH"))

# Extract environmental vectors for ggplot
scores_envvecs <- data.frame(scores(modBest, display = "bp", scaling = "species"))

# Establish scaling factor for site and habitat scores
mult_vars_hab   <- 7
mult_scores_hab <- 2.5

# Establish color scheme
cbPalette1 <- brewer.pal(4, "Dark2")
cbPalette2 <- brewer.pal(8, "Paired")[c(2, 3, 6, 7)]  # Lighter
cbPalette2 <- brewer.pal(8, "Paired")[c(2, 4, 6, 8)]  # Darker

modBest_hab <- modBest


rdaAll_hab <- ggplot(scores_sites) +
  geom_vline(xintercept = 0,
             colour     = "grey70",
             size       = 1) +
  geom_hline(yintercept = 0,
             colour     = "grey70",
             size       = 1) +
  scale_x_continuous(limits       = symmetric_range,
                     breaks       = c(-6,-3,0,3,6),
                     minor_breaks = NULL) +
  scale_y_continuous(limits       = symmetric_range,
                     breaks       = c(-6,-3,0,3,6),
                     minor_breaks = NULL) +
  labs(title = NULL,
       x     = paste("CA1 (", round(100*modBest_hab$CCA$eig[1]/modBest$tot.chi,2), "%)", sep = ""),
       y     = paste("CA1 (", round(100*modBest_hab$CCA$eig[2]/modBest$tot.chi,2), "%)", sep = ""),
       fill  = NULL,
       color = NULL) +
  scale_color_manual(values = cbPalette1,
                     labels = c("AB","CK","TB","CH")) +
  theme_bw() +
  theme(legend.text       = element_text(size=rel(1)),
        legend.position   = c(0.051,0.75),
        legend.spacing    = unit(-.15,"in"),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        legend.title      = element_blank(),
        panel.grid        = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 1),position = c(0.1,0.5)),
         fill  = guide_legend(override.aes = list(size = 1),position = c(0.2,0.5))) +
  geom_point(aes(x     = mult_scores_hab*dbRDA1,
                 y     = mult_scores_hab*dbRDA2),
             size   = 1,
             stroke = 0.1,
             alpha  = 0.1) +
  scale_fill_manual(values = c(cbPalette2),
                    labels = c("Wi","Sp","Su","Fa")) +
  geom_segment(data = scores_envvecs,
               aes(x    = 0, 
                   xend = mult_vars_hab*dbRDA1, 
                   y    = 0, 
                   yend = mult_vars_hab*dbRDA2),
               arrow = arrow(length = unit(.1,"inches")),
               color = "grey1",
               size  = 1,
               alpha = 0.6) + 
  geom_label_repel(data = scores_spp_subset,
                   aes(x     = dbRDA1,
                       y     = dbRDA2,
                       label = rownames(scores_spp_subset)),
                   force         = 0.01,
                   alpha         = 0.8,
                   size          = 3.5,
                   fontface      = "bold",
                   label.padding = 0.25,
                   colour        = "darkgreen",
                   fill          ="white",
                   label.r       = 0.25,
                   label.size    = 0.25,
                   segment.alpha = 0) +
  geom_label_repel(data  = scores_envvecs,
                   aes(x     = mult_vars_hab*dbRDA1, 
                       y     = mult_vars_hab*dbRDA2, 
                       label = c("Bottom\nVegetation","Temperature","Salinity","Water\nClarity","pH","DO")),
                   segment.alpha = 0,
                   size          = 4,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   box.padding   = .25,
                   label.padding = .25,
                   lineheight    = 0.75,
                   label.r       = 0.25,
                   label.size    = 0.5,
                   nudge_x       = ifelse(scores_envvecs$dbRDA1>0,0.1,-0.1),
                   nudge_y       = ifelse(scores_envvecs$dbRDA2>0,0.2,-0.2)) +
  geom_polygon(data = coords_ellipse_season,
               aes(x     = dbRDA1,
                   y     = dbRDA2,
                   fill  = Season),
               alpha = 0.6,
               size  = 1) +
  geom_polygon(data = coords_ellipse_bay,
             aes(x     = dbRDA1,
                 y     = dbRDA2,
                 color = Bay),
             size  = 1,
             alpha = 0)
ggsave("RDA_plot_all.tiff", rdaAll_hab, width = 8, height = 5, dpi = 300)



save(rdaAll_hab, coords_ellipse_season, coords_ellipse_bay, scores_spp_subset, scores_envvecs, modBest_hab, file = "HabRDA_Plot.RData")
