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

# Import and parse data
load("TaxaID.RData")
PredPreyComps  <- read_csv("PredPreyComps.csv")
fact           <- PredPreyComps[, 2:6]
pred           <- PredPreyComps[, 8:60]
pred_samp      <- PredPreyComps[, 7]
prey           <- PredPreyComps[, 63:112]
prey_samp      <- PredPreyComps[, 61:62]
abio           <- PredPreyComps[, 113:118]
prey           <- prey[, taxID]
colnames(prey) <- as.character(1:ncol(prey))
colnames(pred) <- c("White Catfish",
                    "Ocellated Flounder",
                    "Hardhead Catfish",
                    "Silver Perch",
                    "Blacknose Shark",
                    "Crevalle Jack",
                    "Finetooth Shark",
                    "Bull Shark", 
                    "Blacktip Shark",
                    "Bank Sea Bass",
                    "Rock Sea Baa",
                    "Black Sea Bass",
                    "Common Snook",
                    "Mayan Cichlid",
                    "Sand Seatrout",
                    "Spotted Seatrout",
                    "Ladyfish",
                    "Red Grouper",
                    "Nurse Shark",
                    "Blue Catfish",
                    "Channel Catfish",
                    "Hogfish",
                    "Longnose Gar",
                    "Florida Gar",        
                    "Mutton Snapepr",
                    "Gray Snapper",
                    "Lane Snapepr",
                    "Tarpon",
                    "Largemouth Bass",
                    "Hybrid Bass",
                    "Striped Bass",
                    "Gag",
                    "Lemon Shark",
                    "Yellowtail Snapper",
                    "Gulf Toadfish",
                    "Gulf Flounder",
                    "Bluefish",
                    "Smalltooth Sawfish",
                    "Cobia",
                    "Cownose Ray",
                    "Atlantic Sharpnose Shark",
                    "Red Drum",
                    "Spanish Mackerel",
                    "Lookdown",
                    "Great Barracuda",
                    "Northern Sennet",
                    "Bonnethead",
                    "Atlantic Needlefish",
                    "Redfin Needlefish",
                    "Timucu",
                    "Dusky Flounder",
                    "Inshore Lizardfish",
                    "Southern Hake")

# De-trend predator data by sample size, isolate residual axes
detrend_Pred_mod      <- dbrda(pred^0.25 ~ pred_samp, sqrt.dist = T, dist = "bray")
detrend_scores        <- detrend_Pred_mod$CA$u %*% diag(sqrt(detrend_Pred_mod$CA$eig))
detrend_pctVar        <- detrend_Pred_mod$CA$eig/sum(detrend_Pred_mod$CA$eig)
detrend_cumsum        <- cumsum(detrend_pctVar)
detrend_scoreRetained <- data.frame(detrend_scores[, 1:26])

# Concatenate model vars
vars <- data.frame(abio, prey_samp, detrend_scoreRetained)


# Establish null dbRDA model with environmental vars condition, sqrt Bray Curtis on 4th root transform
modNull <- dbrda(prey^0.25 ~ Condition(nseines + ntrawls + bottomveg + temp + pH + sal + DO + vis), sqrt.dist = T, dist = "bray", data = vars)

# Use stepwise variable selection to produce optimal predator model controlling for environment
modBest <- step(modNull, scope = terms.formula(prey^0.25 ~ . + Condition(nseines + ntrawls + bottomveg + temp + pH + sal + DO + vis), data = vars), direction = "both", trace = 2)
modBest

# Hypothesis test of optimal dbRDA model including by axis and by term
anova_best_Full <- anova(modBest, parallel = 2)
anova_best_Full
anova_best_Axis <- anova(modBest, by = "axis", parallel = 2, model = "reduced", cutoff = 0.05)
anova_best_Axis
anova_best_Term <- anova(modBest, by = "terms", parallel = 2, model = "reduced")
anova_best_Term

# Calculate total variability in pred comm selected in model
predVar_expl <- sum(detrend_pctVar[c(1, 6, 16, 4, 3, 13, 20, 9, 2, 10, 19, 5, 11, 8, 25, 14, 17, 7)])

# Calculate top 15 prey and pred
clust_all               <- rep(1, nrow(prey))
indval_all              <- indval(prey, clust_all)
avg_abund               <- apply(prey, 2, mean)
rel_abund               <- avg_abund/sum(avg_abund)
indvals_all             <- data.frame(indval = indval_all$relfrq * rel_abund, 
                                      relfrq = indval_all$relfrq, 
                                      relabn = rel_abund)
names(indvals_all)      <- c("indval", "relfrq", "relabn")
indval_top              <- indvals_all[order(indvals_all$indval, decreasing = TRUE), ]
indval_top15            <- indval_top[1:15, ]
prey_subset             <- rownames(indval_top15)

clust_all_pred          <- rep(1, nrow(pred))
indval_all_pred         <- indval(pred, clust_all_pred)
avg_abund_pred          <- apply(pred, 2, mean)
rel_abund_pred          <- avg_abund_pred/sum(avg_abund_pred)
indvals_all_pred        <- data.frame(indval = indval_all_pred$relfrq * rel_abund_pred, relfrq = indval_all_pred$relfrq, relabn = rel_abund_pred)
names(indvals_all_pred) <- c("indval", "relfrq", "relabn")
indval_top_pred         <- indvals_all_pred[order(indvals_all_pred$indval, decreasing = TRUE), ]
indval_top15_pred       <- indval_top_pred[1:15, ]
pred_subset             <- rownames(indval_top15_pred)

# Convert factor columns to factors
fact$Bay    <- factor(fact$Bay, levels = c("AP", "CK", "TB", "CH"))
fact$Zone   <- as.factor(fact$Zone)
fact$Year   <- as.factor(fact$Year)
fact$Season <- as.factor(fact$Season)

# Calculate correlation vectors for prey and pred
spp_vecs              <- envfit(modBest, decostand(prey^0.25, "hellinger"))
scores_spp            <- data.frame(scores(spp_vecs, display = "vectors"), sub = "prey")

scores_spp_subset     <- scores_spp[prey_subset, ]

rda_envvecs           <- envfit(modBest, decostand(pred^0.25, "hellinger"))
scores_envvecs        <- data.frame(scores(rda_envvecs, display = "vectors"), sub = "pred")
scores_envvecs_subset <- scores_envvecs[pred_subset, ]

vecs_all              <- rbind(scores_envvecs_subset, scores_spp_subset)

# Extract sample scores
scores_sites <- scores(modBest, display = "wa")
scores_sites <- data.frame(scores_sites, Bay = fact$Bay, Season = fact$Season)
modBest_pred <- modBest

# Set scaling facotrs for pred and prey
mult_vars_pred <- 9
mult_scores_pred <- 1

# Set color scales
cbPalette1 <- brewer.pal(4, "Dark2")
cbPalette2 <- brewer.pal(8, "Paired")[c(2, 3, 6, 7)]  # Lighter
cbPalette2 <- brewer.pal(8, "Paired")[c(2, 4, 6, 8)]  # Darker

# Construct predator dbRDA plot
rdaAll_pred <- ggplot(scores_sites) +
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
       x     = paste("CA1 (", round(100*modBest_pred$CCA$eig[1]/modBest_pred$tot.chi,2), "%)", sep = ""),
       y     = paste("CA1 (", round(100*modBest_pred$CCA$eig[2]/modBest_pred$tot.chi,2), "%)", sep = ""))+
  scale_color_manual(values = c("blue","darkgreen"),
                     guide  = guide_none()) +
  theme_bw() +
  geom_point(aes(x = mult_scores_pred*dbRDA1,
                 y = mult_scores_pred*dbRDA2),
             size  = 1,
             alpha = 0.1) +
  geom_segment(data = vecs_all,
               aes(x     = 0, 
                   xend  = mult_vars_pred*dbRDA1, 
                   y     = 0, 
                   yend  = mult_vars_pred*dbRDA2, 
                   color = sub),
               arrow = arrow(length = unit(.1,"inches")),
               size  = 1,
               alpha = 0.6) +
  geom_label_repel(data = vecs_all,
                   aes(x      = mult_vars_pred*dbRDA1,
                       y      = mult_vars_pred*dbRDA2,
                       label  = rownames(vecs_all),
                       colour = sub,
                       size   = sub),
                   force         = 0.1,
                   alpha         = 0.8,
                   fontface      = "bold",
                   label.padding = 0.25,
                   box.padding   = 0.1,
                   fill          ="white",
                   label.r       = 0.25,
                   label.size    = 0.25,
                   segment.alpha =0,
                   nudge_x       = ifelse(vecs_all$dbRDA1>0,0.1,-0.1),
                   nudge_y       = ifelse(vecs_all$dbRDA2>0,0.2,-0.2)) +
  scale_size_manual(values = c(2.5,3.5),
                    guide  = guide_none()) +
  theme(panel.grid = element_blank())
ggsave("RDA_plot_all.tiff", rdaAll_pred, width = 8, height = 5, dpi = 300)

# Save plot and model data
save(rdaAll_pred, vecs_all, modBest_pred, file = "PredRDA_Plot.RData")
