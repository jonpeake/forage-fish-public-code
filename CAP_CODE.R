# Load Required Packages
library(BiodiversityR)
library(readr)
library(labdsv)
library(ggplot2)
library(colortools)
library(ggrepel)
library(RColorBrewer)
library(egg)

# Parse Data
dat       <- read_csv("Forage_Permanova.csv")
tax       <- dat[, 9:58]
tax       <- tax[, taxID]
dis       <- vegdist(tax^0.25)

samp_info <- dat[, 7:8]

fact      <- as.data.frame(dat[, 1:6])

# Clear original variable
rm(dat)

# Detrend for sample size
detrend_mod    <- dbrda(dis ~ nseines + ntrawls, dat = samp_info, sqrt.dist = T)
detrend_scores <- detrend_mod$CA$u %*% diag(sqrt(detrend_mod$CA$eig))



# Convert factor columns to factors
fact$Bay    <- as.factor(fact$Bay)
fact$Zone   <- as.factor(fact$Zone)
fact$Year   <- as.factor(fact$Year)
fact$Season <- as.factor(fact$Season)

# Calculate indvals and isolate top taxa

clust_all          <- rep(1, nrow(tax))
indval_all         <- indval(tax, clust_all)
avg_abund          <- apply(tax, 2, mean)
rel_abund          <- avg_abund/sum(avg_abund)
indvals_all        <- data.frame(indval = indval_all$relfrq * rel_abund, 
                                 relfrq = indval_all$relfrq, 
                                 relabn = rel_abund)
indval_top         <- indvals_all[order(indvals_all$indval, decreasing = TRUE), ]
indval_top15       <- indval_top[1:15, ]


load("TaxaID.RData")
colnames(tax) <- as.character(1:dim(tax)[2])

cbPalette1    <- brewer.pal(8, "Paired")[c(2, 3, 6, 7)]
cbPalette2    <- brewer.pal(4, "Dark2")

CAP <- list(Season = list(mods = list(), plots = list()), Bay = list(mods = list(), plots = list()))


# Initialize for loop for Season CAPs by bay
for (i in levels(fact$Bay)) {
    print(i)
    # Subset Data
    dis_temp            <- dist((detrend_scores[fact$Bay == i, ]))
    fact_temp           <- fact[fact$Bay == i, ]
    tax_temp            <- tax[fact$Bay == i, ]
    tax_temp            <- tax_temp[, colSums(tax_temp) > 0]
    
    clust_temp          <- rep(1, nrow(tax_temp))
    indval_temp         <- indval(tax_temp, clust_temp)
    avg_abund_temp      <- apply(tax_temp, 2, mean)
    rel_abund_temp      <- avg_abund_temp/sum(avg_abund_temp)
    indvals_temp        <- data.frame(indval = indval_temp$relfrq * rel_abund_temp, 
                                      relfrq = indval_temp$relfrq, 
                                      relabn = rel_abund_temp)
    indval_temp_top     <- indvals_temp[order(indvals_temp$indval, decreasing = TRUE), ]
    indval_temp_top10   <- indval_temp_top[1:10, ]
    tax_subset          <- tax_temp[, rownames(indval_temp_top10)]
    
    # Conduct analysis
    CAP$Season$mods[[i]] <- CAPdiscrim(dis_temp ~ Season, data = fact_temp, permutations = 0, mmax = 500)
    CAP$Season$mods[[i]] <- add.spec.scores(CAP$Season$mods[[i]], tax_subset)
    
}



# Initialize for loop for Bay CAPs by season
for (i in levels(fact$Season)) {
    print(i)
    # Subset Data
    dis_temp            <- dist((detrend_scores[fact$Season == i, ]))
    fact_temp           <- fact[fact$Season == i, ]
    tax_temp            <- tax[fact$Season == i, ]
    tax_temp            <- tax_temp[, colSums(tax_temp) > 0]
    
    clust_temp          <- rep(1, nrow(tax_temp))
    indval_temp         <- indval(tax_temp, clust_temp)
    avg_abund_temp      <- apply(tax_temp, 2, mean)
    rel_abund_temp      <- avg_abund_temp/sum(avg_abund_temp)
    indvals_temp        <- data.frame(indval = indval_temp$relfrq * rel_abund_temp, 
                                      relfrq = indval_temp$relfrq, 
                                      relabn = rel_abund_temp)
    indval_temp_top     <- indvals_temp[order(indvals_temp$indval, decreasing = TRUE), ]
    indval_temp_top10   <- indval_temp_top[1:10, ]
    tax_subset          <- tax_temp[, rownames(indval_temp_top10)]
    
    # Conduct analysis
    CAP$Bay$mods[[i]] <- CAPdiscrim(dis_temp ~ Bay, data = fact_temp, permutations = 0, mmax = 500)
    CAP$Bay$mods[[i]] <- add.spec.scores(CAP$Bay$mods[[i]], tax_subset)
}

# Initialize for loop for Season CAP plots by bay
for (i in levels(fact$Bay)) {
    # Parse variation
    var_expl <- round(100 * sum(CAP$Season$mods[[i]][["manova"]][["SS"]][["y[, group]"]])/CAP$Season$mods[[i]][["tot"]], 2)
    var_CA <- round(100 * CAP$Season$mods[[i]][["lda.other"]][["svd"]]^2/sum(CAP$Season$mods[[i]][["lda.other"]][["svd"]]^2), 2)
    
    # Parse data
    site_scores <- data.frame(CAP$Season$mods[[i]]$x)
    site_scores$Season <- factor(CAP$Season$mods[[i]]$group, levels = c("1", "2", "3", "4"))
    spp_vecs <- data.frame(CAP$Season$mods[[i]]$cproj)
    mult <- 12
    
    # Plot data
    CAP$Season$plots[[i]] <- ggplot(site_scores) +
        geom_vline(xintercept = 0,
                   colour     = "grey70",
                   size       = .25) +
        geom_hline(yintercept = 0,
                   colour     = "grey70",
                   size       = .25) +
        scale_x_continuous(limits       = symmetric_range((3+mult)*spp_vecs$LD1),
                           breaks       = c(-6,-3,0,3,6),
                           minor_breaks = NULL) +
        scale_y_continuous(limits       = symmetric_range,
                           breaks       = c(-6,-3,0,3,6),
                           minor_breaks = NULL) +
        geom_point(aes(x    = LD1, 
                       y    = LD2, 
                       fill = Season), 
                   size   = 1, 
                   stroke = 0.1,
                   pch    = 21, 
                   colour = "black") +
        labs(title = NULL,
             x     = paste("CA1 (",var_CA[1],"%)",sep=""),
             y     = paste("CA2 (",var_CA[2],"%)",sep=""),
             fill  = NULL) +
        scale_fill_manual(values = cbPalette1,
                          labels = c("Wi","Sp","Su","Fa")) +
        theme_bw() +
        theme(legend.text       = element_text(size=rel(0.8)),
              legend.position   = c(0.051,0.89),
              legend.background = element_blank(),
              legend.key        = element_blank(),
              panel.grid        = element_blank()) +
        geom_segment(data = spp_vecs,
                     aes(x    = 0,
                         xend = mult*LD1,
                         y    = 0,
                         yend = mult*LD2),
                     arrow = arrow(length = unit(0.1,"inches")),
                     color = "grey1",
                     size  = .75,
                     alpha = 0.6) +
        geom_label_repel(data = spp_vecs,
                         aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                             y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                             label = rownames(spp_vecs)),
                         segment.alpha = 0,
                         size          = 3.25,
                         color         = "blue",
                         fontface      = "bold",
                         fill          = "white",
                         alpha         = 0.7,
                         box.padding   = .25,
                         lineheight    = 0.4,
                         label.size    = 0.25,
                         nudge_x       = ifelse(spp_vecs$LD1>0,0.05,-0.05),
                         nudge_y       = ifelse(spp_vecs$LD2>-0.02,0.05,-0.05),
                         force         = 0.3) +
        guides(fill = guide_legend(override.aes = list(size = 2.5),
                                   keyheight    = 0.15,
                                   keywidth     = 0.2))
    ggsave(paste("CAP_Season_", i, ".tiff", sep = ""), CAP$Season$plots[[i]], width = 4, height = 2.5, dpi = 600)
}

# Initialize for loop for Bay CAP plots by season
for (i in levels(fact$Season)) {
    # Parse variation
    var_expl <- round(100 * sum(CAP$Bay$mods[[i]][["manova"]][["SS"]][["y[, group]"]])/CAP$Bay$mods[[i]][["tot"]], 2)
    var_CA <- round(100 * CAP$Bay$mods[[i]][["lda.other"]][["svd"]]^2/sum(CAP$Bay$mods[[i]][["lda.other"]][["svd"]]^2), 2)
    
    # Parse data
    site_scores <- data.frame(CAP$Bay$mods[[i]]$x)
    site_scores$Bay <- factor(CAP$Bay$mods[[i]]$group, levels = c("AP", "CK", "TB", "CH"))
    spp_vecs <- data.frame(CAP$Bay$mods[[i]]$cproj)
    mult <- 12
    
    # Plot data
    CAP$Bay$plots[[i]] <- ggplot(site_scores) + 
        geom_vline(xintercept = 0,
                   colour     = "grey70",
                   size       = .25) +
        geom_hline(yintercept = 0,
                   colour     = "grey70",
                   size       = .25) +
        scale_x_continuous(limits       = symmetric_range((1+mult)*spp_vecs$LD1),
                           breaks       = c(-6,-3,0,3,6),
                           minor_breaks = NULL) +
        scale_y_continuous(limits       = symmetric_range,
                           breaks       = c(-6,-3,0,3,6),
                           minor_breaks = NULL) +
        geom_point(aes(x    = LD1, 
                       y    = LD2, 
                       fill = Bay),
                   size   = 1, 
                   stroke = 0.1,
                   pch    = 21, 
                   colour = "black") +
        labs(title = NULL,
             x     = paste("CA1 (",var_CA[1],"%)",sep=""),
             y     = paste("CA2 (",var_CA[2],"%)",sep=""),
             fill  = NULL) +
        scale_fill_manual(values = cbPalette2,
                          labels = c("AB","CK","TB","CH")) +
        theme_bw() +
        theme(legend.text       = element_text(size=rel(0.8)),
              legend.position   = c(0.055,0.89),
              legend.background = element_blank(),
              legend.key        = element_blank(),
              panel.grid        = element_blank()) +
        geom_segment(data = spp_vecs,
                     aes(x    = 0,
                         xend = mult*LD1,
                         y    = 0,
                         yend = mult*LD2),
                     arrow = arrow(length = unit(0.1,"inches")),
                     color = "grey1",
                     size  = .75,
                     alpha = 0.6)+
        geom_label_repel(data = spp_vecs,
                         aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                             y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                             label = rownames(spp_vecs)),
                         segment.alpha = 0,
                         size          = 3.25,
                         color         = "blue",
                         fontface      = "bold",
                         fill          = "white",
                         alpha         = 0.7,
                         box.padding   = .25,
                         lineheight    = 0.4,
                         label.size    = 0.25,
                         nudge_x       = ifelse(spp_vecs$LD1>0,0.05,-0.05),
                         nudge_y       = ifelse(spp_vecs$LD2>-0.02,0.05,-0.05),
                         force         = 0.3) +
        guides(fill = guide_legend(override.aes = list(size = 2.5),
                                   keyheight    = 0.15,
                                   keywidth     = 0.2))
    ggsave(paste("CAP_Bay_", c("Wi", "Sp", "Su", "Fa")[as.numeric(i)], ".tiff", sep = ""), CAP$Bay$plots[[i]], width = 4, height = 2.5, dpi = 600)
    
}

# Save plot data
save(CAP, file = "CAP_Plots.RData")








