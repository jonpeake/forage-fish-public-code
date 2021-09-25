#### Load Libraries ####
library(ggplot2)
library(ggrepel)
library(colorspace)
library(colortools)
library(RColorBrewer)
library(egg)
library(patchwork)
library(tidyverse)

#### Load Plot Data ####
files <- list.files(pattern = "*.RData")
lapply(files, load, .GlobalEnv)

#### Figure 2: Shadeplot ####

Fig2 <- forage_shade


Fig2

ggsave("Fig2.pdf", Fig2, width = 4.5, height = 9)

#### Figure S1: Abundance and Frequency of Occurrence ####

FigS1a <- rankabun_plots$AB + 
         labs(x     = NULL, 
              title = "Apalachicola Bay")

FigS1b <- rankfreq_plots$AB +
         labs(x = NULL)

FigS1c <- rankabun_plots$CK + 
         labs(x     = NULL, 
              title = "Cedar Key")

FigS1d <- rankfreq_plots$CK + 
         labs(x = NULL)

FigS1e <- rankabun_plots$TB + 
         labs(x     = NULL, 
              title = "Tampa Bay")

FigS1f <- rankfreq_plots$TB + 
         labs(x = NULL)

FigS1g <- rankabun_plots$CH + 
         labs(title = "Charlotte Harbor")

FigS1h <- rankfreq_plots$CH

FigS1  <- ((FigS1a + FigS1b) +
           plot_layout(ncol = 2))/
         ((FigS1c + FigS1d) +
           plot_layout(ncol = 2))/
         (FigS1e + FigS1f +
           plot_layout(ncol = 2))/
         (FigS1g + FigS1h +
           plot_layout(ncol = 2))/
         (guide_area()) +
         plot_layout(nrow = 5) +
         plot_layout(guides  = 'collect',
                     heights = c(1,1,1,1,.01)) +
         plot_annotation(tag_levels = "a")
  

FigS1  <- FigS1 & theme(plot.title.position = "plot",
                      plot.title          = element_text(vjust = 4, face = "bold", size = 16),
                      plot.tag.position   = c(0.1, .95),
                      legend.position     = "bottom",
                      legend.key.height   = unit(0,"inches"),
                      legend.title        = element_blank(),
                      axis.title.y        = element_text(hjust = 0.5, vjust = 0),
                      plot.margin         = unit(c(.1,.1,.1,.1), "inches"))
FigS1

ggsave("FigS1.pdf", FigS1, height = 8, width = 15)

#### Figure S2: Richness and Evenness ####

FigS2a <- RichBay_plot + 
         theme(axis.text.x = element_blank()) + 
         labs(title = "a")

FigS2b <- RichSeason_plot + 
         theme(axis.text  = element_blank(), 
               axis.title = element_blank()) + 
         labs(title = "b")

FigS2c <- EvenBay_plot + 
         labs(title = "c")


FigS2d <- EvenSeason_plot + theme(axis.text.y = element_blank(), axis.title = element_blank()) + labs(title = "d")

FigS2  <- (FigS2a | FigS2b)/
         (FigS2c | FigS2d)

FigS2

ggsave("FigS2.pdf", FigS2, height = 6, width = 6)

#### Figure 3: CAP Plots by Bay and Season ####

# Initialize for loop for Season CAP plots by bay
CAP_plots  <- list(Bay = list(), Season = list())
spp_vecs   <- list(Bay = list(), Season = list())
cbPalette1 <- brewer.pal(8, "Paired")[c(2, 3, 6, 7)]
cbPalette2 <- brewer.pal(4, "Dark2")

cap_pctvars <- function(cap_out){
  # A function that takes the output of the CAPdiscrim function and 
  # produces the percent of among-group variability explained by each
  # canonical axis.
  #
  # Inputs: 
  # cap_out = 'CAPDiscrim' function list output
  # 
  # Output: 
  # pct_var = A numerical array with dimensions 1 x p, where p is the
  #           number of canonical axes, corresponding to the percent
  #           of among-group variability explained by each axis.
  #         
  
  # Extract eigenvalues from manova sublist
  eig     <- cap_out[["manova"]][["Eigenvalues"]]
  
  # Convert eigenvalues to analog of canonical variability
  vars    <- eig*cap_out$tot*(cap_out$varm/100)/(eig+1)
  
  # Calculate percent of variability explained by each axis
  pct_var <- t(100*vars/sum(vars))
  
  # Trim rows not corresponding to canonical axes
  pct_var <- pct_var[1:ncol(cap_out$x),]
  return(pct_var)
}

for (i in c("AP", "CK", "TB", "CH")) {
    # Parse variation
    var_expl <- round(100 * sum(CAP$Season$mods[[i]][["manova"]][["SS"]][["y[, group]"]])/CAP$Season$mods[[i]][["tot"]], 2)
    var_CA   <- round(cap_pctvars(CAP$Season$mods[[i]]), 2)
    
    # Parse data
    site_scores          <- data.frame(CAP$Season$mods[[i]]$x)
    site_scores$Season   <- factor(CAP$Season$mods[[i]]$group, levels = c("1", "2", "3", "4"))
    spp_vecs$Season[[i]] <- data.frame(CAP$Season$mods[[i]]$cproj, labs = rownames(CAP$Season$mods[[i]]$cproj))
    mult                 <- 12
    
    # Plot data
    CAP_plots$Season[[i]] <- ggplot(site_scores) + 
        geom_vline(xintercept = 0,
                   colour     = "grey70",
                   size       = .25) +
        geom_hline(yintercept = 0,
                   colour     = "grey70",
                   size       = .25) +
        scale_x_continuous(limits       = symmetric_range((3+mult)*spp_vecs$Season[[i]]$LD1),
                           breaks       = c(-6,-3,0,3,6),
                           minor_breaks = NULL) +
        scale_y_continuous(limits       = symmetric_range,
                           breaks       = c(-6,-3,0,3,6),
                           minor_breaks = NULL) +
        geom_point(aes(x    = LD1, 
                       y    = LD2, 
                       fill = Season), 
                   size   = 1.5, 
                   stroke = 0.001,
                   pch    = 21, 
                   colour = "black") +
        labs(title = NULL,
             x     = paste("CA1 (",var_CA[1],"%)",sep=""),
             y     = paste("CA2 (",var_CA[2],"%)",sep=""),
             fill  = NULL) +
        scale_fill_manual(values = cbPalette1,
                          labels = c("Wi","Sp","Su","Fa")) +
        theme_bw() +
        theme(legend.text       = element_text(size = rel(0.8)),
              legend.position   = c(0.051,0.89),
              legend.background = element_blank(),
              legend.key        = element_blank(),
              panel.grid        = element_blank()) +
        geom_segment(data = spp_vecs$Season[[i]],
                     aes(x    = 0,
                         xend = mult*LD1,
                         y    = 0,
                         yend = mult*LD2),
                     arrow = arrow(length = unit(0.05,"inches")),
                     color = "grey1",
                     size  = 0.5,
                     alpha = 0.6) +
        geom_label_repel(data = spp_vecs$Season[[i]],
                         aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                             y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                             label = labs),
                         segment.alpha = 0,
                         size          = 2.5,
                         color         = "blue",
                         fontface      = "bold",
                         fill          = "white",
                         alpha         = 0.7,
                         box.padding   = .25,
                         lineheight    = 0.4,
                         label.size    = 0.25,
                         nudge_x       = ifelse(spp_vecs$Season[[i]]$LD1>0,0.05,-0.05),
                         nudge_y       = ifelse(spp_vecs$Season[[i]]$LD2>-0.02,0.05,-0.05),
                         force         = 0.3) +
        guides(fill = guide_legend(override.aes = list(size = 2.5),
                                   keyheight    = 0.15,
                                   keywidth     = 0.2))
}

# Initialize for loop for Bay CAP plots by season
for (i in c("1", "2", "3", "4")) {
    # Parse variation
    var_expl <- round(100 * sum(CAP$Bay$mods[[i]][["manova"]][["SS"]][["y[, group]"]])/CAP$Bay$mods[[i]][["tot"]], 2)
    var_CA   <- round(cap_pctvars(CAP$Bay$mods[[i]]), 2)
    
    # Parse data
    site_scores       <- data.frame(CAP$Bay$mods[[i]]$x)
    site_scores$Bay   <- factor(CAP$Bay$mods[[i]]$group, levels = c("AP", "CK", "TB", "CH"))
    spp_vecs$Bay[[i]] <- data.frame(CAP$Bay$mods[[i]]$cproj, labs = rownames(CAP$Bay$mods[[i]]$cproj))
    mult              <- 12
    
    # Plot data
    CAP_plots$Bay[[i]] <- ggplot(site_scores) +
        geom_vline(xintercept = 0,
                   colour     = "grey70",
                   size       = .25) +
        geom_hline(yintercept = 0,
                   colour     = "grey70",
                   size       = .25) +
        scale_x_continuous(limits       = symmetric_range((1+mult)*spp_vecs$Bay[[i]]$LD1),
                           breaks       = c(-6,-3,0,3,6),
                           minor_breaks = NULL) +
        scale_y_continuous(limits       = symmetric_range,
                           breaks       = c(-6,-3,0,3,6),
                           minor_breaks = NULL) +
        geom_point(aes(x    = LD1, 
                       y    = LD2, 
                       fill = Bay),
                   size   = 1.5, 
                   stroke = 0.001,
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
        geom_segment(data = spp_vecs$Bay[[i]],
                     aes(x    = 0,
                         xend = mult*LD1,
                         y    = 0,
                         yend = mult*LD2),
                     arrow = arrow(length = unit(0.05,"inches")),
                     color = "grey1",
                     size  = 0.5,
                     alpha = 0.6)+
        geom_label_repel(data = spp_vecs$Bay[[i]],
                         aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                             y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                             label = labs),
                         segment.alpha = 0,
                         size          = 2.5,
                         color         = "blue",
                         fontface      = "bold",
                         fill          = "white",
                         alpha         = 0.7,
                         box.padding   = .25,
                         lineheight    = 0.4,
                         label.size    = 0.25,
                         nudge_x       = ifelse(spp_vecs$Bay[[i]]$LD1>0,0.05,-0.05),
                         nudge_y       = ifelse(spp_vecs$Bay[[i]]$LD2>-0.02,0.05,-0.05),
                         force         = 0.3) +
        guides(fill = guide_legend(override.aes = list(size = 2.5),
                                   keyheight    = 0.15,
                                   keywidth     = 0.2))    
}

Fig3a <- CAP_plots$Bay$`1` + 
         labs(title = "a) Winter (15.855)")
Fig3a

Fig3b <- CAP_plots$Bay$`2` + 
         labs(title = "b) Spring (6.56%)")
Fig3b

Fig3c <- CAP_plots$Bay$`3` + 
         labs(title = "c) Summer (2.05%)")
Fig3c

Fig3d <- CAP_plots$Bay$`4` + 
         labs(title = "d) Fall (10.63%)")
Fig3d

Fig3e <- CAP_plots$Season$AP + 
         labs(title = "e) Apalachicola Bay (6.23%)")
Fig3e

Fig3f <- CAP_plots$Season$CK + 
         labs(title = "f) Cedar Key (15.28%)")
Fig3f

Fig3g <- CAP_plots$Season$TB + 
         labs(title = "g) Tampa Bay (0.26%)")
Fig3g

Fig3h <- CAP_plots$Season$CH + 
         labs(title = "h) Charlotte Harbor (7.55%)")
Fig3h

Fig3 <- (Fig3a + Fig3b  + 
         Fig3c + Fig3d  + 
         Fig3e + Fig3f  + 
         Fig3g + Fig3h) + 
         plot_layout(ncol = 2)
Fig3

ggsave("Fig3.pdf", Fig3, width = 8.5, height = 11)


#### Figure 4 and S2S1: RDA Plot ####

# Set scaling factors
mult_vars_pred   <- 9
mult_scores_pred <- 1

mult_vars_hab    <- 7
mult_scores_hab  <- 2.5

Fig4    <- rdaAll_hab + 
         labs(title = "Total variability explained: 11.82%")

FigS2S1 <- rdaAll_pred + 
         labs(title = "Total variability explained: 3.64%")

ggsave("Fig4.pdf", Fig4, width = 8.5, height = 5.5)
ggsave("FigS2S1.pdf", FigS2S2, width = 8.5, height = 5.5)

#### Figure S3: Time Series ####
FigS3a <- rich_timePlot + 
         theme(axis.title.y = element_text(hjust = 1, vjust = 1))

FigS3b <- rich_spec_plot

FigS3c <- even_timePlot + 
         theme(axis.title.y = element_text(hjust = 1, vjust = 1))

FigS3d <- even_spec_plot

FigS3e <- pco_plots[[1]] + 
         theme(axis.title.y = element_text(hjust = 0, vjust = 1))

FigS3f <- spec_plots[[1]]

FigS3g <- pco_plots[[2]] + 
         theme(axis.title.y = element_text(hjust = 0, vjust = 1))

FigS3h <- spec_plots[[2]]

FigS3i <- pco_plots[[3]] + 
         theme(axis.title.y = element_text(hjust = 0, vjust = 1))

FigS3j <- spec_plots[[3]]

FigS3k <- pco_plots[[4]] + 
         theme(axis.title.y = element_text(hjust = 0, vjust = 1), 
               axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.7)))

FigS3l <- spec_plots[[4]] + 
         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.7)))

FigS3 <- 
  (guide_area())/
  (FigS3a + FigS3b +
   FigS3c + FigS3d +
   FigS3e + FigS3f +
   FigS3g + FigS3h +
   FigS3i + FigS3j +
   FigS3k + FigS3l +
   plot_layout(ncol = 2)) +
   plot_layout(guides = 'collect',heights = c(.05,1)) +
   plot_annotation(tag_levels = "a")

FigS3 <- FigS3 & 
        theme(plot.tag.position = c(0, 1.2),
              legend.position   = "top",
              legend.key.height = unit(0,"inches"),
              axis.title.y      = element_text(hjust = 0.5,vjust = 0),
              plot.margin       = unit(c(.15,.1,.15,.1),"inches"))

FigS3

ggsave("FigS3.pdf", FigS3, height = 6, width = 8.5)



