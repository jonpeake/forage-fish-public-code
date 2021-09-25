# Load required packages
library(vegan)
library(ggplot2)
library(readxl)
library(readr)
library(reshape2)
library(colortools)
library(RColorBrewer)
library(spectral)
library(rsoi)

# Load and parse data
dat                  <- read_csv("Forage_Permanova.csv")
RichEvenSamp         <- read_excel("RichEvenSamp.xlsx")
RichEvenSamp$nseines <- Forage_Permanova$nseines
RichEvenSamp$ntrawls <- Forage_Permanova$ntrawls
tax                  <- dat[, 9:58]
dis                  <- vegdist(tax^0.25)


samp_info <- dat[, 7:8]

fact <- as.data.frame(dat[, 1:6])

# Clear original variable
rm(dat)

# Detrend richness and evenness data for sample effort
detrendMod_Rich        <- lm(S ~ nseines + ntrawls, data = RichEvenSamp)
S_detrend              <- detrendMod_Rich$residuals + detrendMod_Rich$coefficients[1]
RichEvenSamp$S_detrend <- S_detrend

detrendMod_Even                                          <- lm(`J'` ~ nseines, data = RichEvenSamp)
J_detrend                                                <- detrendMod_Even$residuals + detrendMod_Even$coefficients[1]
RichEvenSamp[is.na(RichEvenSamp$`J'`) == 0, "J_detrend"] <- J_detrend

# Combine detrended data with factors, remove excess years
RichEvenSamp$t   <- RichEvenSamp$Year + (RichEvenSamp$Season - 1)/4
RichEvenSamp$Bay <- factor(RichEvenSamp$Bay, levels = c("AP", "CK", "TB", "CH"))
RichEvenSamp     <- RichEvenSamp[!((RichEvenSamp$Bay == "TB" | RichEvenSamp$Bay == "CH") & RichEvenSamp$t < 2005), ]

# Detrend distance matrix for sample effort, isolate residuals and eigenvalues
detrend_mod    <- dbrda(dis ~ nseines + ntrawls, dat = samp_info, sqrt.dist = T)
detrend_scores <- detrend_mod$CA$u %*% diag(sqrt(detrend_mod$CA$eig))
detrend_pctVar <- detrend_mod$CA$eig/sum(detrend_mod$CA$eig)
detrend_cumsum <- cumsum(detrend_pctVar)
detrend_cumsum

# Combine detrended scores with factors, remove excess years
PCOA_Scores     <- data.frame(fact, detrend_scores)
PCOA_Scores$Bay <- as.factor(PCOA_Scores$Bay)
PCOA_Scores$t   <- PCOA_Scores$Year + (PCOA_Scores$Season - 1)/4
PCOA_Scores     <- PCOA_Scores[!((PCOA_Scores$Bay == "TB" | PCOA_Scores$Bay == "CH") & PCOA_Scores$t < 2005), ]

# Calculate summary stats for time series
pco_timeseries_means           <- recast(PCOA_Scores, Bay + t ~ variable, id.var = c("Bay", "t"), measure.var = 7:106, fun.aggregate = mean)
pco_timeseries_se              <- recast(PCOA_Scores, Bay + t ~ variable, id.var = c("Bay", "t"), measure.var = 7:106, fun.aggregate = function(x) {
    sd(x)/sqrt(length(x))
})

rich_even_bay_timeseries_means <- recast(RichEvenSamp, Bay + Year + Season + t ~ variable, id.var = c(6, 9, 10, 15), measure.var = 13:14, fun.aggregate = function(x) {
    mean(x, na.rm = TRUE)
})
rich_even_bay_timeseries_se    <- recast(RichEvenSamp, Bay + Year + Season + t ~ variable, id.var = c(6, 9, 10, 15), measure.var = 13:14, fun.aggregate = function(x) {
    sd(x, na.rm = TRUE)/sqrt(length(x))
})


# Set up lists to drop all analyses and plots
pco_bay <- list()
for (i in levels(pco_timeseries_means$Bay)) {
    pco_bay[[i]][["means"]] <- pco_timeseries_means[pco_timeseries_means$Bay == i, 2:102]
    pco_bay[[i]][["se"]]    <- pco_timeseries_se[pco_timeseries_se$Bay == i, 2:102]
}

rich_bay <- list()

even_bay <- list()


#### Spectral analysis ####

# Isolate pco scores, richness, and evenness by estuary
pco_bay_all <- list()
for (i in levels(PCOA_Scores$Bay)) {
    pco_bay_all[[i]] <- PCOA_Scores[PCOA_Scores$Bay == i, c(which(colnames(PCOA_Scores) == "t"), 7:106)]
}

for (i in levels(RichEvenSamp$Bay)) {
    rich_bay[[i]] <- RichEvenSamp[RichEvenSamp$Bay == i, c("t", "S_detrend")]
}

for (i in levels(RichEvenSamp$Bay)) {
    even_bay[[i]] <- RichEvenSamp[RichEvenSamp$Bay == i, c("t", "J_detrend")]
    even_bay[[i]] <- even_bay[[i]][is.na(even_bay[[i]]$J_detrend) == FALSE, ]
}

pco_bay_spec_all <- list()
pco_sig_per_all  <- list()

# Isolate and remove linear trend
pco_bay_detrend_all <- lapply(pco_bay_all, FUN = function(x) {
    t      <- x$t - min(x$t)
    resids <- data.frame(t = t)
    mod    <- list()
    for (i in 1:20) {
        pco             <- x[, i + 1]
        mod[[i]]        <- lm(pco ~ t)
        resids[, i + 1] <- resid(mod[[i]])
    }
    return(list(resids = resids, mod = mod))
})

rich_bay_detrend <- lapply(rich_bay, FUN = function(x) {
    mod <- lm(x$S_detrend ~ x$t)
    data.frame(t = x$t - min(x$t), rich = resid(mod))
})

even_bay_detrend <- lapply(even_bay, FUN = function(x) {
    mod <- lm(x$J_detrend ~ x$t)
    data.frame(t = x$t - min(x$t), even = resid(mod))
})

# Conduct lomb-scargle spectral analysis on first 20 pco axes by estuary, isolate significant periodicity
for (i in 1:20) {
    pco_bay_spec_all[[paste("PCO", i, sep = "")]] <- lapply(pco_bay_detrend_all, FUN = function(x) spec.lomb(x$resids$t, x$resids[, i + 1], f = seq(0.01, 2, 0.01), mode = "generalized"))
    pco_sig_per_all[[paste("PCO", i, sep = "")]]  <- lapply(pco_bay_spec_all[[paste("PCO", i, sep = "")]], FUN = function(x) {
        fAp     <- data.frame(f = x$f, A = x$A, p = x$p)
        fAp_sig <- fAp[fAp$p < 0.05, ]
        f_sig   <- fAp_sig$f
        per_sig <- 1/f_sig
        return(per_sig)
    })
}

# Conduct spectral analysis for richness and evennness by estuary, isolate significant periodicity
rich_bay_spec <- lapply(rich_bay_detrend, FUN = function(x) spec.lomb(x$t, x$rich, f = seq(0.01, 2, 0.01), mode = "generalized"))
rich_sig_per  <- lapply(rich_bay_spec, FUN = function(x) {
    f_PSD_p     <- data.frame(f = x$f, A = x$PSD, p = x$p)
    f_PSD_p_sig <- f_PSD_p[f_PSD_p$p < 0.05, ]
    f_sig       <- f_PSD_p_sig$f
    per_sig     <- 1/f_sig
    return(per_sig)
})

even_bay_spec <- lapply(even_bay_detrend, FUN = function(x) spec.lomb(x$t, x$even, f = seq(0.01, 2, 0.01), mode = "generalized"))
even_sig_per  <- lapply(even_bay_spec, FUN = function(x) {
    f_PSD_p     <- data.frame(f = x$f, A = x$PSD, p = x$p)
    f_PSD_p_sig <- f_PSD_p[f_PSD_p$p < 0.05, ]
    f_sig       <- f_PSD_p_sig$f
    per_sig     <- 1/f_sig
    return(per_sig)
})

#### Plots ####

# Prepare richness and evenness time series for plotting
rich_time     <- data.frame(t = rich_even_bay_timeseries_means$t, mean = rich_even_bay_timeseries_means$S_detrend, se = rich_even_bay_timeseries_se$S_detrend, bay = rich_even_bay_timeseries_means$Bay)
rich_time$bay <- factor(rich_time$bay, levels = c("AP", "CK", "TB", "CH"))
even_time     <- data.frame(t = rich_even_bay_timeseries_means$t, mean = rich_even_bay_timeseries_means$J_detrend, se = rich_even_bay_timeseries_se$J_detrend, bay = rich_even_bay_timeseries_means$Bay)
even_time     <- even_time[-148, ]
even_time$bay <- factor(even_time$bay, levels = c("AP", "CK", "TB", "CH"))

# Set color scheme
cbPalette  <- brewer.pal(4, "Dark2")

# Plot richness and evennness time series
rich_timePlot <- ggplot(rich_time) +
  geom_ribbon(aes(x    = t,
                  ymin = mean - se,
                  ymax = mean + se,
                  fill = bay),
    alpha       = 0.3,
    show.legend = FALSE)+
  scale_fill_manual(values = cbPalette) +
  geom_line(aes(x      = t,
                y      = mean,
                colour = bay),
    size = 0.2) +
  scale_x_continuous(limits = c(1998,2018),
                     breaks = c(1998:2018),
                     expand = c(0.005,0.005)) +
  theme_bw() +
  theme(axis.text.x       = element_blank(),
        legend.key.width  = unit(5,"lines"),
        axis.title.y      = element_text(vjust = 0.5, hjust = 0.5),
        axis.text.y       = element_text(size = rel(0.8)),
        legend.key.height = unit(4,"lines")) +
  labs(title  = NULL,
       x      = NULL,
       y      = "\nRichness",
       colour = NULL) +
  scale_color_manual(values = cbPalette,
                     labels = c("AB","CK","TB","CH")) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Rich_time.tiff",rich_timePlot,width = 4, height = .8, dpi = 600)

even_timePlot <- ggplot(even_time) +
  geom_ribbon(aes(x    = t, 
                  ymin = mean - se,
                  ymax = mean + se,
                  fill = bay),
    alpha       = 0.3,
    show.legend = FALSE) +
  scale_fill_manual(values = cbPalette) +
  geom_line(aes(x      = t,
                y      = mean,
                colour = bay),
    size = .2,
    show.legend = FALSE) +
  scale_x_continuous(limits = c(1998,2018),
                     breaks = c(1998:2018),
                     expand = c(0.005,0.005)) +
  theme_bw() +
  theme(axis.text.x       = element_blank(),
        axis.title.y      = element_text(vjust = 0.5,hjust = 0.5),
        axis.text.y       = element_text(size = rel(0.8)),
        legend.key.width  = unit(5,"lines"),
        legend.key.height = unit(4,"lines"))+
  labs(title = NULL,
       x = NULL,
       y = "\nEvenness",
       colour = NULL) +
  scale_color_manual(values = cbPalette,
                     labels = c("AB","CK","TB","CH")) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Even_time.tiff",even_timePlot,width = 4, height = .8, dpi = 600)

# Parse richness and evenness spectrograms
rich                                <- rich_bay_spec
rich_per_PSD                        <- do.call(rbind, lapply(seq_along(rich), FUN = function(j) {
    data.frame(bay = names(rich)[[j]], per = 1/rich[[j]]$f, psd = rich[[j]]$PSD, p = rich[[j]]$p)
}))
rich_per_sig                        <- rich_per_PSD$per
rich_psd_sig                        <- rich_per_PSD$psd
rich_psd_sig[rich_per_PSD$p > 0.05] <- -200
rich_per_PSD_sig                    <- data.frame(bay = rich_per_PSD$bay, rich_per_sig, rich_psd_sig)
rich_per_PSD$bay                    <- factor(rich_per_PSD$bay, levels = c("AP", "CK", "TB", "CH"))
rich_per_PSD_sig$bay                <- factor(rich_per_PSD_sig$bay, levels = c("AP", "CK", "TB", "CH"))

even                                <- even_bay_spec
even_per_PSD                        <- do.call(rbind, lapply(seq_along(even), FUN = function(j) {
    data.frame(bay = names(even)[[j]], per = 1/even[[j]]$f, psd = even[[j]]$PSD, p = even[[j]]$p)
}))
even_per_sig                        <- even_per_PSD$per
even_psd_sig                        <- even_per_PSD$psd
even_psd_sig[even_per_PSD$p > 0.05] <- -200
even_per_PSD_sig                    <- data.frame(bay = even_per_PSD$bay, even_per_sig, even_psd_sig)
even_per_PSD$bay                    <- factor(even_per_PSD$bay, levels = c("AP", "CK", "TB", "CH"))
even_per_PSD_sig$bay                <- factor(even_per_PSD_sig$bay, levels = c("AP", "CK", "TB", "CH"))

# Plot richness and evenness spectrograms
rich_spec_plot <- ggplot(rich_per_PSD) +
  geom_ribbon(dat = rich_per_PSD_sig,
              aes(x    = rich_per_sig,
                  ymin = 0,
                  ymax = rich_psd_sig,
                  fill = bay),
              alpha       = 0.5,
              show.legend = FALSE) +
  geom_line(aes(x      = per,
                y      = psd,
                colour = bay),
            size        = 0.2,
            show.legend = FALSE) +
  scale_fill_manual(values = cbPalette) +
  scale_x_log10(limits = c(0.5,10),
                breaks = c(0.5,1,2,3,4,5,6,7,8,9,10),
                expand = c(0.001,0.001)) +
  scale_y_continuous(limits = c(0,NA)) +
  theme_bw() +
  theme(axis.text.x       = element_blank(),
        legend.key.width  = unit(5,"lines"),
        axis.text.y       = element_text(size = rel(0.8)),
        axis.title.y      = element_text(vjust = 0.5, hjust = 0.5),
        legend.key.height = unit(4,"lines")) +
  labs(title  = NULL,
       x      = NULL, 
       y      = "\nPSD",
       colour = NULL) +
  scale_color_manual(values = cbPalette,
                     labels = c("AB","CK","TB","CH")) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Rich_spec.tiff",rich_spec_plot,width = 4, height = .8, dpi = 600)

even_spec_plot <- ggplot(even_per_PSD) +
  geom_ribbon(dat = even_per_PSD_sig,
              aes(x    = even_per_sig,
                  ymin = 0,
                  ymax = even_psd_sig,
                  fill = bay),
              alpha       = 0.5,
              show.legend = FALSE) +
  geom_line(aes(x      = per,
                y      = psd,
                colour = bay),
            size        = 0.2,
            show.legend = FALSE) +
  scale_fill_manual(values = cbPalette) +
  scale_x_log10(limits = c(0.5,10),
                breaks = c(0.5,1,2,3,4,5,6,7,8,9,10),
                expand = c(0.001,0.001)) +
  scale_y_continuous(limits = c(0,NA)) +
  theme_bw() +
  theme(axis.text.x       = element_blank(),
        legend.key.width  = unit(5,"lines"),
        axis.text.y       = element_text(size = rel(0.8)),
        axis.title.y      = element_text(vjust = 0.5, hjust = 0.5),
        legend.key.height = unit(4,"lines")) +
  labs(title  = NULL,
       x      = NULL,
       y      = "\nPSD",
       colour = NULL) +
  scale_color_manual(values = cbPalette,
                     labels = c("AB","CK","TB","CH")) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Even_spec.tiff",even_spec_plot,width = 4, height = 0.8, dpi = 600)


pco_time   <- list()
pco_plots  <- list()
spec_plots <- list()

for (i in 1:20) {
    # Parse pco
    pco_time[[i]] <- data.frame(t = pco_timeseries_means$t, mean = pco_timeseries_means[, i + 2], se = pco_timeseries_se[, i + 2], bay = pco_timeseries_means$Bay)
    
    # Reorder factor
    pco_time[[i]]$bay <- factor(pco_time[[i]]$bay, levels = c("AP", "CK", "TB", "CH"))
    
    # Plot pco time series
    pco_plots[[i]] <- ggplot(pco_time[[i]]) +
        geom_ribbon(aes(x    = t,
                        ymin = mean-se,
                        ymax = mean+se,
                        fill = bay),
            alpha       = 0.3,
            show.legend = FALSE) +
        geom_line(aes(x      = t,
                      y      = mean,
                      colour = bay),
            size        = .2,
            show.legend = FALSE)+
        scale_fill_manual(values = cbPalette) +
        scale_x_continuous(limits = c(1998,2018),
                           breaks = c(1998:2018),
                           expand = c(0.005,0.005)) +
        theme_bw() +
        theme(axis.text.x       = element_blank(),
              axis.text.y       = element_text(size = rel(0.8)),
              legend.key.width  = unit(5,"lines"),
              axis.title.y      = element_text(vjust = 0.5, hjust = 0.5),
              legend.key.height = unit(4,"lines")) +
        labs(title  = NULL,
             x      = NULL,
             y      = paste("PCO",i,"\n(", round(100*detrend_pctVar[i],2),"%)",sep=""),
             colour = NULL)+
        scale_color_manual(values = cbPalette,
                           labels = c("AB","CK","TB","CH"))+
        guides(colour = guide_legend(override.aes = list(size = 2)))
    ggsave(paste("PCO", i, "_time.tiff", sep = ""), pco_plots[[i]], width = 4, height = 0.8, dpi = 600)
    
    # Parse spectrograms
    x                         <- pco_bay_spec_all[[i]]
    per_PSD                   <- do.call(rbind, lapply(seq_along(x), FUN = function(j) {
        data.frame(bay = names(x)[[j]], per = 1/x[[j]]$f, psd = x[[j]]$PSD, p = x[[j]]$p)
    }))
    per_sig                   <- per_PSD$per
    psd_sig                   <- per_PSD$psd
    psd_sig[per_PSD$p > 0.05] <- -200
    per_PSD_sig               <- data.frame(bay = per_PSD$bay, per_sig, psd_sig)
    per_PSD$bay               <- factor(per_PSD$bay, levels = c("AP", "CK", "TB", "CH"))
    per_PSD_sig$bay           <- factor(per_PSD_sig$bay, levels = c("AP", "CK", "TB", "CH"))
    
    # Plot spectrograms
    spec_plots[[i]] <- ggplot(per_PSD) +
        geom_ribbon(dat = per_PSD_sig,
                    aes(x     = per_sig,
                        ymin  = 0,
                        ymax  = psd_sig,
                        fill  = bay,
                        color = bay),
            alpha       = 0.5,
            show.legend = FALSE,
            size        = 0) +
        geom_line(aes(x      = per,
                      y      = psd,
                      colour = bay),
            size        = .2,
            show.legend = FALSE) +    
        scale_fill_manual(values = cbPalette) +
        scale_x_log10(limits = c(0.5,10),
                      breaks = c(0.5,1,2,3,4,5,6,7,8,9,10),
                      expand = c(0.001,0.001)) +
        scale_y_continuous(limits = c(0,NA)) +
        theme_bw() +
        theme(axis.text.x       = element_blank(),
              axis.text.y       = element_text(size = rel(0.8)),
              legend.key.width  = unit(5,"lines"),
              axis.title.y      = element_text(vjust=0.5,hjust=0.5),
              legend.key.height = unit(4,"lines")) +
        labs(title  = NULL,
             x      = NULL,
             y      = "\nPSD",
             colour = NULL) +
        scale_color_manual(values = cbPalette,
                           labels = c("AB","CK","TB","CH")) +
        guides(colour = guide_legend(override.aes = list(size=2)))          
    ggsave(paste("PCO", i, "_spec.tiff", sep = ""), spec_plots[[i]], width = 4, height = 0.8, dpi = 600)
}



#### Synchrony ####

# Establish variables for synchrony
ts_sync              <- list(inter = list(), intra = list())
table_inter_RichEven <- data.frame()
table_intra_RichEven <- data.frame()
table_inter_all      <- data.frame()
table_intra_all      <- data.frame()

# Split into estuaries
rich_bay_time <- list()
for (i in levels(rich_time$bay)) {
    rich_bay_time[[i]] <- rich_time[rich_time$bay == i, c("t", "mean")]
}

even_bay_time <- list()
for (i in levels(even_time$bay)) {
    even_bay_time[[i]] <- even_time[even_time$bay == i, c("t", "mean")]
    even_bay_time[[i]] <- even_bay_time[[i]][is.na(even_bay_time[[i]]$mean) == FALSE, ]
}

# Decompose into inter- and intra-annual signals
series_rich <- lapply(rich_bay_time, FUN = function(x) {
    timeser       <- ts(x$mean, frequency = 4, start = min(x$t))
    timeser_comps <- decompose(timeser)
    interTS       <- x$mean - timeser_comps$seasonal
    intraTS       <- ts(timeser_comps$seasonal, frequency = 4, start = min(x$t))
    out           <- list(interTS = interTS, intraTS = intraTS)
    out
})

series_even <- lapply(even_bay_time, FUN = function(x) {
    timeser       <- ts(x$mean, frequency = 4, start = min(x$t))
    timeser_comps <- decompose(timeser)
    interTS       <- x$mean - timeser_comps$seasonal
    intraTS       <- ts(timeser_comps$seasonal, frequency = 4, start = min(x$t))
    out           <- list(interTS = interTS, intraTS = intraTS)
    out
})

# Establish pairwise companion list
compList    <- c("AP_CK", "AP_TB", "AP_CH", "CK_CH", "CK_TB", "TB_CH")
table_inter <- data.frame()
table_intra <- data.frame()

# Conduct pairwise synchrony tests
for (j in c("AP", "CK", "TB", "CH")) {
    for (k in c("AP", "CK", "TB", "CH")) {
        if (paste(j, k, sep = "_") %in% compList) {
          
            # Run ccf for inter-annual richness and isolate results, significance cutoffs
            corr_inter <- ccf(series_rich[[j]][["interTS"]], series_rich[[k]][["interTS"]], na.action = na.pass, main = paste(j, "_", k, " Rich inter", sep = ""), lag.max = 2)
            max_inter  <- corr_inter$acf[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
            lag_inter  <- corr_inter$lag[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
            CI_inter   <- qnorm((1 + 0.95)/2)/sqrt(corr_inter$n.used)
            sig_inter  <- ifelse(abs(max_inter) > CI_inter, 1, 0)
            out_inter  <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_inter), lag = as.numeric(lag_inter), sig = as.numeric(sig_inter))
            
            ts_sync$inter$Rich[[paste(j, k, sep = "_")]][["corr"]]  <- corr_inter
            ts_sync$inter$Rich[[paste(j, k, sep = "_")]][["table"]] <- out_inter
            table_inter                                             <- rbind(table_inter, cbind(stat = "rich", out_inter))
            
            # Run ccf for inter-annual evenness and isolate results, significance cutoffs
            corr_inter <- ccf(series_even[[j]][["interTS"]], series_even[[k]][["interTS"]], na.action = na.pass, main = paste(j, "_", k, " Even inter", sep = ""), lag.max = 2)
            max_inter  <- corr_inter$acf[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
            lag_inter  <- corr_inter$lag[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
            CI_inter   <- qnorm((1 + 0.95)/2)/sqrt(corr_inter$n.used)
            sig_inter  <- ifelse(abs(max_inter) > CI_inter, 1, 0)
            out_inter  <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_inter), lag = as.numeric(lag_inter), sig = as.numeric(sig_inter))
            
            ts_sync$inter$Even[[paste(j, k, sep = "_")]][["corr"]]  <- corr_inter
            ts_sync$inter$Even[[paste(j, k, sep = "_")]][["table"]] <- out_inter
            table_inter                                             <- rbind(table_inter, cbind(stat = "even", out_inter))
            
            # Run ccf for intra-annual richness and isolate results, significance cutoffs
            corr_intra <- ccf(series_rich[[j]][["intraTS"]], series_rich[[k]][["intraTS"]], na.action = na.pass, main = paste(j, "_", k, " Even intra", sep = ""), lag.max = 1)
            max_intra  <- corr_intra$acf[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
            lag_intra  <- corr_intra$lag[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
            CI_intra   <- qnorm((1 + 0.95)/2)/sqrt(corr_intra$n.used)
            sig_intra  <- ifelse(abs(max_intra) > CI_intra, 1, 0)
            out_intra  <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_intra), lag = as.numeric(lag_intra), sig = as.numeric(sig_intra))
            
            ts_sync$intra$Rich[[paste(j, k, sep = "_")]][["corr"]]  <- corr_intra
            ts_sync$intra$Rich[[paste(j, k, sep = "_")]][["table"]] <- out_intra
            table_intra                                             <- rbind(table_intra, cbind(stat = "rich", out_intra))
            
            # Run ccf for intra-annual evenness and isolate results, significance cutoffs
            corr_intra <- ccf(series_even[[j]][["intraTS"]], series_even[[k]][["intraTS"]], na.action = na.pass, main = paste(j, "_", k, " Even intra", sep = ""), lag.max = 0)
            max_intra  <- corr_intra$acf[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
            lag_intra  <- corr_intra$lag[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
            CI_intra   <- qnorm((1 + 0.95)/2)/sqrt(corr_intra$n.used)
            sig_intra  <- ifelse(abs(max_intra) > CI_intra, 1, 0)
            out_intra  <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_intra), lag = as.numeric(lag_intra), sig = as.numeric(sig_intra))
            
            ts_sync$intra$Even[[paste(j, k, sep = "_")]][["corr"]]  <- corr_intra
            ts_sync$intra$Even[[paste(j, k, sep = "_")]][["table"]] <- out_intra
            table_intra                                             <- rbind(table_intra, cbind(stat = "even", out_intra))
        }
    }
}

# Run through first 20 PCOs
for (i in 1:20) {
    
    # Decompose time series into inter- and intra-annual signals
    series_pco <- lapply(pco_bay, FUN = function(x) {
        timeser       <- ts(x$means[, i + 1], frequency = 4, start = min(x$means$t))
        timeser_comps <- decompose(timeser)
        interTS       <- x$means[, i + 1] - timeser_comps$seasonal
        intraTS       <- ts(timeser_comps$seasonal, frequency = 4, start = min(x$means$t))
        out           <- list(interTS = interTS, intraTS = intraTS)
        out
    })
    
    # Establish pairwise comparison list
    compList    <- c("AP_CK", "AP_TB", "AP_CH", "CK_CH", "CK_TB", "TB_CH")
    table_inter <- data.frame()
    table_intra <- data.frame()
    
    # Conduct pairwise synchrony tests
    for (j in c("AP", "CK", "TB", "CH")) {
        for (k in c("AP", "CK", "TB", "CH")) {
            if (paste(j, k, sep = "_") %in% compList) {
                
                # Run ccf for inter-annual and isolate results, significance cutoffs
                corr_inter <- ccf(series_pco[[j]][["interTS"]], series_pco[[k]][["interTS"]], na.action = na.pass, main = paste(j, "_", k, " PCO", i, " inter", sep = ""), lag.max = 2)
                max_inter  <- corr_inter$acf[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
                lag_inter  <- corr_inter$lag[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
                CI_inter   <- qnorm((1 + 0.95)/2)/sqrt(corr_inter$n.used)
                sig_inter  <- ifelse(abs(max_inter) > CI_inter, 1, 0)
                out_inter  <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_inter), lag = as.numeric(lag_inter), sig = as.numeric(sig_inter))
                
                # Fill list items with correlation and significance
                ts_sync$inter[[paste("PCO", i, sep = "")]][[paste(j, k, sep = "_")]][["corr"]]  <- corr_inter
                ts_sync$inter[[paste("PCO", i, sep = "")]][[paste(j, k, sep = "_")]][["table"]] <- out_inter
                table_inter                                                                     <- rbind(table_inter, cbind(stat = paste("PCO",i,sep = ""), out_inter))
                
                
                # Run ccf for intra-annual and isolate results, significance cutoffs
                corr_intra <- ccf(series_pco[[j]][["intraTS"]], series_pco[[k]][["intraTS"]], na.action = na.pass, main = paste(j, "_", k, " PCO", i, " intra", sep = ""), lag.max = 0)
                max_intra  <- corr_intra$acf[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
                lag_intra  <- corr_intra$lag[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
                CI_intra   <- qnorm((1 + 0.95)/2)/sqrt(corr_intra$n.used)
                sig_intra  <- ifelse(abs(max_intra) > CI_intra, 1, 0)
                out_intra  <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_intra), lag = as.numeric(lag_intra), sig = as.numeric(sig_intra))
                
                # Fill list items with correlation and significance
                ts_sync$intra[[paste("PCO", i, sep = "")]][[paste(j, k, sep = "_")]][["corr"]]  <- corr_intra
                ts_sync$intra[[paste("PCO", i, sep = "")]][[paste(j, k, sep = "_")]][["table"]] <- out_intra
                table_intra                                                                     <- rbind(table_intra, cbind(stat = paste("PCO",i,sep = ""), out_intra))
            }
        }
    }
    
    # Fill list items with results tables
    ts_sync$inter[[paste("PCO", i, sep = "")]][["table"]] <- table_inter
    ts_sync$intra[[paste("PCO", i, sep = "")]][["table"]] <- table_intra
    table_inter_all <- rbind(table_inter_all, cbind(type = "inter", table_inter))
    table_intra_all <- rbind(table_intra_all, cbind(type = "intra", table_intra))
}

table_inter_RichEven <- rbind(table_inter_RichEven, cbind(type = "inter", table_inter))
table_intra_RichEven <- rbind(table_intra_RichEven, cbind(type = "intra", table_intra))

table_inter <- rbind(table_inter_RichEven,table_inter_all)
table_intra <- rbind(table_intra_RichEven,table_intra_all)

# Fill list items with results tables
ts_sync$inter[["table"]] <- table_inter
ts_sync$intra[["table"]] <- table_intra
ts_sync[["table"]]       <- rbind(table_inter, table_intra)
ts_sync$table$type       <- as.factor(ts_sync$table$type)
ts_sync$table$stat       <- as.factor(ts_sync$table$stat)
ts_sync$table$bay1       <- factor(ts_sync$table$bay1, levels = c("AP", "CK", "TB", "CH"))
ts_sync$table$bay2       <- factor(ts_sync$table$bay2, levels = c("AP", "CK", "TB", "CH"))
ts_sync$table$maxcorr    <- as.numeric(ts_sync$table$maxcorr)
ts_sync$table$lag        <- as.numeric(ts_sync$table$lag)
ts_sync$table$sig        <- as.numeric(ts_sync$table$sig)
ts_sync$table
cast_table               <- recast(ts_sync$table, formula = type + stat + bay2 ~ variable + bay1, measure.var = c("maxcorr", "lag", "sig"), fun.aggregate = mean)

# Export plot data
save(pco_plots, spec_plots, rich_timePlot, rich_spec_plot, even_timePlot, even_spec_plot, file = "TimeSeries_Plots.RData")


### MEI Correlations ###

# Downlaod MEI data using rsoi package
mei <- download_mei()
nao <- download_nao()

# Clean dataset and aggregate data seasonally
yearstart  <- min(mei$Year)
mei        <- mei[order(mei$Date), "MEI"]
mei        <- as.numeric(unlist(mei))
series_mei <- ts(na.omit(mei), frequency = 12, start = yearstart)
series_mei <- aggregate(series_mei, nfrequency = 4, FUN = mean)

yearstart  <- min(nao$Year)
nao        <- nao[, "NAO"]
nao        <- as.numeric(unlist(nao))
series_nao <- ts(na.omit(nao), frequency = 12, start = yearstart)
series_nao <- aggregate(series_nao, nfrequency = 4, FUN = mean)


# Run through each of first 20 PCO axes
for (i in 1:20) {
    
    # Decompose PCOs and isolate inter-annual signal
    series_pco_mei    <- lapply(pco_bay, FUN = function(x) {
        timeser       <- ts(x$means[, i + 1], frequency = 4, start = min(x$means$t))
        timeser_comps <- decompose(timeser)
        interTS       <- x$means[, i + 1] - timeser_comps$seasonal
        out           <- interTS
        out
    })
    
    table_mei <- data.frame()
    table_nao <- data.frame()
    
    # Calculate cross correlations for each estuary
    for (j in c("AP", "CK", "TB", "CH")) {
        pco_mei  <- ts.intersect(series_mei, series_pco_mei[[j]])
        corr_mei <- ccf(series_mei, series_pco_mei[[j]], lag.max = 1)
        max_mei  <- corr_mei$acf[which(abs(corr_mei$acf) == max(abs(corr_mei$acf)))]
        lag_mei  <- corr_mei$lag[which(abs(corr_mei$acf) == max(abs(corr_mei$acf)))]
        CI_mei   <- qnorm((1 + 0.95)/2)/sqrt(corr_mei$n.used)
        sig_mei  <- ifelse(abs(max_mei) > CI_mei, 1, 0)
        out_mei  <- cbind(bay = j, maxcorr = as.numeric(max_mei), lag = as.numeric(lag_mei), sig = as.numeric(sig_mei))
        
        # Fill list items
        ts_sync$mei[[paste("PCO", i, sep = "")]][[j]][["corr"]]  <- corr_mei
        ts_sync$mei[[paste("PCO", i, sep = "")]][[j]][["table"]] <- out_mei
        
        table_mei <- rbind(table_mei, cbind(PCO = i, out_mei))
        
        pco_nao  <- ts.intersect(series_nao, series_pco_mei[[j]])
        corr_nao <- ccf(series_nao, series_pco_mei[[j]], lag.max = 1)
        max_nao  <- corr_nao$acf[which(abs(corr_nao$acf) == max(abs(corr_nao$acf)))]
        lag_nao  <- corr_nao$lag[which(abs(corr_nao$acf) == max(abs(corr_nao$acf)))]
        CI_nao   <- qnorm((1 + 0.95)/2)/sqrt(corr_nao$n.used)
        sig_nao  <- ifelse(abs(max_nao) > CI_nao, 1, 0)
        out_nao  <- cbind(bay = j, maxcorr = as.numeric(max_nao), lag = as.numeric(lag_nao), sig = as.numeric(sig_nao))
        
        # Fill list items
        ts_sync$nao[[paste("PCO", i, sep = "")]][[j]][["corr"]]  <- corr_nao
        ts_sync$nao[[paste("PCO", i, sep = "")]][[j]][["table"]] <- out_nao
        
        table_nao <- rbind(table_nao, cbind(PCO = i, out_nao))
        
    }
    
    # Fill list items
    ts_sync$mei[[paste("PCO", i, sep = "")]][["table"]] <- table_mei
    ts_sync$nao[[paste("PCO", i, sep = "")]][["table"]] <- table_nao
}
