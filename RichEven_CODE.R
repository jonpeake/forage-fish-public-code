#### Load required packages ####
library(vegan)
library(ggplot2)
library(readxl)
library(readr)
library(reshape2)
library(colortools)
library(RColorBrewer)
library(spectral)
library(reshape2)
library(rcompanion)
library(lmPerm)
library(RVAideMemoire)

#### Massage data ####
RichEvenSamp <- read_excel("RichEvenSamp.xlsx")
Forage_Permanova <- read_csv("Forage_Permanova.csv")
RichEvenSamp$nseines <- Forage_Permanova$nseines
RichEvenSamp$ntrawls <- Forage_Permanova$ntrawls

detrendMod_Rich <- lm(S ~ nseines + ntrawls, data = RichEvenSamp)
S_detrend <- detrendMod_Rich$residuals + detrendMod_Rich$coefficients[1]

RichEvenSamp$S_detrend <- S_detrend

detrendMod_Even <- lm(`J'` ~ nseines, data = RichEvenSamp)
J_detrend <- detrendMod_Even$residuals + detrendMod_Even$coefficients[1]
RichEvenSamp[is.na(RichEvenSamp$`J'`) == 0, "J_detrend"] <- J_detrend
RichEvenSamp$t <- RichEvenSamp$Year + (RichEvenSamp$Season - 1)/4
RichEvenSamp$Bay <- factor(RichEvenSamp$Bay, levels = c("AP", "CK", "TB", "CH"))
RichEvenSamp <- RichEvenSamp[!((RichEvenSamp$Bay == "TB" | RichEvenSamp$Bay == "CH") & RichEvenSamp$t < 2005), ]

rich_even_bay_timeseries_means <- recast(RichEvenSamp, Bay + Year + Season + t ~ variable, id.var = c(6, 9, 10, 15), measure.var = 13:14, fun.aggregate = function(x) {
    mean(x, na.rm = TRUE)
})
rich_even_bay_timeseries_se <- recast(RichEvenSamp, Bay + Year + Season + t ~ variable, id.var = c(6, 9, 10, 15), measure.var = 13:14, fun.aggregate = function(x) {
    sd(x, na.rm = TRUE)/sqrt(length(x))
})

#### PSD Analysis ####

rich_bay <- list()
for (i in levels(RichEvenSamp$Bay)) {
    rich_bay[[i]] <- RichEvenSamp[RichEvenSamp$Bay == i, c("t", "S_detrend")]
}

rich_bay_detrend <- lapply(rich_bay, FUN = function(x) {
    mod <- lm(x$S_detrend ~ x$t)
    data.frame(t = x$t - min(x$t), rich = resid(mod))
})


rich_bay_spec <- lapply(rich_bay_detrend, FUN = function(x) spec.lomb(x$t, x$rich, f = seq(0.01, 2, 0.01), mode = "generalized"))
rich_sig_per <- lapply(rich_bay_spec, FUN = function(x) {
    f_PSD_p <- data.frame(f = x$f, A = x$PSD, p = x$p)
    f_PSD_p_sig <- f_PSD_p[f_PSD_p$p < 0.05, ]
    f_sig <- f_PSD_p_sig$f
    per_sig <- 1/f_sig
    return(per_sig)
})

even_bay <- list()
for (i in levels(RichEvenSamp$Bay)) {
    even_bay[[i]] <- RichEvenSamp[RichEvenSamp$Bay == i, c("t", "J_detrend")]
    even_bay[[i]] <- even_bay[[i]][is.na(even_bay[[i]]$J_detrend) == FALSE, ]
}

even_bay_detrend <- lapply(even_bay, FUN = function(x) {
    mod <- lm(x$J_detrend ~ x$t)
    data.frame(t = x$t - min(x$t), even = resid(mod))
})


even_bay_spec <- lapply(even_bay_detrend, FUN = function(x) spec.lomb(x$t, x$even, f = seq(0.01, 2, 0.01), mode = "generalized"))
even_sig_per <- lapply(even_bay_spec, FUN = function(x) {
    f_PSD_p <- data.frame(f = x$f, A = x$PSD, p = x$p)
    f_PSD_p_sig <- f_PSD_p[f_PSD_p$p < 0.05, ]
    f_sig <- f_PSD_p_sig$f
    per_sig <- 1/f_sig
    return(per_sig)
})



#### Plots ####
cbPalette <- brewer.pal(4, "Dark2")

rich_time <- data.frame(t = rich_even_bay_timeseries_means$t, mean = rich_even_bay_timeseries_means$S_detrend, se = rich_even_bay_timeseries_se$S_detrend, bay = rich_even_bay_timeseries_means$Bay)
rich_time$bay <- factor(rich_time$bay, levels = c("AP", "CK", "TB", "CH"))
even_time <- data.frame(t = rich_even_bay_timeseries_means$t, mean = rich_even_bay_timeseries_means$J_detrend, se = rich_even_bay_timeseries_se$J_detrend, bay = rich_even_bay_timeseries_means$Bay)
even_time <- even_time[-148, ]
even_time$bay <- factor(even_time$bay, levels = c("AP", "CK", "TB", "CH"))

# Plot TS

rich_timePlot <- ggplot(rich_time) + geom_ribbon(aes(x = t, ymin = mean - se, ymax = mean + se, fill = bay), alpha = 0.3, show.legend = FALSE) + scale_fill_manual(values = cbPalette) + 
    geom_line(aes(x = t, y = mean, colour = bay), size = 0.2) + scale_x_continuous(limits = c(1998, 2018), breaks = c(1998:2018), expand = c(0.005, 0.005)) + theme_bw() + theme(axis.text.x = element_blank(), 
    legend.key.width = unit(5, "lines"), axis.title.y = element_text(vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = rel(0.8)), legend.key.height = unit(4, "lines")) + 
    labs(title = NULL, x = NULL, y = "\nRichness", colour = NULL) + scale_color_manual(values = cbPalette, labels = c("AB", "CK", "TB", "CH")) + guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Rich_time.tiff", rich_timePlot, width = 4, height = 0.8, dpi = 600)

even_timePlot <- ggplot(even_time) + geom_ribbon(aes(x = t, ymin = mean - se, ymax = mean + se, fill = bay), alpha = 0.3, show.legend = FALSE) + scale_fill_manual(values = cbPalette) + 
    geom_line(aes(x = t, y = mean, colour = bay), size = 0.2, show.legend = FALSE) + scale_x_continuous(limits = c(1998, 2018), breaks = c(1998:2018), expand = c(0.005, 0.005)) + theme_bw() + 
    theme(axis.text.x = element_blank(), axis.title.y = element_text(vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = rel(0.8)), legend.key.width = unit(5, "lines"), legend.key.height = unit(4, 
        "lines")) + labs(title = NULL, x = NULL, y = "\nEvenness", colour = NULL) + scale_color_manual(values = cbPalette, labels = c("AB", "CK", "TB", "CH")) + guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Even_time.tiff", even_timePlot, width = 4, height = 0.8, dpi = 600)

# Plot specs

rich <- rich_bay_spec
rich_per_PSD <- do.call(rbind, lapply(seq_along(rich), FUN = function(j) {
    data.frame(bay = names(rich)[[j]], per = 1/rich[[j]]$f, psd = rich[[j]]$PSD, p = rich[[j]]$p)
}))
rich_per_sig <- rich_per_PSD$per
rich_psd_sig <- rich_per_PSD$psd
rich_psd_sig[rich_per_PSD$p > 0.05] <- -200
rich_per_PSD_sig <- data.frame(bay = rich_per_PSD$bay, rich_per_sig, rich_psd_sig)
rich_per_PSD$bay <- factor(rich_per_PSD$bay, levels = c("AP", "CK", "TB", "CH"))
rich_per_PSD_sig$bay <- factor(rich_per_PSD_sig$bay, levels = c("AP", "CK", "TB", "CH"))

rich_spec_plot <- ggplot(rich_per_PSD) + geom_ribbon(dat = rich_per_PSD_sig, aes(x = rich_per_sig, ymin = 0, ymax = rich_psd_sig, fill = bay), alpha = 0.5, show.legend = FALSE) + geom_line(aes(x = per, 
    y = psd, colour = bay), size = 0.2, show.legend = FALSE) + scale_fill_manual(values = cbPalette) + scale_x_log10(limits = c(0.5, 10), breaks = c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
    10), expand = c(0.001, 0.001)) + scale_y_continuous(limits = c(0, NA)) + theme_bw() + theme(axis.text.x = element_blank(), legend.key.width = unit(5, "lines"), axis.text.y = element_text(size = rel(0.8)), 
    axis.title.y = element_text(vjust = 0.5, hjust = 0.5), legend.key.height = unit(4, "lines")) + labs(title = NULL, x = NULL, y = "\nPSD", colour = NULL) + scale_color_manual(values = cbPalette, 
    labels = c("AB", "CK", "TB", "CH")) + guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Rich_spec.tiff", rich_spec_plot, width = 4, height = 0.8, dpi = 600)

even <- even_bay_spec
even_per_PSD <- do.call(rbind, lapply(seq_along(even), FUN = function(j) {
    data.frame(bay = names(even)[[j]], per = 1/even[[j]]$f, psd = even[[j]]$PSD, p = even[[j]]$p)
}))
even_per_sig <- even_per_PSD$per
even_psd_sig <- even_per_PSD$psd
even_psd_sig[even_per_PSD$p > 0.05] <- -200
even_per_PSD_sig <- data.frame(bay = even_per_PSD$bay, even_per_sig, even_psd_sig)
even_per_PSD$bay <- factor(even_per_PSD$bay, levels = c("AP", "CK", "TB", "CH"))
even_per_PSD_sig$bay <- factor(even_per_PSD_sig$bay, levels = c("AP", "CK", "TB", "CH"))

even_spec_plot <- ggplot(even_per_PSD) + geom_ribbon(dat = even_per_PSD_sig, aes(x = even_per_sig, ymin = 0, ymax = even_psd_sig, fill = bay), alpha = 0.5, show.legend = FALSE) + geom_line(aes(x = per, 
    y = psd, colour = bay), size = 0.2, show.legend = FALSE) + scale_fill_manual(values = cbPalette) + scale_x_log10(limits = c(0.5, 10), breaks = c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
    10), expand = c(0.001, 0.001)) + scale_y_continuous(limits = c(0, NA)) + theme_bw() + theme(axis.text.x = element_blank(), legend.key.width = unit(5, "lines"), axis.text.y = element_text(size = rel(0.8)), 
    axis.title.y = element_text(vjust = 0.5, hjust = 0.5), legend.key.height = unit(4, "lines")) + labs(title = NULL, x = NULL, y = "\nPSD", colour = NULL) + scale_color_manual(values = cbPalette, 
    labels = c("AB", "CK", "TB", "CH")) + guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Even_spec.tiff", even_spec_plot, width = 4, height = 0.8, dpi = 600)

#### Synchrony ####

ts_sync <- list(inter = list(), intra = list())
pdf(file = "allDiagPlots.pdf")
table_inter_RichEven <- data.frame()
table_intra_RichEven <- data.frame()

rich_bay_time <- list()
for (i in levels(rich_time$bay)) {
    rich_bay_time[[i]] <- rich_time[rich_time$bay == i, c("t", "mean")]
}

even_bay_time <- list()
for (i in levels(even_time$bay)) {
    even_bay_time[[i]] <- even_time[even_time$bay == i, c("t", "mean")]
    even_bay_time[[i]] <- even_bay_time[[i]][is.na(even_bay_time[[i]]$mean) == FALSE, ]
}

series_rich <- lapply(rich_bay_time, FUN = function(x) {
    timeser <- ts(x$mean, frequency = 4, start = min(x$t))
    timeser_comps <- decompose(timeser)
    interTS <- x$mean - timeser_comps$seasonal
    intraTS <- ts(timeser_comps$seasonal, frequency = 4, start = min(x$t))
    out <- list(interTS = interTS, intraTS = intraTS)
    out
})

series_even <- lapply(even_bay_time, FUN = function(x) {
    timeser <- ts(x$mean, frequency = 4, start = min(x$t))
    timeser_comps <- decompose(timeser)
    interTS <- x$mean - timeser_comps$seasonal
    intraTS <- ts(timeser_comps$seasonal, frequency = 4, start = min(x$t))
    out <- list(interTS = interTS, intraTS = intraTS)
    out
})

compList <- c("AP_CK", "AP_TB", "AP_CH", "CK_CH", "CK_TB", "TB_CH")
table_inter <- data.frame()
table_intra <- data.frame()

for (j in c("AP", "CK", "TB", "CH")) {
    for (k in c("AP", "CK", "TB", "CH")) {
        if (paste(j, k, sep = "_") %in% compList) {
            corr_inter <- ccf(series_rich[[j]][["interTS"]], series_rich[[k]][["interTS"]], na.action = na.pass, main = paste(j, "_", k, " Rich inter", sep = ""), lag.max = 2)
            max_inter <- corr_inter$acf[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
            lag_inter <- corr_inter$lag[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
            CI_inter <- qnorm((1 + 0.95)/2)/sqrt(corr_inter$n.used)
            sig_inter <- ifelse(abs(max_inter) > CI_inter, 1, 0)
            out_inter <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_inter), lag = as.numeric(lag_inter), sig = as.numeric(sig_inter))
            ts_sync$inter$Rich[[paste(j, k, sep = "_")]][["corr"]] <- corr_inter
            ts_sync$inter$Rich[[paste(j, k, sep = "_")]][["table"]] <- out_inter
            
            table_inter <- rbind(table_inter, cbind(stat = "rich", out_inter))
            
            corr_inter <- ccf(series_even[[j]][["interTS"]], series_even[[k]][["interTS"]], na.action = na.pass, main = paste(j, "_", k, " Even inter", sep = ""), lag.max = 2)
            max_inter <- corr_inter$acf[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
            lag_inter <- corr_inter$lag[which(abs(corr_inter$acf) == max(abs(corr_inter$acf)))]
            CI_inter <- qnorm((1 + 0.95)/2)/sqrt(corr_inter$n.used)
            sig_inter <- ifelse(abs(max_inter) > CI_inter, 1, 0)
            out_inter <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_inter), lag = as.numeric(lag_inter), sig = as.numeric(sig_inter))
            ts_sync$inter$Even[[paste(j, k, sep = "_")]][["corr"]] <- corr_inter
            ts_sync$inter$Even[[paste(j, k, sep = "_")]][["table"]] <- out_inter
            
            table_inter <- rbind(table_inter, cbind(stat = "even", out_inter))
            
            
            corr_intra <- ccf(series_rich[[j]][["intraTS"]], series_rich[[k]][["intraTS"]], na.action = na.pass, main = paste(j, "_", k, " Even intra", sep = ""), lag.max = 1)
            max_intra <- corr_intra$acf[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
            lag_intra <- corr_intra$lag[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
            CI_intra <- qnorm((1 + 0.95)/2)/sqrt(corr_intra$n.used)
            sig_intra <- ifelse(abs(max_intra) > CI_intra, 1, 0)
            out_intra <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_intra), lag = as.numeric(lag_intra), sig = as.numeric(sig_intra))
            ts_sync$intra$Rich[[paste(j, k, sep = "_")]][["corr"]] <- corr_intra
            ts_sync$intra$Rich[[paste(j, k, sep = "_")]][["table"]] <- out_intra
            
            table_intra <- rbind(table_intra, cbind(stat = "rich", out_intra))
            
            corr_intra <- ccf(series_even[[j]][["intraTS"]], series_even[[k]][["intraTS"]], na.action = na.pass, main = paste(j, "_", k, " Even intra", sep = ""), lag.max = 0)
            max_intra <- corr_intra$acf[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
            lag_intra <- corr_intra$lag[which(abs(corr_intra$acf) == max(abs(corr_intra$acf)))]
            CI_intra <- qnorm((1 + 0.95)/2)/sqrt(corr_intra$n.used)
            sig_intra <- ifelse(abs(max_intra) > CI_intra, 1, 0)
            out_intra <- cbind(bay1 = j, bay2 = k, maxcorr = as.numeric(max_intra), lag = as.numeric(lag_intra), sig = as.numeric(sig_intra))
            ts_sync$intra$Even[[paste(j, k, sep = "_")]][["corr"]] <- corr_intra
            ts_sync$intra$Even[[paste(j, k, sep = "_")]][["table"]] <- out_intra
            
            table_intra <- rbind(table_intra, cbind(stat = "even", out_intra))
            
        }
    }
}
# ts_sync$inter[['table']] = table_inter ts_sync$intra[['table']] = table_intra
table_inter_RichEven <- rbind(table_inter_RichEven, cbind(type = "inter", table_inter))
table_intra_RichEven <- rbind(table_intra_RichEven, cbind(type = "intra", table_intra))

dev.off()
ts_sync$inter[["table"]] <- table_inter_RichEven
ts_sync$intra[["table"]] <- table_intra_RichEven
ts_sync[["table"]] <- rbind(table_inter_RichEven, table_intra_RichEven)
ts_sync$table$type <- as.factor(ts_sync$table$type)
ts_sync$table$stat <- as.factor(ts_sync$table$stat)
ts_sync$table$bay1 <- factor(ts_sync$table$bay1, levels = c("AP", "CK", "TB", "CH"))
ts_sync$table$bay2 <- factor(ts_sync$table$bay2, levels = c("AP", "CK", "TB", "CH"))
ts_sync$table$maxcorr <- as.numeric(ts_sync$table$maxcorr)
ts_sync$table$lag <- as.numeric(ts_sync$table$lag)
ts_sync$table$sig <- as.numeric(ts_sync$table$sig)
ts_sync$table
cast_table <- recast(ts_sync$table, formula = type + stat + bay2 ~ variable + bay1, measure.var = c("maxcorr", "lag", "sig"), fun.aggregate = mean)

save(rich_timePlot, rich_spec_plot, even_timePlot, even_spec_plot, file = "RichEven_Plots.RData")
