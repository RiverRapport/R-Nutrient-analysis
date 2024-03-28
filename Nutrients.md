Nutrient Indicator - Great River Rapport
================
Erin Smith
2024-03-28

    ## Warning: package 'knitr' was built under R version 4.2.3

## Load packages

``` r
library(Kendall)
library(corrplot)
library(trend)
library(segmented)
library(modifiedmk)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(EnvStats)
library(lubridate)
library(lmtest)
library(dplyr)
library(tidyr)
library(openair)
library(mgcv)
library(mgcViz)
library(gratia)
library(itsadug)
library(FSA)
library(ggsignif)
library(viridis)
library(RColorBrewer)
library(forecast)
library(gganimate)
library(gifski)
library(nlme)
library(ggmap)
library(Hmisc)
library(formatR)
library(FSA)
```

## Load dataframes and format data

``` r
GN <- read.csv("Gov_nut_data_depth.csv", header = TRUE)

AD <- read.csv("Ac_nut_data_depth.csv", header = TRUE)

clim <- read.csv("Monthly_climate_USLR_1960.2022.csv", header = TRUE)
clim$River_section <- factor(clim$River_section, levels = c("TI", "BR", "LSL", "CA",
    "LSF"))

Ac_gov <- rbind(GN, AD)
Ac_gov$Sample_date <- as.Date(Ac_gov$Sample_date, "%Y-%m-%d")
Ac_gov$River_section <- factor(Ac_gov$River_section, levels = c("TI", "BR", "LSL",
    "CA", "LSF"))
Ac_gov$Site <- factor(Ac_gov$Site)
Ac_gov$Data_source <- factor(Ac_gov$Data_source)
Ac_gov$year <- year(ymd(Ac_gov$Sample_date))
Ac_gov$month <- month(ymd(Ac_gov$Sample_date))
Ac_gov$day <- day(ymd(Ac_gov$Sample_date))
Ac_gov$DOY <- yday(Ac_gov$Sample_date)
Ac_gov$NS_MC <- factor(Ac_gov$NS_MC)
Ac_gov$NP <- with(Ac_gov, TN/TP)

# add climate data
Nut_precip <- merge(x = Ac_gov, y = clim, by = c("River_section", "year", "month"),
    all = FALSE)
# add season to df
Nut_precip <- Nut_precip %>%
    mutate(season = case_when(month > 2 & month < 6 ~ "spring", month > 5 & month <
        9 ~ "summer", month > 8 & month < 12 ~ "fall", month > 11 | month < 3 ~ "winter"),
        cont_year = decimal_date(Sample_date))

write.csv(Nut_precip, "USLR_nutrients.csv", row.names = FALSE)

Meta_data <- Nut_precip %>%
    group_by(Data_source) %>%
    summarise(Yearmin = min(year), Yearmax = max(year), RS = unique(River_section))

write.csv(Meta_data, "RR_Nutrients_Metadata.csv", row.names = FALSE)

GLIP_USGS <- subset(Nut_precip, Data_source == "OMECP_GLIP" | Data_source == "USGS-NY")

GLIP_Ann <- GLIP_USGS %>%
    group_by(year, River_section) %>%
    summarise(TP.mg.L = mean(TP, na.rm = T)) %>%
    pivot_wider(names_from = River_section, values_from = TP.mg.L)

write.csv(GLIP_Ann, "TP_TI.BR.CA_Annual.csv", row.names = F)
```

## Nutrient Analysis

### Correlation Matrix

``` r
res <- cor(Nut_precip[, c(15, 16, 18, 23, 24, 20)], use = "pairwise.complete.obs")
corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
```

![](Nutrients_files/figure-gfm/Correlation%20matrix-1.png)<!-- -->

``` r
Avg <- Nut_precip %>%
    mutate(era = ifelse(year < 1981, "1966-1980", if_else(year < 1991, "1981-1990",
        if_else(year < 2001, "1991-2000", if_else(year < 2011, "2001-2010", "2011-2022"))))) %>%
    group_by(River_section, era) %>%
    summarise(TPx = mean(TP, na.rm = TRUE), TPsd = sd(TP, na.rm = TRUE), TNx = mean(TN,
        na.rm = TRUE), TNsd = sd(TN, na.rm = TRUE))
write.csv(Avg, "Avg_TP_USLR.csv", row.names = F)
```

### Nutrients by Depth

``` r
TP_aov <- aov(TP ~ NS_MC, data = Nut_precip)
par(mfrow = c(2, 2))
plot(TP_aov)
```

![](Nutrients_files/figure-gfm/Kruskal%20Wallis%20for%20depth-1.png)<!-- -->

``` r
kruskal.test(TN ~ NS_MC, data = Nut_precip)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  TN by NS_MC
    ## Kruskal-Wallis chi-squared = 181.57, df = 1, p-value < 2.2e-16

``` r
p1 <- ggplot(data = Nut_precip, aes(NS_MC, TN, fill = NS_MC)) + geom_boxplot() +
    xlab("Depth") + ylab("TN (mg*L-1)") + scale_fill_viridis(discrete = TRUE) + scale_x_discrete(labels = c("Main Channel",
    "Nearshore")) + theme(axis.title = element_text(size = 15), panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 12), legend.position = "none") + geom_signif(comparisons = list(c("MC",
    "NS")), map_signif_level = TRUE)

ggsave("TN_depth_boxplot.svg", p1, dpi = 600, width = 10, height = 8)

p2 <- ggplot(data = Nut_precip, aes(NS_MC, TP, fill = NS_MC)) + geom_boxplot() +
    xlab("Depth") + ylab("TP (ug*L-1)") + scale_fill_viridis(discrete = TRUE) + scale_x_discrete(labels = c("Main Channel",
    "Nearshore")) + theme(axis.title = element_text(size = 15), panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 12), legend.position = "none") + geom_signif(comparisons = list(c("MC",
    "NS")), map_signif_level = TRUE)

ggsave("TP_depth_boxplot.svg", p2, dpi = 600, width = 10, height = 8)

ggarrange(p1, p2, ncol = 2)
```

![](Nutrients_files/figure-gfm/Kruskal%20Wallis%20for%20depth-2.png)<!-- -->
\### Nutrients by River Section

``` r
TP_aov <- aov(TP ~ River_section, data = Nut_precip)
par(mfrow = c(2, 2))
plot(TP_aov)
```

![](Nutrients_files/figure-gfm/Kruskal%20Wallis%20for%20river%20section-1.png)<!-- -->

``` r
TP_RS <- kruskal.test(TP ~ River_section, data = Nut_precip)
TP_RS
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  TP by River_section
    ## Kruskal-Wallis chi-squared = 345.35, df = 4, p-value < 2.2e-16

``` r
dunnTest(TP ~ River_section, data = Nut_precip)
```

    ## Dunn (1964) Kruskal-Wallis multiple comparison

    ##   p-values adjusted with the Holm method.

    ##    Comparison          Z      P.unadj        P.adj
    ## 1     BR - CA  -5.765079 8.161937e-09 4.897162e-08
    ## 2    BR - LSF -14.214809 7.415299e-46 6.673769e-45
    ## 3    CA - LSF  -5.418616 6.006224e-08 2.402490e-07
    ## 4    BR - LSL -10.433194 1.749075e-25 1.399260e-24
    ## 5    CA - LSL  -6.840777 7.876490e-12 5.513543e-11
    ## 6   LSF - LSL  -4.078725 4.528329e-05 1.358499e-04
    ## 7     BR - TI -15.955216 2.620666e-57 2.620666e-56
    ## 8     CA - TI  -3.820681 1.330835e-04 2.661669e-04
    ## 9    LSF - TI   3.375218 7.375713e-04 7.375713e-04
    ## 10   LSL - TI   5.586803 2.312883e-08 1.156442e-07

``` r
p2 <- ggplot(data = Nut_precip, aes(River_section, TP, fill = River_section)) + geom_boxplot() +
    xlab("River Section") + scale_fill_viridis(discrete = TRUE, labels = c("Thousand Islands",
    "Brockville Narrows", "Lake St. Lawrence", "Cornwall-Akwesasne", "Lake St. Francis")) +
    ylab("TP (µg*L-1)") + labs(fill = "River Section") + theme(axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14), axis.text = element_text(size = 12),
    legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))  #+
# geom_signif(comparisons = list(c('BR',
# 'LSF'),c('CA','LSF'),c('BR','TI'),c('LSF','TI')), y_position =
# c(0.4,0.5,0.4,0.55), map_signif_level=TRUE)
p2
```

![](Nutrients_files/figure-gfm/Kruskal%20Wallis%20for%20river%20section-2.png)<!-- -->

## Mann Kendall time series analysis

### Phosphorus

``` r
Seas_Nut <- Nut_precip %>%
    filter(year < 2021) %>%
    group_by(year, season) %>%
    summarise(TP = mean(TP, na.rm = TRUE), TDP = mean(TDP, na.rm = TRUE), SRP = mean(SRP,
        na.rm = TRUE), NH3.4 = mean(NH3.4, na.rm = TRUE), NO2.3 = mean(NO2.3, na.rm = TRUE),
        TN = mean(TN, na.rm = TRUE))

mod_raw <- lm(TP ~ year, data = Seas_Nut)
par(mar = c(3, 3, 3, 0), mfrow = c(1, 2))
acf(resid(mod_raw), plot = T, main = paste("ACF Data_source"), lag.max = 35)
pacf(resid(mod_raw), plot = T, main = paste("pACF Data_source"), lag.max = 35)
```

![](Nutrients_files/figure-gfm/Seasonal%20Mann%20Kendall-1.png)<!-- -->

``` r
SMKresultsP <- NULL
for (parameter in c("TP", "SRP", "TDP")) {
    P_SMK <- kendallSeasonalTrendTest(Seas_Nut[[parameter]] ~ season + year, data = Seas_Nut,
        na.action = na.pass, alternative = "two.sided", independent.obs = TRUE)

    SMKresultsP = rbind(SMKresultsP, data.frame(parameter, data.frame(as.list(P_SMK$estimate),
        as.list(P_SMK$p.value))))
}
write.table(SMKresultsP, "Seasonal Mann Kendall_P.txt", sep = "\t")
SMKresultsP
```

    ##   parameter        tau         slope  intercept Chi.Square..Het.    z..Trend.
    ## 1        TP -0.2444252 -1.533796e-04 0.31552670       0.09633296 2.195878e-07
    ## 2       SRP -0.1692979 -4.979541e-05 0.10696670       0.57989831 8.843957e-04
    ## 3       TDP -0.1195830 -3.846154e-05 0.06333248       0.03659001 3.289438e-02

#### Test for breakpoints

``` r
icefree <- Nut_precip %>%
    filter(month > 4 & month < 11) %>%
    group_by(year) %>%
    summarise(TP = mean(TP, na.rm = TRUE), TDP = mean(TDP, na.rm = TRUE), SRP = mean(SRP,
        na.rm = TRUE), NH3.4 = mean(NH3.4, na.rm = TRUE), NO2.3 = mean(NO2.3, na.rm = TRUE),
        TN = mean(TN, na.rm = TRUE))

for (parameter in c("TP", "SRP")) {
    y <- icefree[[parameter]]
    x <- icefree$year
    break.value <- segmented(lm(y ~ x), seg.Z = ~x)$psi[, 2]
    break.value
    seg.res <- segmented(lm(y ~ x), seg.Z = ~x)
    p1 <- pscore.test(seg.res, seg.Z = ~x, k = 10, alternative = c("two.sided"),
        values = NULL, dispersion = NULL, df.t = NULL, more.break = FALSE, n.break = 1)
    plot(seg.res, conf.level = 0.95, shade = TRUE, xlab = "Year", ylab = paste(parameter,
        "ug*L-1"), ylim = c(0, 0.1), xlim = c(1965, 2022))
    points(icefree[[parameter]] ~ icefree$year, pch = 19)
    lines(icefree[[parameter]] ~ icefree$year)
    lines(seg.res, col = 2, pch = 19, lwd = 2)
    points(seg.res, col = 4, link = FALSE)
    abline(v = break.value, col = "blue")
    print(p1)
}
```

![](Nutrients_files/figure-gfm/Check%20for%20Breakpoints-1.png)<!-- -->

    ## 
    ##  Score test for one/two changes in the slope
    ## 
    ## data:  formula = y ~ x 
    ## breakpoint for variable = x 
    ## model = gaussian , link = identity , method = segmented.lm
    ## observed value = 7.062, n.points = 10, p-value = 2.986e-09
    ## alternative hypothesis: two.sided   (1 breakpoint)

![](Nutrients_files/figure-gfm/Check%20for%20Breakpoints-2.png)<!-- -->

    ## 
    ##  Score test for one/two changes in the slope
    ## 
    ## data:  formula = y ~ x 
    ## breakpoint for variable = x 
    ## model = gaussian , link = identity , method = segmented.lm
    ## observed value = 17.968, n.points = 10, p-value < 2.2e-16
    ## alternative hypothesis: two.sided   (1 breakpoint)

### Plot temporal N:P trend

``` r
Annual_nut <- Nut_precip %>%
    group_by(year, River_section) %>%
    summarise(TP = mean(TP, na.rm = TRUE), NP = mean(NP, na.rm = TRUE), NO3 = mean(NO3,
        na.rm = TRUE), TN = mean(TN, na.rm = TRUE))

NP <- ggplot(Annual_nut, aes(year, NP)) + geom_point() + geom_smooth(method = "gam",
    span = 0.5) + xlab("Year") + ylab("N:P") + geom_hline(yintercept = 29, color = "orange",
    size = 1, linetype = "dashed") + theme_bw()
NP
```

![](Nutrients_files/figure-gfm/N:P%20temporal%20plot-1.png)<!-- -->

``` r
ggsave("Plots/NP_Ratio_SLR.svg", NP, dpi = 600, height = 8, width = 10)

NP_MK <- kendallTrendTest(NP ~ year, data = Annual_nut, na.action = na.pass, alternative = "two.sided")
NP_MK$p.value
```

    ##            z 
    ## 5.376252e-17

``` r
NP_MK$estimate
```

    ##           tau         slope     intercept 
    ##     0.4661262     1.2722364 -2477.3418860

``` r
y <- Annual_nut$NP
x <- Annual_nut$year
break.value <- segmented(lm(y ~ x), seg.Z = ~x, psi = 2000)$psi[, 2]
break.value
```

    ## [1] 2008

``` r
seg.res <- segmented(lm(y ~ x), seg.Z = ~x, psi = 2000)
p1 <- pscore.test(seg.res, seg.Z = ~x, k = 10, alternative = c("two.sided"), values = NULL,
    dispersion = NULL, df.t = NULL, more.break = FALSE, n.break = 1)
plot(seg.res, conf.level = 0.95, shade = TRUE, xlab = "Year", ylab = paste("N:P"),
    xlim = c(1970, 2022))
points(Annual_nut$NP ~ Annual_nut$year, pch = 19)
lines(Annual_nut$NP ~ Annual_nut$year)
lines(seg.res, col = 2, pch = 19, lwd = 2)
points(seg.res, col = 4, link = FALSE)
abline(v = break.value, col = "blue")
```

![](Nutrients_files/figure-gfm/N:P%20temporal%20plot-2.png)<!-- -->

``` r
print(p1)
```

    ## 
    ##  Score test for one/two changes in the slope
    ## 
    ## data:  formula = y ~ x 
    ## breakpoint for variable = x 
    ## model = gaussian , link = identity , method = segmented.lm
    ## observed value = -4.6303, n.points = 10, p-value = 8.125e-06
    ## alternative hypothesis: two.sided   (1 breakpoint)

``` r
NP_gam <- gam(NP ~ s(year, bs = "tp"), method = "REML", family = tw(link = "log"),
    data = Annual_nut)

k.check(NP_gam)
```

    ##         k'      edf   k-index p-value
    ## s(year)  9 5.904057 0.9940024  0.4925

``` r
appraise(NP_gam)
```

![](Nutrients_files/figure-gfm/N:P%20temporal%20plot-3.png)<!-- -->

``` r
par(mfrow = c(1, 2))
acf(resid(NP_gam), lag.max = 36, main = "ACF")
pacf(resid(NP_gam), lag.max = 36, main = "pACF")
```

![](Nutrients_files/figure-gfm/N:P%20temporal%20plot-4.png)<!-- -->

``` r
summary(NP_gam)
```

    ## 
    ## Family: Tweedie(p=1.99) 
    ## Link function: log 
    ## 
    ## Formula:
    ## NP ~ s(year, bs = "tp")
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   4.1798     0.0262   159.5   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##           edf Ref.df     F p-value    
    ## s(year) 5.904  7.084 37.84  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.371   Deviance explained = 62.2%
    ## -REML = 660.83  Scale est. = 0.10524   n = 147

``` r
m2.d <- derivatives(NP_gam, term = "s(year)", n = 200)
m2.d <- m2.d %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))

sm_year <- smooth_estimates(NP_gam, smooth = "s(year)", dist = 0.05)
sm_year$upper <- sm_year$est + 1.96 * sm_year$se
sm_year$lower <- sm_year$est - 1.96 * sm_year$se

TP_year_plot <- ggplot(sm_year, aes(x = year, y = est)) + geom_ribbon(alpha = 0.5,
    aes(ymin = lower, ymax = upper, fill = "confidence interval")) + annotate("rect",
    xmin = c(2009.543), xmax = c(2016.613), ymin = -Inf, ymax = Inf, fill = "purple",
    alpha = 0.2) + annotate("rect", xmin = c(1967.121, 1992.035), xmax = c(1987.322,
    2004.492), ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.2) + geom_line(data = sm_year,
    aes(x = year, y = est, color = "GAM"), linewidth = 1) + scale_fill_manual(values = "lightblue",
    name = NULL) + scale_color_manual(values = "blue", name = NULL) + ylab("Partial effect of year on N:P") +
    scale_x_continuous(name = "Year", breaks = c(1970, 1980, 1990, 2000, 2010, 2020)) +
    theme_minimal(base_size = 16) + theme(legend.position = "none")
TP_year_plot
```

![](Nutrients_files/figure-gfm/N:P%20temporal%20plot-5.png)<!-- -->

### Nitrogen Mann Kendall

``` r
SMKresultsN <- NULL
for (parameter in c("TN", "NO2.3", "NH3.4")) {
    N_SMK <- kendallSeasonalTrendTest(Seas_Nut[[parameter]] ~ season + year, data = Seas_Nut,
        na.action = na.pass, alternative = "two.sided", independent.obs = TRUE)

    SMKresultsN = rbind(SMKresultsN, data.frame(parameter, data.frame(as.list(N_SMK$estimate),
        as.list(N_SMK$p.value))))
}
write.table(SMKresultsN, "Seasonal Mann Kendall_N.txt", sep = "\t")
SMKresultsN
```

    ##   parameter         tau         slope   intercept Chi.Square..Het.    z..Trend.
    ## 1        TN  0.41324310  4.324612e-03  -8.7535220       0.46081714 1.778177e-16
    ## 2     NO2.3  0.63481986  5.590093e-03 -11.4718456       0.69755077 5.298653e-40
    ## 3     NH3.4 -0.07863211 -9.198677e-05   0.3478187       0.01514836 7.782626e-02

### pscore test to test for breakpoints

``` r
for (parameter in c("TN", "NO2.3", "NH3.4")) {
    y <- icefree[[parameter]]
    x <- icefree$year
    break.value <- segmented(lm(y ~ x), seg.Z = ~x, psi = 2000)$psi[, 2]
    break.value
    seg.res <- segmented(lm(y ~ x), seg.Z = ~x, psi = 2000)
    p1 <- pscore.test(seg.res, seg.Z = ~x, k = 10, alternative = c("two.sided"),
        values = NULL, dispersion = NULL, df.t = NULL, more.break = FALSE, n.break = 1)
    plot(seg.res, conf.level = 0.95, shade = TRUE, xlab = "Year", ylab = paste(parameter,
        "mg*L-1"), ylim = c(0, 1), xlim = c(1960, 2022))
    points(icefree[[parameter]] ~ icefree$year, pch = 19)
    lines(icefree[[parameter]] ~ icefree$year)
    lines(seg.res, col = 2, pch = 19, lwd = 2)
    points(seg.res, col = 4, link = FALSE)
    abline(v = break.value, col = "blue")
    print(p1)
}
```

![](Nutrients_files/figure-gfm/N%20breakpoints-1.png)<!-- -->

    ## 
    ##  Score test for one/two changes in the slope
    ## 
    ## data:  formula = y ~ x 
    ## breakpoint for variable = x 
    ## model = gaussian , link = identity , method = segmented.lm
    ## observed value = 0.53622, n.points = 10, p-value = 0.5944
    ## alternative hypothesis: two.sided   (1 breakpoint)

![](Nutrients_files/figure-gfm/N%20breakpoints-2.png)<!-- -->

    ## 
    ##  Score test for one/two changes in the slope
    ## 
    ## data:  formula = y ~ x 
    ## breakpoint for variable = x 
    ## model = gaussian , link = identity , method = segmented.lm
    ## observed value = -1.7114, n.points = 10, p-value = 0.09297
    ## alternative hypothesis: two.sided   (1 breakpoint)

![](Nutrients_files/figure-gfm/N%20breakpoints-3.png)<!-- -->

    ## 
    ##  Score test for one/two changes in the slope
    ## 
    ## data:  formula = y ~ x 
    ## breakpoint for variable = x 
    ## model = gaussian , link = identity , method = segmented.lm
    ## observed value = 6.5509, n.points = 10, p-value = 2.38e-08
    ## alternative hypothesis: two.sided   (1 breakpoint)

### Plot temporal N trends

``` r
par(mfrow = c(1, 3))
for (parameter in c("TN", "NO2.3", "NH3.4")) {
    print(ggplot(icefree, aes(year, icefree[[parameter]])) + geom_point() + geom_smooth(method = "loess",
        span = 0.5) + xlab("Year") + ylab(parameter) + theme_bw())
}
```

![](Nutrients_files/figure-gfm/Seasonal%20N%20plots-1.png)<!-- -->![](Nutrients_files/figure-gfm/Seasonal%20N%20plots-2.png)<!-- -->![](Nutrients_files/figure-gfm/Seasonal%20N%20plots-3.png)<!-- -->

## Let’s Try GAMs

### Data exploration

``` r
for (parameter in c("TP", "TN")) {
    hist(Nut_precip[[parameter]], main = paste("Histogram of", parameter), xlab = parameter)
}
```

![](Nutrients_files/figure-gfm/Data%20exploration%20for%20GAMs-1.png)<!-- -->![](Nutrients_files/figure-gfm/Data%20exploration%20for%20GAMs-2.png)<!-- -->

``` r
# all parameters left skewed
```

### Need to check for autocorrelation

``` r
# try models for different subsets of the data
spring <- Nut_precip %>%
    filter(season == "spring")

# check on basic model
M <- list(c(1, 0.5), NA)
TPM1 <- gam(TP ~ s(year, m = 2) + s(year, by = River_section, m = 1) + s(River_section,
    bs = "re") + s(Data_source, bs = "re") + s(DOY, bs = "tp") + ti(year, DOY, bs = "tp",
    m = c(1, 0.5)) + s(Longitude, Latitude, bs = "ds", m = c(1, 0.5), k = 10) + ti(Longitude,
    Latitude, year, d = c(2, 1), bs = c("ds", "tp"), m = c(1, 0.5)) + ti(Longitude,
    Latitude, DOY, d = c(2, 1), bs = c("ds", "tp"), m = c(1, 0.5)) + s(NS_MC, bs = "re"),
    data = spring, method = "REML", family = tw(link = "log"))
k.check(TPM1)
```

    ##                             k'          edf   k-index p-value
    ## s(year)                      9  1.004553180 0.9517196  0.1850
    ## s(year):River_sectionTI      8  5.605691183 0.9517196  0.1750
    ## s(year):River_sectionBR      8  5.145864251 0.9517196  0.2075
    ## s(year):River_sectionLSL     8  0.003413749 0.9517196  0.1850
    ## s(year):River_sectionCA      8  0.013307279 0.9517196  0.1650
    ## s(year):River_sectionLSF     8  4.394017078 0.9517196  0.1925
    ## s(River_section)             5  0.011136485        NA      NA
    ## s(Data_source)              17 11.246891228        NA      NA
    ## s(DOY)                       9  7.168558512 0.9475017  0.1250
    ## ti(year,DOY)                16  5.594209048 0.8324931  0.0000
    ## s(Longitude,Latitude)        9  1.794300779 0.9607488  0.4175
    ## ti(Longitude,Latitude,year) 79 28.847427931 0.8415672  0.0000
    ## ti(Longitude,Latitude,DOY)  71 14.917526711 0.8418798  0.0000
    ## s(NS_MC)                     2  0.913434568        NA      NA

``` r
appraise(TPM1)
```

![](Nutrients_files/figure-gfm/Subset%20and%20autocorrelation%20checks-1.png)<!-- -->

``` r
par(mfrow = c(1, 2))
acf(resid(TPM1), lag.max = 36, main = "ACF")
pacf(resid(TPM1), lag.max = 36, main = "pACF")
```

![](Nutrients_files/figure-gfm/Subset%20and%20autocorrelation%20checks-2.png)<!-- -->

``` r
gam_errors <- residuals(TPM1, type = "response")
error_mod <- auto.arima(gam_errors)
error_mod
```

    ## Series: gam_errors 
    ## ARIMA(2,0,2) with zero mean 
    ## 
    ## Coefficients:
    ##           ar1      ar2     ma1     ma2
    ##       -0.8621  -0.5736  0.8377  0.6267
    ## s.e.   0.2982   0.2280  0.2874  0.2001
    ## 
    ## sigma^2 = 0.0002954:  log likelihood = 6661.36
    ## AIC=-13312.71   AICc=-13312.69   BIC=-13283.55

``` r
# check recommended order for AR (p) and/or MA (q)
checkresiduals(error_mod, theme = theme_bw())
```

![](Nutrients_files/figure-gfm/Subset%20and%20autocorrelation%20checks-3.png)<!-- -->

    ## 
    ##  Ljung-Box test
    ## 
    ## data:  Residuals from ARIMA(2,0,2) with zero mean
    ## Q* = 19.903, df = 6, p-value = 0.002882
    ## 
    ## Model df: 4.   Total lags used: 10

``` r
# Make base model with gamm
TPM1_ARMA0 <- gamm(TP ~ s(year, k = 10, bs = "tp") + s(year, River_section, bs = "fs") +
    s(River_section, bs = "re") + s(Data_source, bs = "re") + s(DOY, bs = "tp") +
    ti(year, DOY, bs = "tp", m = c(1, 0.5)) + s(Longitude, Latitude, bs = "ds", m = c(1,
    0.5), k = 15) + ti(Longitude, Latitude, year, d = c(2, 1), bs = c("ds", "tp"),
    m = c(1, 0.5)) + ti(Longitude, Latitude, DOY, d = c(2, 1), bs = c("ds", "tp"),
    m = c(1, 0.5)) + s(NS_MC, bs = "re"), data = spring, family = "tw")
```

    ## 
    ##  Maximum number of PQL iterations:  20

``` r
arma_res <- auto.arima(resid(TPM1_ARMA0$lme, type = "normalized"), stationary = TRUE,
    seasonal = FALSE)
arma_res
```

    ## Series: resid(TPM1_ARMA0$lme, type = "normalized") 
    ## ARIMA(1,0,1) with zero mean 
    ## 
    ## Coefficients:
    ##          ar1      ma1
    ##       0.9497  -0.9123
    ## s.e.  0.0165   0.0208
    ## 
    ## sigma^2 = 0.9735:  log likelihood = -3538.17
    ## AIC=7082.35   AICc=7082.36   BIC=7099.84

``` r
# indicates summer model should be p=1, q=1

# add autoregressive moving average (ARMA)
TPM1_ARMA <- gamm(TP ~ s(year, k = 10, bs = "tp") + s(year, River_section, bs = "fs") +
    s(River_section, bs = "re") + s(Data_source, bs = "re") + s(DOY, bs = "tp") +
    ti(year, DOY, bs = "tp", m = c(1, 0.5)) + s(Longitude, Latitude, bs = "ds", m = c(1,
    0.5), k = 15) + ti(Longitude, Latitude, year, d = c(2, 1), bs = c("ds", "tp"),
    m = c(1, 0.5)) + s(NS_MC, bs = "re"), data = spring, correlation = corARMA(form = ~1 |
    year, p = 1), family = "tw")
```

    ## 
    ##  Maximum number of PQL iterations:  20

``` r
appraise(TPM1_ARMA$gam)
```

![](Nutrients_files/figure-gfm/Subset%20and%20autocorrelation%20checks-4.png)<!-- -->

``` r
par(mfrow = c(1, 2))
res <- resid(TPM1_ARMA$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
```

![](Nutrients_files/figure-gfm/Subset%20and%20autocorrelation%20checks-5.png)<!-- -->

``` r
# ACF is worse with autocorrelation terms - use original model
```

### Build TP GAM model

``` r
TPBM <- gam(TP ~ s(year, m = 2) + s(year, by = River_section, m = 1) + s(River_section,
    bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") + s(Longitude, Latitude,
    bs = "ds", m = c(1, 0.5), k = 20) + s(Data_source, bs = "re"), data = spring,
    method = "REML", family = tw(link = "log"))
gam.check(TPBM)
```

![](Nutrients_files/figure-gfm/TP%20GAM%20building-1.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 13 iterations.
    ## Gradient range [-0.001949752,0.001435482]
    ## (score -8953.368 & scale 0.1752608).
    ## Hessian positive definite, eigenvalue range [5.281127e-05,1926.455].
    ## Model rank =  102 / 102 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                                k'      edf k-index p-value  
    ## s(year)                  9.00e+00 1.01e+00    0.94   0.040 *
    ## s(year):River_sectionTI  8.00e+00 6.71e+00    0.94   0.050 *
    ## s(year):River_sectionBR  8.00e+00 7.11e+00    0.94   0.045 *
    ## s(year):River_sectionLSL 8.00e+00 2.89e-04    0.94   0.050 *
    ## s(year):River_sectionCA  8.00e+00 1.75e-04    0.94   0.075 .
    ## s(year):River_sectionLSF 8.00e+00 4.59e+00    0.94   0.055 .
    ## s(River_section)         5.00e+00 8.26e-04      NA      NA  
    ## s(DOY)                   9.00e+00 7.08e+00    0.93   0.015 *
    ## s(NS_MC)                 2.00e+00 8.69e-01      NA      NA  
    ## s(Longitude,Latitude)    1.90e+01 1.52e+01    0.96   0.340  
    ## s(Data_source)           1.70e+01 1.16e+01      NA      NA  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TPM1 <- gam(TP ~ s(year, m = 2) + s(year, by = River_section, m = 1) + s(River_section,
    bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") + s(Longitude, Latitude,
    bs = "ds", m = c(1, 0.5), k = 20) + s(Data_source, bs = "re") + ti(Longitude,
    Latitude, year, d = c(2, 1), bs = c("ds", "tp"), m = list(c(1, 0.5)), k = 20),
    data = spring, method = "REML", family = tw(link = "log"))
gam.check(TPM1)
```

![](Nutrients_files/figure-gfm/TP%20GAM%20building-2.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 12 iterations.
    ## Gradient range [-0.000860771,0.00373658]
    ## (score -8969.181 & scale 0.165598).
    ## Hessian positive definite, eigenvalue range [0.0001580858,1921.531].
    ## Model rank =  255 / 255 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                                   k'      edf k-index p-value    
    ## s(year)                     9.00e+00 2.27e+00    0.95   0.160    
    ## s(year):River_sectionTI     8.00e+00 4.25e+00    0.95   0.120    
    ## s(year):River_sectionBR     8.00e+00 3.44e-03    0.95   0.120    
    ## s(year):River_sectionLSL    8.00e+00 2.03e-03    0.95   0.130    
    ## s(year):River_sectionCA     8.00e+00 5.79e-03    0.95   0.135    
    ## s(year):River_sectionLSF    8.00e+00 9.13e-01    0.95   0.120    
    ## s(River_section)            5.00e+00 4.39e-04      NA      NA    
    ## s(DOY)                      9.00e+00 1.00e+00    0.93   0.025 *  
    ## s(NS_MC)                    2.00e+00 7.68e-01      NA      NA    
    ## s(Longitude,Latitude)       1.90e+01 1.48e+01    0.97   0.520    
    ## s(Data_source)              1.70e+01 1.16e+01      NA      NA    
    ## ti(Longitude,Latitude,year) 1.53e+02 5.57e+01    0.86  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TPM2 <- gam(TP ~ s(year, m = 2) + s(year, by = River_section, m = 1) + s(River_section,
    bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") + s(Longitude, Latitude,
    bs = "ds", m = c(1, 0.5), k = 20) + s(Data_source, bs = "re") + ti(Longitude,
    Latitude, year, d = c(2, 1), bs = c("ds", "tp"), m = list(c(1, 0.5)), k = 20) +
    ti(Longitude, Latitude, DOY, d = c(2, 1), bs = c("ds", "tp"), m = list(c(1, 0.5)),
        k = 20), data = spring, method = "REML", family = tw(link = "log"))
gam.check(TPM2)
```

![](Nutrients_files/figure-gfm/TP%20GAM%20building-3.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 13 iterations.
    ## Gradient range [-0.002212052,0.004896142]
    ## (score -8992.143 & scale 0.1573139).
    ## Hessian positive definite, eigenvalue range [0.0001013637,1925.306].
    ## Model rank =  425 / 425 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                                   k'      edf k-index p-value    
    ## s(year)                     9.00e+00 1.47e+00    0.95    0.16    
    ## s(year):River_sectionTI     8.00e+00 5.21e+00    0.95    0.11    
    ## s(year):River_sectionBR     8.00e+00 4.22e+00    0.95    0.14    
    ## s(year):River_sectionLSL    8.00e+00 5.60e-04    0.95    0.17    
    ## s(year):River_sectionCA     8.00e+00 4.79e-03    0.95    0.15    
    ## s(year):River_sectionLSF    8.00e+00 1.00e-02    0.95    0.14    
    ## s(River_section)            5.00e+00 4.12e-04      NA      NA    
    ## s(DOY)                      9.00e+00 1.00e+00    0.95    0.18    
    ## s(NS_MC)                    2.00e+00 9.03e-01      NA      NA    
    ## s(Longitude,Latitude)       1.90e+01 1.54e+01    0.97    0.56    
    ## s(Data_source)              1.70e+01 1.08e+01      NA      NA    
    ## ti(Longitude,Latitude,year) 1.53e+02 4.38e+01    0.86  <2e-16 ***
    ## ti(Longitude,Latitude,DOY)  1.70e+02 4.18e+01    0.81  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TPM3 <- gam(TP ~ s(year, m = 2) + s(year, by = River_section, m = 1) + s(River_section,
    bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") + s(Longitude, Latitude,
    bs = "ds", m = c(1, 0.5), k = 20) + s(Data_source, bs = "re") + ti(Longitude,
    Latitude, year, d = c(2, 1), bs = c("ds", "tp"), m = list(c(1, 0.5)), k = 20) +
    ti(Longitude, Latitude, DOY, d = c(2, 1), bs = c("ds", "tp"), m = list(c(1, 0.5)),
        k = 20) + ti(year, DOY, bs = "tp", k = 20), data = spring, method = "REML",
    family = tw(link = "log"))
gam.check(TPM3)
```

![](Nutrients_files/figure-gfm/TP%20GAM%20building-4.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 11 iterations.
    ## Gradient range [-0.002951588,0.006072483]
    ## (score -9032.355 & scale 0.1387946).
    ## Hessian positive definite, eigenvalue range [0.0002845247,1924.793].
    ## Model rank =  746 / 746 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                                   k'      edf k-index p-value    
    ## s(year)                     9.00e+00 1.03e+00    0.97    0.47    
    ## s(year):River_sectionTI     8.00e+00 5.74e+00    0.97    0.50    
    ## s(year):River_sectionBR     8.00e+00 5.92e+00    0.97    0.52    
    ## s(year):River_sectionLSL    8.00e+00 1.11e-03    0.97    0.49    
    ## s(year):River_sectionCA     8.00e+00 1.13e-02    0.97    0.45    
    ## s(year):River_sectionLSF    8.00e+00 1.35e-03    0.97    0.48    
    ## s(River_section)            5.00e+00 1.75e-03      NA      NA    
    ## s(DOY)                      9.00e+00 1.00e+00    0.98    0.70    
    ## s(NS_MC)                    2.00e+00 8.89e-01      NA      NA    
    ## s(Longitude,Latitude)       1.90e+01 1.58e+01    0.98    0.68    
    ## s(Data_source)              1.70e+01 1.06e+01      NA      NA    
    ## ti(Longitude,Latitude,year) 1.29e+02 3.51e+01    0.86  <2e-16 ***
    ## ti(Longitude,Latitude,DOY)  1.54e+02 2.67e+01    0.83  <2e-16 ***
    ## ti(year,DOY)                3.61e+02 1.03e+02    0.89  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TPM4 <- gam(TP ~ s(year, m = 2) + s(year, by = River_section, m = 1) + s(River_section,
    bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") + s(Longitude, Latitude,
    bs = "ds", m = c(1, 0.5), k = 20) + s(Data_source, bs = "re") + ti(Longitude,
    Latitude, year, d = c(2, 1), bs = c("ds", "tp"), m = list(c(1, 0.5)), k = 20) +
    ti(Longitude, Latitude, DOY, d = c(2, 1), bs = c("ds", "tp"), m = list(c(1, 0.5)),
        k = 20) + ti(year, DOY, bs = "tp", k = 20) + s(Total_Precip, bs = "tp", k = 10),
    data = spring, method = "REML", family = tw(link = "log"))
gam.check(TPM4)
```

![](Nutrients_files/figure-gfm/TP%20GAM%20building-5.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 13 iterations.
    ## Gradient range [-0.002233202,0.008386099]
    ## (score -8887.602 & scale 0.1372224).
    ## Hessian positive definite, eigenvalue range [0.000169028,1904.685].
    ## Model rank =  747 / 747 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                                   k'      edf k-index p-value    
    ## s(year)                     9.00e+00 1.06e+00    0.98    0.70    
    ## s(year):River_sectionTI     8.00e+00 5.60e+00    0.98    0.66    
    ## s(year):River_sectionBR     8.00e+00 6.04e+00    0.98    0.67    
    ## s(year):River_sectionLSL    8.00e+00 4.21e-04    0.98    0.72    
    ## s(year):River_sectionCA     8.00e+00 1.65e-03    0.98    0.67    
    ## s(year):River_sectionLSF    8.00e+00 8.31e-03    0.98    0.67    
    ## s(River_section)            5.00e+00 6.21e-04      NA      NA    
    ## s(DOY)                      9.00e+00 1.00e+00    0.98    0.71    
    ## s(NS_MC)                    2.00e+00 8.84e-01      NA      NA    
    ## s(Longitude,Latitude)       1.90e+01 1.57e+01    0.95    0.22    
    ## s(Data_source)              1.70e+01 1.09e+01      NA      NA    
    ## ti(Longitude,Latitude,year) 1.26e+02 3.40e+01    0.86  <2e-16 ***
    ## ti(Longitude,Latitude,DOY)  1.49e+02 2.58e+01    0.83  <2e-16 ***
    ## ti(year,DOY)                3.61e+02 9.91e+01    0.89  <2e-16 ***
    ## s(Total_Precip)             9.00e+00 6.44e+00    0.96    0.26    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# compare models
AIC(TPBM, TPM1, TPM2, TPM3, TPM4)
```

    ##             df       AIC
    ## TPBM  61.69908 -18036.77
    ## TPM1 104.18801 -18112.02
    ## TPM2 144.09987 -18158.47
    ## TPM3 227.20839 -18315.60
    ## TPM4 229.39138 -18031.27

``` r
summary(TPM3)
```

    ## 
    ## Family: Tweedie(p=1.829) 
    ## Link function: log 
    ## 
    ## Formula:
    ## TP ~ s(year, m = 2) + s(year, by = River_section, m = 1) + s(River_section, 
    ##     bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") + s(Longitude, 
    ##     Latitude, bs = "ds", m = c(1, 0.5), k = 20) + s(Data_source, 
    ##     bs = "re") + ti(Longitude, Latitude, year, d = c(2, 1), bs = c("ds", 
    ##     "tp"), m = list(c(1, 0.5)), k = 20) + ti(Longitude, Latitude, 
    ##     DOY, d = c(2, 1), bs = c("ds", "tp"), m = list(c(1, 0.5)), 
    ##     k = 20) + ti(year, DOY, bs = "tp", k = 20)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -4.1257     0.2687  -15.36   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf  Ref.df      F  p-value    
    ## s(year)                     1.034e+00   1.040 17.511  2.5e-05 ***
    ## s(year):River_sectionTI     5.736e+00   8.000  3.389  < 2e-16 ***
    ## s(year):River_sectionBR     5.924e+00   8.000  2.526  < 2e-16 ***
    ## s(year):River_sectionLSL    1.108e-03   3.000  0.000 0.667464    
    ## s(year):River_sectionCA     1.133e-02   8.000  0.001 0.299725    
    ## s(year):River_sectionLSF    1.354e-03   8.000  0.000 0.227168    
    ## s(River_section)            1.752e-03   4.000  0.000 0.475800    
    ## s(DOY)                      1.001e+00   1.002  0.006 0.943755    
    ## s(NS_MC)                    8.886e-01   1.000 15.999 0.000802 ***
    ## s(Longitude,Latitude)       1.578e+01  19.000 69.914 0.009873 ** 
    ## s(Data_source)              1.062e+01  16.000  2.982  < 2e-16 ***
    ## ti(Longitude,Latitude,year) 3.507e+01 129.000  0.850  < 2e-16 ***
    ## ti(Longitude,Latitude,DOY)  2.672e+01 154.000  0.413  < 2e-16 ***
    ## ti(year,DOY)                1.026e+02 141.738  1.799  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.349   Deviance explained = 61.6%
    ## -REML = -9032.4  Scale est. = 0.13879   n = 2518

### TP GAM model checking

``` r
# Basic diagnostics:
k.check(TPM3)
```

    ##                              k'          edf   k-index p-value
    ## s(year)                       9 1.034399e+00 0.9723072  0.5025
    ## s(year):River_sectionTI       8 5.735998e+00 0.9723072  0.4300
    ## s(year):River_sectionBR       8 5.923688e+00 0.9723072  0.4650
    ## s(year):River_sectionLSL      8 1.108338e-03 0.9723072  0.4700
    ## s(year):River_sectionCA       8 1.133145e-02 0.9723072  0.4600
    ## s(year):River_sectionLSF      8 1.353598e-03 0.9723072  0.4700
    ## s(River_section)              5 1.752391e-03        NA      NA
    ## s(DOY)                        9 1.001262e+00 0.9839665  0.6825
    ## s(NS_MC)                      2 8.885640e-01        NA      NA
    ## s(Longitude,Latitude)        19 1.577736e+01 0.9839904  0.6050
    ## s(Data_source)               17 1.061892e+01        NA      NA
    ## ti(Longitude,Latitude,year) 129 3.507427e+01 0.8599132  0.0000
    ## ti(Longitude,Latitude,DOY)  154 2.671825e+01 0.8346472  0.0000
    ## ti(year,DOY)                361 1.026303e+02 0.8877018  0.0000

``` r
appraise(TPM3)
```

![](Nutrients_files/figure-gfm/TP%20GAM%20checking-1.png)<!-- -->

``` r
# Autocorrelation check:
par(mar = c(3, 3, 3, 0), mfrow = c(1, 2))
acf(resid(TPM3), lag.max = 36, main = "ACF")
pacf(resid(TPM3), lag.max = 36, main = "pACF")
```

![](Nutrients_files/figure-gfm/TP%20GAM%20checking-2.png)<!-- -->

``` r
# concurvity check:
concurvity(TPM3, full = TRUE)  #if any of the values at worst are higher than 0.8
```

    ##          para s(year) s(year):River_sectionTI s(year):River_sectionBR
    ## worst       1       1                       1               1.0000000
    ## observed    1       1                       1               0.9999999
    ## estimate    1       1                       1               1.0000000
    ##          s(year):River_sectionLSL s(year):River_sectionCA
    ## worst                    1.160654               1.0000000
    ## observed                 1.000000               1.0000000
    ## estimate                 1.000000               0.9999964
    ##          s(year):River_sectionLSF s(River_section)    s(DOY)  s(NS_MC)
    ## worst                           1                1 0.9999997 1.0000000
    ## observed                        1                1 0.9997352 0.9981252
    ## estimate                        1                1 0.9997247 0.9990626
    ##          s(Longitude,Latitude) s(Data_source) ti(Longitude,Latitude,year)
    ## worst                1.0000000      1.0000000                   1.0000000
    ## observed             0.9999997      0.9915513                   0.7700428
    ## estimate             1.0000000      0.9937514                   0.5693328
    ##          ti(Longitude,Latitude,DOY) ti(year,DOY)
    ## worst                     1.0000000    1.0000000
    ## observed                  0.6667135    0.2665683
    ## estimate                  0.5970696    0.1477582

### Plot TP GAM using gratia

``` r
# First look at year partial effects for each river section
draw(TPM3, select = c(1:6), residuals = FALSE) & theme_bw()
```

![](Nutrients_files/figure-gfm/Gratia%20TP%20GAM%20plots-1.png)<!-- -->

``` r
# Year partial effects with long + lat
draw(TPM3, select = 10, residuals = FALSE, continuous_fill = ggplot2::scale_fill_distiller(palette = "PuOr",
    type = "div"))
```

![](Nutrients_files/figure-gfm/Gratia%20TP%20GAM%20plots-2.png)<!-- -->

``` r
draw(TPM3, select = 12, residuals = TRUE, continuous_fill = ggplot2::scale_fill_distiller(palette = "PuOr",
    type = "div", na.value = "transparent"))
```

![](Nutrients_files/figure-gfm/Gratia%20TP%20GAM%20plots-3.png)<!-- -->

### Let’s customize the plots

``` r
# customized partial effect plots
m2.d <- derivatives(TPM3, term = "s(year)", n = 200)
m2.d <- m2.d %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))

sm_year <- smooth_estimates(TPM3, smooth = "s(year)", dist = 0.05)
sm_year$upper <- sm_year$est + 1.96 * sm_year$se
sm_year$lower <- sm_year$est - 1.96 * sm_year$se

TP_year_plot <- ggplot(sm_year, aes(x = year, y = est)) + geom_ribbon(alpha = 0.5,
    aes(ymin = lower, ymax = upper, fill = "confidence interval")) + geom_line(data = sm_year,
    aes(x = year, y = est, color = "GAM"), linewidth = 1) + scale_fill_manual(values = "lightblue",
    name = NULL) + scale_color_manual(values = "blue", name = NULL) + ylab("Partial effect of year on TP") +
    scale_x_continuous(name = "Year", breaks = c(1970, 1980, 1990, 2000, 2010, 2020)) +
    theme_minimal(base_size = 16) + theme(legend.position = "none")

TP_year_plot
```

![](Nutrients_files/figure-gfm/TP%20GAM%20plots-1.png)<!-- -->

``` r
ggsave("TP_year_plot.svg", TP_year_plot, width = 12, height = 8, dpi = 600)

sm_year_TI <- smooth_estimates(TPM3, smooth = "s(year):River_sectionTI")
sm_year_TI$upper <- sm_year_TI$est + 1.96 * sm_year_TI$se
sm_year_TI$lower <- sm_year_TI$est - 1.96 * sm_year_TI$se

sm_year_TI <- sm_year_TI %>%
    mutate(sig = ifelse(upper < 0 | lower > 0, "sig", "ns"))

d_tp_TI <- derivatives(TPM3, term = "s(year):River_sectionTI", n = 200)
d_tp_TI <- d_tp_TI %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))

p_year_TI <- ggplot(sm_year_TI, aes(x = year, y = est)) + annotate("rect", xmin = c(1975.201,
    1992.035), xmax = c(1980.925, 1995.739), ymin = -Inf, ymax = Inf, fill = "purple",
    alpha = 0.2) + annotate("rect", xmin = c(1966.111, 1984.291), xmax = c(1966.784,
    1987.995), ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.2) + geom_ribbon(aes(ymin = lower,
    ymax = upper, fill = "confidence interval"), alpha = 0.5) + geom_line(data = sm_year_TI,
    aes(x = year, y = est, color = "GAM"), linewidth = 1) + scale_fill_manual(values = "lightblue",
    name = NULL) + scale_color_manual(values = "darkblue", name = NULL) + ylab("Partial effect \n of year on TP") +
    scale_x_continuous(name = "Year", breaks = c(1970, 1980, 1990, 2000, 2010, 2020)) +
    theme_minimal(base_size = 16) + theme(legend.position = "none", axis.title.x = element_blank())

sm_year_BR <- smooth_estimates(TPM3, smooth = "s(year):River_sectionBR")
sm_year_BR$upper <- sm_year_BR$est + 1.96 * sm_year_BR$se
sm_year_BR$lower <- sm_year_BR$est - 1.96 * sm_year_BR$se

sm_year_BR <- sm_year_BR %>%
    mutate(sig = ifelse(upper < 0 | lower > 0, "sig", "ns"))

d_tp_BR <- derivatives(TPM3, term = "s(year):River_sectionBR", n = 200)
d_tp_BR <- d_tp_BR %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))

p_year_BR <- ggplot(sm_year_BR, aes(x = year, y = est)) + annotate("rect", xmin = c(1976.211,
    1992.035), xmax = c(1980.925, 1992.709), ymin = -Inf, ymax = Inf, fill = "purple",
    alpha = 0.2) + annotate("rect", xmin = c(1966.111, 1999.106, 2014.256), xmax = c(1970.824,
    2001.799, 2017.96), ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "confidence interval"), alpha = 0.5) +
    geom_line(data = sm_year_BR, aes(x = year, y = est, color = "GAM"), linewidth = 1) +
    scale_fill_manual(values = "lightblue", name = NULL) + scale_color_manual(values = "darkblue",
    name = NULL) + ylab("Partial effect \n of year on TP") + scale_x_continuous(name = "Year",
    breaks = c(1970, 1980, 1990, 2000, 2010, 2020)) + theme_minimal(base_size = 16) +
    theme(legend.position = "none")

ggarrange(p_year_TI, p_year_BR, labels = c("Thousand Islands", "Brockville Narrows"),
    nrow = 2, vjust = 0.2) + theme(plot.margin = margin(10, 1, 1, 1))
```

![](Nutrients_files/figure-gfm/TP%20GAM%20plots-2.png)<!-- -->

``` r
sm <- smooth_estimates(TPM3, smooth = "s(Longitude,Latitude)", dist = 0.05)
# ggplot(sm, aes(x = Longitude, y = Latitude)) + geom_raster(aes(fill = est)) +
# geom_point(data = spring, alpha = 0.2) + # add a point layer for original
# data scale_fill_distiller(palette='PuOr',type='div',na.value='white')+
# theme_minimal()

library(raster)
sm1 <- as.data.frame(sm[, c(6:7, 4)])
dfr <- rasterFromXYZ(sm1, crs = "+datum=WGS84")
plot(dfr)
```

![](Nutrients_files/figure-gfm/TP%20GAM%20plots-3.png)<!-- -->

``` r
writeRaster(dfr, "TP_GAM.tif", overwrite = TRUE)
# We will visualize the background in GIS software
```

## TN GAMs

``` r
# Try base model
TNM1 <- gam(TN ~ s(year, k = 20, m = 2, bs = "tp") + s(year, by = River_section,
    m = 1) + s(River_section, bs = "re") + s(Data_source, bs = "re") + s(DOY, bs = "tp") +
    s(Longitude, Latitude, bs = "ds", m = c(1, 0.5), k = 10) + s(NS_MC, bs = "re"),
    data = spring, method = "REML", family = tw(link = "log"))
k.check(TNM1)
```

    ##                          k'          edf   k-index p-value
    ## s(year)                  19 4.279329e+00 0.9355260  0.0025
    ## s(year):River_sectionTI   8 2.214023e+00 0.9355260  0.0150
    ## s(year):River_sectionBR   8 2.801673e-03 0.9355260  0.0025
    ## s(year):River_sectionCA   8 5.197892e+00 0.9355260  0.0000
    ## s(year):River_sectionLSF  8 3.088715e+00 0.9355260  0.0125
    ## s(River_section)          4 4.929627e-07        NA      NA
    ## s(Data_source)           11 7.098161e+00        NA      NA
    ## s(DOY)                    9 3.054379e+00 0.9967393  0.5625
    ## s(Longitude,Latitude)     9 7.107281e+00 0.9831844  0.3775
    ## s(NS_MC)                  2 9.675193e-01        NA      NA

``` r
appraise(TNM1)
```

![](Nutrients_files/figure-gfm/TN%20GAM%20model%20building-1.png)<!-- -->

``` r
par(mar = c(3, 3, 3, 0), mfrow = c(2, 2))
acf(resid(TNM1), lag.max = 36, main = "ACF")
pacf(resid(TNM1), lag.max = 36, main = "pACF")

# check if need to account for autocorrelation or moving average
gam_errors <- residuals(TNM1, type = "response")
error_mod <- auto.arima(gam_errors)
error_mod
```

    ## Series: gam_errors 
    ## ARIMA(1,0,3) with non-zero mean 
    ## 
    ## Coefficients:
    ##          ar1      ma1     ma2      ma3   mean
    ##       0.8931  -0.8596  0.1807  -0.2262  8e-04
    ## s.e.  0.0674   0.0684  0.0310   0.0230  5e-03
    ## 
    ## sigma^2 = 0.05735:  log likelihood = 21.47
    ## AIC=-30.95   AICc=-30.9   BIC=2.14

``` r
# check recommended order for AR (p) and/or MA (q)
checkresiduals(error_mod, theme = theme_bw())
```

    ## 
    ##  Ljung-Box test
    ## 
    ## data:  Residuals from ARIMA(1,0,3) with non-zero mean
    ## Q* = 23.879, df = 6, p-value = 0.0005496
    ## 
    ## Model df: 4.   Total lags used: 10

``` r
# suggests p=4 for summer data

# add autoregressive moving average (ARMA)
TNM1_ARMA0 <- gamm(TN ~ s(year, k = 10, bs = "tp") + s(year, River_section, bs = "fs") +
    s(River_section, bs = "re") + s(Data_source, bs = "re") + s(DOY, bs = "tp") +
    s(Longitude, Latitude, bs = "ds", m = c(1, 0.5), k = 10) + s(NS_MC, bs = "re"),
    data = spring, family = "tw")
```

    ## 
    ##  Maximum number of PQL iterations:  20

``` r
auto.arima(resid(TNM1_ARMA0$lme, type = "normalized"), stationary = TRUE, seasonal = FALSE)
```

    ## Series: resid(TNM1_ARMA0$lme, type = "normalized") 
    ## ARIMA(2,0,3) with zero mean 
    ## 
    ## Coefficients:
    ##           ar1     ar2     ma1      ma2      ma3
    ##       -0.0071  0.8186  0.0541  -0.7509  -0.0382
    ## s.e.   0.0785  0.0679  0.0825   0.0788   0.0259
    ## 
    ## sigma^2 = 0.9725:  log likelihood = -2574.26
    ## AIC=5160.53   AICc=5160.57   BIC=5193.61

``` r
# indicates summer model should be p=1, q=3
TNM1_AR1 <- gamm(TN ~ s(year, k = 10, bs = "tp") + s(year, River_section, bs = "fs") +
    s(River_section, bs = "re") + s(Data_source, bs = "re") + s(DOY, bs = "tp") +
    s(Longitude, Latitude, bs = "ds", m = c(1, 0.5), k = 10) + s(NS_MC, bs = "re"),
    data = spring, correlation = corARMA(form = ~1 | year, p = 1, q = 3), family = "tw")
```

    ## 
    ##  Maximum number of PQL iterations:  20

``` r
appraise(TNM1_AR1$gam)
layout(matrix(1:2, ncol = 2))
res <- resid(TNM1_AR1$lme, type = "normalized")  #need to do on normalized resid
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
```

![](Nutrients_files/figure-gfm/TN%20GAM%20model%20building-2.png)<!-- -->

``` r
layout(1)

anova(TNM1_ARMA0$lme, TNM1_AR1$lme)
```

    ##                Model df      AIC      BIC    logLik   Test  L.Ratio p-value
    ## TNM1_ARMA0$lme     1 13 151.5756 223.2609 -62.78781                        
    ## TNM1_AR1$lme       2 17 134.3747 228.1170 -50.18733 1 vs 2 25.20097  <.0001

``` r
# factor-smooth interaction of year with river section
TNBM <- gam(TN ~ s(year, k = 10, m = 2, bs = "tp") + s(year, by = River_section,
    m = 1) + s(River_section, bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") +
    s(Longitude, Latitude, bs = "ds", m = c(1, 0.5), k = 10) + s(Data_source, bs = "re"),
    data = spring, method = "REML", family = tw(link = "log"))
gam.check(TNBM)
```

![](Nutrients_files/figure-gfm/TN%20GAM%20building-1.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 17 iterations.
    ## Gradient range [-0.001335813,0.006343199]
    ## (score -1349.07 & scale 0.03819407).
    ## Hessian positive definite, eigenvalue range [2.313155e-07,927.3387].
    ## Model rank =  77 / 77 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                                k'      edf k-index p-value    
    ## s(year)                  9.00e+00 4.08e+00    0.94  <2e-16 ***
    ## s(year):River_sectionTI  8.00e+00 2.23e+00    0.94  <2e-16 ***
    ## s(year):River_sectionBR  8.00e+00 1.50e-02    0.94  <2e-16 ***
    ## s(year):River_sectionCA  8.00e+00 5.20e+00    0.94   0.010 ** 
    ## s(year):River_sectionLSF 8.00e+00 3.09e+00    0.94   0.015 *  
    ## s(River_section)         4.00e+00 5.70e-07      NA      NA    
    ## s(DOY)                   9.00e+00 3.05e+00    1.00   0.500    
    ## s(NS_MC)                 2.00e+00 9.67e-01      NA      NA    
    ## s(Longitude,Latitude)    9.00e+00 7.11e+00    0.98   0.425    
    ## s(Data_source)           1.10e+01 7.10e+00      NA      NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TNM1 <- gam(TN ~ s(year, k = 10, m = 2, bs = "tp") + s(year, by = River_section,
    m = 1) + s(River_section, bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") +
    s(Longitude, Latitude, bs = "ds", m = c(1, 0.5), k = 10) + s(Data_source, bs = "re") +
    s(Total_Precip), data = spring, method = "REML", family = tw(link = "log"))

gam.check(TNM1)
```

![](Nutrients_files/figure-gfm/TN%20GAM%20building-2.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 14 iterations.
    ## Gradient range [-0.0003780494,0.0008245085]
    ## (score -1351.444 & scale 0.03764099).
    ## Hessian positive definite, eigenvalue range [8.175968e-06,923.1323].
    ## Model rank =  86 / 86 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                                k'      edf k-index p-value   
    ## s(year)                  9.00e+00 4.16e+00    0.94   0.015 * 
    ## s(year):River_sectionTI  8.00e+00 1.26e+00    0.94   0.010 **
    ## s(year):River_sectionBR  8.00e+00 1.79e+00    0.94   0.025 * 
    ## s(year):River_sectionCA  8.00e+00 5.29e+00    0.94   0.005 **
    ## s(year):River_sectionLSF 8.00e+00 3.12e+00    0.94   0.015 * 
    ## s(River_section)         4.00e+00 2.04e-05      NA      NA   
    ## s(DOY)                   9.00e+00 3.78e+00    0.99   0.405   
    ## s(NS_MC)                 2.00e+00 9.68e-01      NA      NA   
    ## s(Longitude,Latitude)    9.00e+00 7.12e+00    0.95   0.100 . 
    ## s(Data_source)           1.10e+01 7.07e+00      NA      NA   
    ## s(Total_Precip)          9.00e+00 3.48e+00    0.94   0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TNM2 <- gam(TN ~ s(year, k = 10, m = 2, bs = "tp") + s(year, by = River_section,
    m = 1) + s(River_section, bs = "re") + s(DOY, bs = "tp") + s(NS_MC, bs = "re") +
    s(Longitude, Latitude, bs = "ds", m = c(1, 0.5), k = 10) + s(Data_source, bs = "re") +
    s(Total_Precip) + ti(year, DOY, bs = "tp"), data = spring, method = "REML", family = tw(link = "log"))
gam.check(TNM2)
```

![](Nutrients_files/figure-gfm/TN%20GAM%20building-3.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 22 iterations.
    ## Gradient range [-8.679871e-05,0.000391237]
    ## (score -1357.467 & scale 0.03733658).
    ## Hessian positive definite, eigenvalue range [5.223632e-08,922.5388].
    ## Model rank =  102 / 102 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                                k'      edf k-index p-value   
    ## s(year)                  9.00e+00 4.14e+00    0.95   0.015 * 
    ## s(year):River_sectionTI  8.00e+00 2.51e+00    0.95   0.010 **
    ## s(year):River_sectionBR  8.00e+00 3.15e-04    0.95   0.020 * 
    ## s(year):River_sectionCA  8.00e+00 5.30e+00    0.95   0.055 . 
    ## s(year):River_sectionLSF 8.00e+00 3.13e+00    0.95   0.010 **
    ## s(River_section)         4.00e+00 4.76e-07      NA      NA   
    ## s(DOY)                   9.00e+00 3.31e+00    0.99   0.425   
    ## s(NS_MC)                 2.00e+00 9.66e-01      NA      NA   
    ## s(Longitude,Latitude)    9.00e+00 7.10e+00    0.95   0.175   
    ## s(Data_source)           1.10e+01 7.09e+00      NA      NA   
    ## s(Total_Precip)          9.00e+00 3.23e+00    0.95   0.030 * 
    ## ti(year,DOY)             1.60e+01 2.53e+00    0.94   0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# compare models
AIC(TNBM, TNM1, TNM2)
```

    ##            df       AIC
    ## TNBM 41.15897 -2783.441
    ## TNM1 47.64612 -2797.819
    ## TNM2 49.55178 -2810.306

``` r
summary(TNM2)
```

    ## 
    ## Family: Tweedie(p=1.99) 
    ## Link function: log 
    ## 
    ## Formula:
    ## TN ~ s(year, k = 10, m = 2, bs = "tp") + s(year, by = River_section, 
    ##     m = 1) + s(River_section, bs = "re") + s(DOY, bs = "tp") + 
    ##     s(NS_MC, bs = "re") + s(Longitude, Latitude, bs = "ds", m = c(1, 
    ##     0.5), k = 10) + s(Data_source, bs = "re") + s(Total_Precip) + 
    ##     ti(year, DOY, bs = "tp")
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  -0.2811     0.2576  -1.092    0.275
    ## 
    ## Approximate significance of smooth terms:
    ##                                edf Ref.df        F  p-value    
    ## s(year)                  4.135e+00  4.966   11.950  < 2e-16 ***
    ## s(year):River_sectionTI  2.512e+00  8.000    0.920 0.010729 *  
    ## s(year):River_sectionBR  3.151e-04  8.000    0.000 0.086353 .  
    ## s(year):River_sectionCA  5.295e+00  8.000    3.336 0.000588 ***
    ## s(year):River_sectionLSF 3.125e+00  8.000    3.025 5.52e-06 ***
    ## s(River_section)         4.765e-07  3.000    0.000 0.419914    
    ## s(DOY)                   3.314e+00  4.105    1.840 0.104868    
    ## s(NS_MC)                 9.657e-01  1.000  117.330  < 2e-16 ***
    ## s(Longitude,Latitude)    7.105e+00  9.000 1066.543 0.008159 ** 
    ## s(Data_source)           7.092e+00 10.000   26.288  < 2e-16 ***
    ## s(Total_Precip)          3.226e+00  4.044    5.579 0.000175 ***
    ## ti(year,DOY)             2.531e+00  3.538    5.193 0.000721 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.665   Deviance explained = 77.1%
    ## -REML = -1357.5  Scale est. = 0.037337  n = 1827

### TN GAM model checking

``` r
# Basic diagnostics:
k.check(TNM2)
```

    ##                          k'          edf   k-index p-value
    ## s(year)                   9 4.135307e+00 0.9499143  0.0125
    ## s(year):River_sectionTI   8 2.512062e+00 0.9499143  0.0175
    ## s(year):River_sectionBR   8 3.150509e-04 0.9499143  0.0275
    ## s(year):River_sectionCA   8 5.295134e+00 0.9499143  0.0250
    ## s(year):River_sectionLSF  8 3.125208e+00 0.9499143  0.0300
    ## s(River_section)          4 4.764818e-07        NA      NA
    ## s(DOY)                    9 3.314247e+00 0.9890429  0.4100
    ## s(NS_MC)                  2 9.656716e-01        NA      NA
    ## s(Longitude,Latitude)     9 7.104604e+00 0.9515165  0.1800
    ## s(Data_source)           11 7.091529e+00        NA      NA
    ## s(Total_Precip)           9 3.226361e+00 0.9521725  0.0350
    ## ti(year,DOY)             16 2.530511e+00 0.9392198  0.0175

``` r
appraise(TNM2)
```

![](Nutrients_files/figure-gfm/TN%20GAM%20checking-1.png)<!-- -->

``` r
# Autocorrelation check:
par(mar = c(3, 3, 3, 0), mfrow = c(1, 2))
acf(resid(TNM2), lag.max = 36, main = "ACF")
pacf(resid(TNM2), lag.max = 36, main = "pACF")
```

![](Nutrients_files/figure-gfm/TN%20GAM%20checking-2.png)<!-- -->

``` r
# concurvity check:
concurvity(TNM2, full = TRUE)  #if any of the values at worst are higher than 0.8, you want to look at full=FALSE as well.
```

    ##          para   s(year) s(year):River_sectionTI s(year):River_sectionBR
    ## worst       1 1.0000000               1.0000000               1.0000000
    ## observed    1 0.9998646               0.9999843               0.9993100
    ## estimate    1 0.9994937               0.9973289               0.9951571
    ##          s(year):River_sectionCA s(year):River_sectionLSF s(River_section)
    ## worst                  1.0000000                1.0000000                1
    ## observed               0.9972614                0.9996465                1
    ## estimate               0.9905921                0.9975868                1
    ##              s(DOY) s(NS_MC) s(Longitude,Latitude) s(Data_source)
    ## worst    0.13361967        1             1.0000000      1.0000000
    ## observed 0.10937345        1             0.7227858      0.9998780
    ## estimate 0.09712706        1             0.9990019      0.9915557
    ##          s(Total_Precip) ti(year,DOY)
    ## worst          0.3239505   0.18447255
    ## observed       0.2088697   0.14040845
    ## estimate       0.2303612   0.08576915

### Plot TN GAM

``` r
# Year partial effects at each river section
draw(TNM2, select = c(1:5), residuals = FALSE) & theme_bw()
```

![](Nutrients_files/figure-gfm/TN%20GAM%20plots-1.png)<!-- -->

``` r
# Year partial effects with long + lat
draw(TNM2, select = 9, residuals = TRUE, continuous_fill = ggplot2::scale_fill_distiller(palette = "PuOr",
    type = "div")) & theme_bw()
```

![](Nutrients_files/figure-gfm/TN%20GAM%20plots-2.png)<!-- -->

``` r
# Precipitation partial effects
draw(TNM2, select = 11, residuals = FALSE) & theme_bw()
```

![](Nutrients_files/figure-gfm/TN%20GAM%20plots-3.png)<!-- -->

``` r
# Now let's customize the plots
tn_sm_year <- smooth_estimates(TNM2, smooth = "s(year)", dist = 0.05)
tn_sm_year$upper <- tn_sm_year$est + 1.96 * tn_sm_year$se
tn_sm_year$lower <- tn_sm_year$est - 1.96 * tn_sm_year$se

tn_sm_year <- tn_sm_year %>%
    mutate(sig = ifelse(upper < 0 | lower > 0, "sig", "ns"))

d_tn <- derivatives(TNM2, term = "s(year)", n = 200)
d_tn <- d_tn %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))

tn_gam_plot <- ggplot(tn_sm_year, aes(x = year, y = est)) + annotate("rect", xmin = c(2006.176),
    xmax = c(2022), ymin = -Inf, ymax = Inf, fill = "purple", alpha = 0.2) + annotate("rect",
    xmin = c(1990.015), xmax = c(1995.065), ymin = -Inf, ymax = Inf, fill = "orange",
    alpha = 0.2) + geom_ribbon(alpha = 0.5, aes(ymin = lower, ymax = upper, fill = "confidence interval")) +
    geom_line(data = tn_sm_year, aes(x = year, y = est, color = "GAM"), linewidth = 1) +
    scale_fill_manual(values = "lightblue", name = NULL) + scale_color_manual(values = "darkblue",
    name = NULL) + ylab("Partial effect of year on TN") + scale_x_continuous(name = "Year",
    breaks = c(1970, 1980, 1990, 2000, 2010, 2020)) + theme_minimal(base_size = 16) +
    theme(legend.position = "none")

tn_gam_plot
```

![](Nutrients_files/figure-gfm/TN%20GAM%20plots-4.png)<!-- -->

``` r
ggsave("TN_year_plot.svg", tn_gam_plot, width = 12, height = 8, dpi = 600)


tn_sm_year_TI <- smooth_estimates(TNM2, smooth = "s(year):River_sectionTI", dist = 0.05)
tn_sm_year_TI$upper <- tn_sm_year_TI$est + 1.96 * tn_sm_year_TI$se
tn_sm_year_TI$lower <- tn_sm_year_TI$est - 1.96 * tn_sm_year_TI$se

d_tn_TI <- derivatives(TNM2, term = "s(year):River_sectionTI", n = 200)
d_tn_TI <- d_tn_TI %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))
# no significant changes

tn_sm_year_BR <- smooth_estimates(TNM2, smooth = "s(year):River_sectionBR", dist = 0.05)
tn_sm_year_BR$upper <- tn_sm_year_BR$est + 1.96 * tn_sm_year_BR$se
tn_sm_year_BR$lower <- tn_sm_year_BR$est - 1.96 * tn_sm_year_BR$se

d_tn_BR <- derivatives(TNM2, term = "s(year):River_sectionBR", n = 200)
d_tn_BR <- d_tn_BR %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))
# no significant changes

tn_sm_year_CA <- smooth_estimates(TNM2, smooth = "s(year):River_sectionCA", dist = 0.05)
tn_sm_year_CA$upper <- tn_sm_year_CA$est + 1.96 * tn_sm_year_CA$se
tn_sm_year_CA$lower <- tn_sm_year_CA$est - 1.96 * tn_sm_year_CA$se

d_tn_CA <- derivatives(TNM2, term = "s(year):River_sectionCA", n = 200)
d_tn_CA <- d_tn_CA %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))

t_year_CA <- ggplot(tn_sm_year_CA, aes(x = year, y = est)) + annotate("rect", xmin = c(1986.312),
    xmax = c(1991.698), ymin = -Inf, ymax = Inf, fill = "purple", alpha = 0.2) +
    annotate("rect", xmin = c(1976.211), xmax = c(1981.935), ymin = -Inf, ymax = Inf,
        fill = "orange", alpha = 0.2) + geom_ribbon(alpha = 0.5, aes(ymin = lower,
    ymax = upper, fill = "confidence interval")) + geom_line(data = tn_sm_year_CA,
    aes(x = year, y = est, color = "GAM"), linewidth = 1) + scale_fill_manual(values = "lightblue",
    name = NULL) + scale_color_manual(values = "darkblue", name = NULL) + ylab("Partial effect \n of year on TN") +
    scale_x_continuous(name = "Year", breaks = c(1970, 1980, 1990, 2000, 2010, 2020)) +
    theme_minimal(base_size = 16) + theme(legend.position = "none")

tn_sm_year_LSF <- smooth_estimates(TNM2, smooth = "s(year):River_sectionLSF", dist = 0.05)
tn_sm_year_LSF$upper <- tn_sm_year_LSF$est + 1.96 * tn_sm_year_LSF$se
tn_sm_year_LSF$lower <- tn_sm_year_LSF$est - 1.96 * tn_sm_year_LSF$se

d_tn_LSF <- derivatives(TNM2, term = "s(year):River_sectionLSF", n = 200)
d_tn_LSF <- d_tn_LSF %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))

t_year_LSF <- ggplot(tn_sm_year_LSF, aes(x = year, y = est)) + annotate("rect", xmin = c(2002.136),
    xmax = c(2003.819), ymin = -Inf, ymax = Inf, fill = "purple", alpha = 0.2) +
    annotate("rect", xmin = c(2008.196), xmax = c(2013.92), ymin = -Inf, ymax = Inf,
        fill = "orange", alpha = 0.2) + geom_ribbon(alpha = 0.5, aes(ymin = lower,
    ymax = upper, fill = "confidence interval")) + geom_line(data = tn_sm_year_LSF,
    aes(x = year, y = est, color = "GAM"), linewidth = 1) + scale_fill_manual(values = "lightblue",
    name = NULL) + scale_color_manual(values = "darkblue", name = NULL) + ylab("Partial effect \n of year on TN") +
    scale_x_continuous(name = "Year", breaks = c(1970, 1980, 1990, 2000, 2010, 2020)) +
    theme_minimal(base_size = 16) + theme(legend.position = "none")

tn_temp.gam_sections <- ggarrange(t_year_CA, t_year_LSF, labels = c("Cornwall", "Lake St. Francis"),
    nrow = 2, vjust = 0.2) + theme(plot.margin = margin(10, 1, 1, 1))
ggsave("TN_CA.LSF_temporalGAM.svg", tn_temp.gam_sections, dpi = 600, height = 10,
    width = 12)
tn_temp.gam_sections
```

![](Nutrients_files/figure-gfm/TN%20GAM%20plots-5.png)<!-- -->

``` r
sm_precip <- smooth_estimates(TNM2, smooth = "s(Total_Precip)", dist = 0.05)
sm_precip$upper <- sm_precip$est + 1.96 * sm_precip$se
sm_precip$lower <- sm_precip$est - 1.96 * sm_precip$se

sm_precip <- sm_precip %>%
    mutate(sig = ifelse(upper < 0 | lower > 0, "sig", "ns"))

d_tn_precip <- derivatives(TNM2, term = "s(Total_Precip)", n = 200)
d_tn_precip <- d_tn_precip %>%
    mutate(change = ifelse(upper < 0, "dec", if_else(lower > 0, "inc", "")))

TN_precip_plot <- ggplot(sm_precip, aes(x = Total_Precip, y = est)) + annotate("rect",
    xmin = c(10.4), xmax = c(66.48643), ymin = -Inf, ymax = Inf, fill = "orange",
    alpha = 0.2) + geom_ribbon(alpha = 0.5, aes(ymin = lower, ymax = upper, fill = "confidence interval")) +
    geom_line(data = sm_precip, aes(x = Total_Precip, y = est, color = "GAM"), linewidth = 1) +
    scale_fill_manual(values = "lightblue", name = NULL) + scale_color_manual(values = "darkblue",
    name = NULL) + ylab("Partial effect of precipitation on TN") + scale_x_continuous(name = "Precipitation") +
    theme_minimal(base_size = 16) + theme(legend.position = "none")
ggsave("TN_precip_plot.svg", TN_precip_plot, width = 12, height = 8, dpi = 600)
TN_precip_plot
```

![](Nutrients_files/figure-gfm/TN%20GAM%20plots-6.png)<!-- -->

## Water Quality Interactions

``` r
lm1 <- lm(log1p(CHLA) ~ log1p(TP), data = Ac_gov)
summary(lm1)
```

    ## 
    ## Call:
    ## lm(formula = log1p(CHLA) ~ log1p(TP), data = Ac_gov)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3729 -0.7668 -0.1695  0.2578  4.3987 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.74758    0.01517   49.27   <2e-16 ***
    ## log1p(TP)   19.23519    0.66700   28.84   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8685 on 4931 degrees of freedom
    ##   (6750 observations deleted due to missingness)
    ## Multiple R-squared:  0.1443, Adjusted R-squared:  0.1441 
    ## F-statistic: 831.6 on 1 and 4931 DF,  p-value: < 2.2e-16

``` r
p1 <- ggplot(data = Ac_gov, aes(log1p(TP), log1p(CHLA))) + geom_point()
p1 <- p1 + geom_smooth(method = lm, color = "blue", fill = "#69b3a2", se = TRUE) +
    theme(axis.title = element_text(size = 15), panel.background = element_rect(fill = "white")) +
    xlab("Log TP (mg*L-1)") + ylab("Log Chl a (ug*L-1)")
p1
```

![](Nutrients_files/figure-gfm/CHLA~TP%20relationship-1.png)<!-- -->

``` r
lm1 <- glm(TP ~ TSS, data = Nut_precip, family = Gamma(link = "log"))
summary(lm1)
```

    ## 
    ## Call:
    ## glm(formula = TP ~ TSS, family = Gamma(link = "log"), data = Nut_precip)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -4.8461  -0.9814  -0.5874   0.1482   4.8242  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -3.944217   0.045213  -87.24   <2e-16 ***
    ## TSS          0.033488   0.001893   17.69   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.918159)
    ## 
    ##     Null deviance: 1527.4  on 1070  degrees of freedom
    ## Residual deviance: 1216.9  on 1069  degrees of freedom
    ##   (10530 observations deleted due to missingness)
    ## AIC: -5684.3
    ## 
    ## Number of Fisher Scoring iterations: 25

``` r
par(mfrow = c(2, 2))
plot(lm1)
```

![](Nutrients_files/figure-gfm/TP~TSS%20relationship-1.png)<!-- -->

``` r
p2 <- ggplot(data = Ac_gov, aes(TSS, TP)) + geom_point()
p2 <- p2 + geom_smooth(method = lm, color = "blue", fill = "#69b3a2", se = TRUE) +
    theme(axis.title = element_text(size = 15), panel.background = element_rect(fill = "white")) +
    xlab("TSS (mg*L-1)") + ylab("TP (ug*L-1)") + xlim(0, 450) + stat_cor(aes(label = paste(..rr.label..,
    ..p.label.., sep = "~`,`~")), label.x = 200, label.y = 0.4)

p2
```

![](Nutrients_files/figure-gfm/TP~TSS%20relationship-2.png)<!-- -->

## Invasive Dreissenid impacts

``` r
# select stations with pre and post invasion data, limit to ice-free season
# (may-oct) within 10 years of invasion (1991)
icefree_dreis <- Ac_gov %>%
    filter(month > 4, month < 11 & Site == "20170010" | Site == "20180011") %>%
    group_by(Site, year, River_section, Longitude, Latitude, month) %>%
    summarise(CHLA = mean(CHLA, na.rm = TRUE), TDP = mean(TDP, na.rm = TRUE), SRP = mean(SRP,
        na.rm = TRUE), TP = mean(TP, na.rm = TRUE), NH3.4 = mean(NH3.4, na.rm = TRUE),
        NO2.3 = mean(NO2.3, na.rm = TRUE), TN = mean(TN, na.rm = TRUE))

icefree_dreis <- icefree_dreis %>%
    mutate(invade = if_else(year < 1991, "pre", "post"))
icefree_dreis$invade <- factor(icefree_dreis$invade, levels = c("pre", "post"))

icefree_dreis <- icefree_dreis %>%
    mutate(chltp = CHLA/(TP * 1000))

icefree_dreis$lchltp <- log1p(icefree_dreis$chltp)

king <- subset(icefree_dreis, River_section == "TI" & year < 2001 & year > 1980)
brock <- subset(icefree_dreis, River_section == "BR" & year < 2001 & year > 1980)

MW_king <- wilcox.test(lchltp ~ invade, data = king)
MW_king
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  lchltp by invade
    ## W = 2510, p-value = 0.0001958
    ## alternative hypothesis: true location shift is not equal to 0

``` r
MW_brock <- wilcox.test(lchltp ~ invade, data = brock)
MW_brock
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  lchltp by invade
    ## W = 5087.5, p-value = 7.682e-14
    ## alternative hypothesis: true location shift is not equal to 0

``` r
lm1 <- lm(log1p(CHLA) ~ log1p(TP), data = subset(king, invade == "pre"))
summary(lm1)
```

    ## 
    ## Call:
    ## lm(formula = log1p(CHLA) ~ log1p(TP), data = subset(king, invade == 
    ##     "pre"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.54397 -0.19016 -0.00888  0.17920  1.03552 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1.086      0.130   8.355 1.56e-11 ***
    ## log1p(TP)     12.144     10.873   1.117    0.269    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2956 on 58 degrees of freedom
    ## Multiple R-squared:  0.02105,    Adjusted R-squared:  0.004175 
    ## F-statistic: 1.247 on 1 and 58 DF,  p-value: 0.2687

``` r
par(mfrow = c(2, 2))
plot(lm1)
```

![](Nutrients_files/figure-gfm/pre/post%20dreissenid%20analysis-1.png)<!-- -->

``` r
lm2 <- lm(log1p(CHLA) ~ log1p(TP), data = subset(king, invade == "post"))
summary(lm2)
```

    ## 
    ## Call:
    ## lm(formula = log1p(CHLA) ~ log1p(TP), data = subset(king, invade == 
    ##     "post"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.58561 -0.15913 -0.01923  0.08607  0.75618 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.5357     0.1103   4.857 9.39e-06 ***
    ## log1p(TP)    27.4370    14.0134   1.958   0.0551 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2814 on 58 degrees of freedom
    ## Multiple R-squared:  0.062,  Adjusted R-squared:  0.04582 
    ## F-statistic: 3.833 on 1 and 58 DF,  p-value: 0.05506

``` r
plot(lm2)
```

![](Nutrients_files/figure-gfm/pre/post%20dreissenid%20analysis-2.png)<!-- -->

``` r
lm3 <- lm(log1p(CHLA) ~ log1p(TP), data = subset(brock, invade == "pre"))
summary(lm3)
```

    ## 
    ## Call:
    ## lm(formula = log1p(CHLA) ~ log1p(TP), data = subset(brock, invade == 
    ##     "pre"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.47184 -0.12955 -0.03439  0.10819  0.77654 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.0043     0.1264   7.946 1.78e-11 ***
    ## log1p(TP)    12.9389    11.0305   1.173    0.245    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2622 on 73 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.0185, Adjusted R-squared:  0.005055 
    ## F-statistic: 1.376 on 1 and 73 DF,  p-value: 0.2446

``` r
plot(lm3)
```

![](Nutrients_files/figure-gfm/pre/post%20dreissenid%20analysis-3.png)<!-- -->

``` r
lm4 <- lm(log1p(CHLA) ~ log1p(TP), data = subset(brock, invade == "post"))
summary(lm4)
```

    ## 
    ## Call:
    ## lm(formula = log1p(CHLA) ~ log1p(TP), data = subset(brock, invade == 
    ##     "post"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.55059 -0.25717  0.01098  0.14608  1.24117 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)   0.5349     0.1823   2.935  0.00438 **
    ## log1p(TP)     1.5771    20.0374   0.079  0.93747   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3724 on 78 degrees of freedom
    ## Multiple R-squared:  7.942e-05,  Adjusted R-squared:  -0.01274 
    ## F-statistic: 0.006195 on 1 and 78 DF,  p-value: 0.9375

``` r
plot(lm4)
```

![](Nutrients_files/figure-gfm/pre/post%20dreissenid%20analysis-4.png)<!-- -->

``` r
theme2 <- theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white"), legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white"))

p1 <- ggplot(data = king, aes(invade, lchltp, fill = invade)) + geom_boxplot() +
    ylab("Log Chla: Log TP") + scale_fill_viridis(discrete = TRUE, labels = c("Pre-invasion",
    "Post-invasion")) + labs(fill = "Invasion Status") + scale_x_discrete(labels = c("Pre",
    "Post")) + ylim(0, 0.75) + geom_signif(comparisons = list(c("pre", "post")),
    map_signif_level = TRUE, y_position = 0.6) + theme2 + ggtitle("Kingston")

p2 <- ggplot(data = brock, aes(invade, lchltp, fill = invade)) + geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, labels = c("Pre-invasion", "Post-invasion")) +
    labs(fill = "Invasion Status") + scale_x_discrete(labels = c("Pre", "Post")) +
    ylim(0, 0.75) + geom_signif(comparisons = list(c("pre", "post")), map_signif_level = TRUE,
    y_position = 0.6) + theme2 + ggtitle("Brockville")

chl.tp_plot <- ggarrange(p1, p2 + rremove("ylab"), nrow = 1, common.legend = TRUE,
    legend = "right")
ggsave("chla.tp_boxplot.svg", chl.tp_plot, width = 12, height = 8, dpi = 600)
chl.tp_plot
```

![](Nutrients_files/figure-gfm/dreissenid%20anaylsis%20plots-1.png)<!-- -->

``` r
p4 <- ggplot(data = king, aes(log1p(TP * 1000), log1p(CHLA))) + geom_point() + geom_smooth(method = lm,
    color = "blue", fill = "#69b3a2", se = TRUE) + theme(axis.title = element_text(size = 15),
    panel.background = element_rect(fill = "white")) + xlab("log TP (ug*L-1)") +
    ylab("log Chla (ug*L-1)") + ylim(0, 3.5) + ggtitle("Kingston") + facet_grid(. ~
    invade)

p5 <- ggplot(data = brock, aes(log1p(TP * 1000), log1p(CHLA))) + geom_point() + geom_smooth(method = lm,
    color = "blue", fill = "#69b3a2", se = TRUE) + theme(axis.title = element_text(size = 15),
    panel.background = element_rect(fill = "white")) + xlab("log TP (ug*L-1)") +
    ylab("log Chla (ug*L-1)") + ylim(0, 3.5) + ggtitle("Brockville") + facet_grid(. ~
    invade)

ggarrange(p4, p5 + rremove("ylab"), ncol = 2, common.legend = TRUE)
```

![](Nutrients_files/figure-gfm/dreissenid%20anaylsis%20plots-2.png)<!-- -->

``` r
uslr<-icefree_dreis %>% 
  filter(year<2001&year>1980)%>% 
    group_by(year,Site,Longitude,Latitude,River_section) %>%   # Grouping variable(s)
  summarise(
    CHLA = mean(CHLA),
    TDP = mean(TDP),       # calculate mean of column var in my.df
    SRP = mean(SRP),         
    TP = mean(TP),
    NH3.4 = mean(NH3.4),
    NO2.3 = mean(NO2.3),
    TN = mean(TN),
    lchltp=mean(lchltp,na.rm=TRUE)
    )

p5<-ggplot(data=uslr,aes(year,lchltp,color=River_section))+
  geom_point() +
  geom_smooth(method=loess,se=TRUE)+
  theme(axis.title = element_text(size=15),panel.background = element_rect(fill="white"))+
  xlab("Year")+
  ylab("Log Chla:TP")+
  geom_vline(xintercept=1990,color="red",size=1,linetype=5)
p5
```

![](Nutrients_files/figure-gfm/Annual%20trends-1.png)<!-- -->
