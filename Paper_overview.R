#' This R code contains the analysis and calculations for: 
#' Pesticides are substantially transported in particulate phase, driven by land use, rainfall event and pesticide characteristics 
#' – a runoff and erosion study in a small agricultural catchment. M.C. Commelin et al. (2022)
#' The data on which this research is build can be found at:
#' Data that comes first in the paper also comes first in this code starting from the introduction.

# Initialization ----------------------------------------------------------
# load packages
library(aomisc)
library(tidyverse)
library(lubridate)
library(cowplot)
library(zoo)
library(scales)
library(sf)
library(ggrepel)
library(hms)
library(extrafont)
library(ggforce)
library(broom)
library(pander)
library(ggpubr)
library(rstatix)

# A specific version of Rttf2pt1 is needed to run extrafont
#remotes::install_version("Rttf2pt1", version = "1.3.8")

#save version of packages use to run the code - add to readme.md with pander
pander(sessionInfo())

# Functions ----------------------------------------------------------------
# discharge functions
source("sources/discharge_calculations.R")

#' Give all words in a character column CAPS at start of word.
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

nl_tz <- locale(tz = "Etc/GMT-1")

# Catchment and landuse ------------------------------------
#' For privacy reasons the gpx data is removed from this data. For the results of this research no spatial data is needed.
#' If spatial data is needed for further research, please contact the author.
# load field data
fields <- read_csv("data/fields.csv")
fields_19 <- filter(fields, year == 2019)
fields_20 <- filter(fields, year == 2020)

#' field numbers for three main fields.
fields_large <- tibble(large = c(rep("C", 9), rep("A", 5), rep("B", 7)),
                       field_nr = c(1009, 1020, 1013, 2024, 2019, 2023, 3024, 3012, 3022, 
                                    1021, 2017, 2009, 3014, 3015,
                                    1007, 1011, 2003, 2013, 3019, 3016, 3023))
# soil sample data
soil_samp <- read_csv("data/Soil_samples.csv", lazy = F)

#' CH 2.1 - Summary for study area
# area catchment
filter(fields, year == 9999)$area
landuse_area <- fields_19 %>%
  group_by(crop_type) %>%
  summarise(area = sum(area))
# arable landuse area
arable_area <- landuse_area %>%
  filter(!str_detect(crop_type, "Grass|Appel|Conifer")) %>%
  summarise(area = sum(area))
# area of orchard
apple_area <- landuse_area %>%
  filter(str_detect(crop_type, "Appel")) %>%
  summarise(area = sum(area))
# area extensive grassland
grass_area <- landuse_area %>%
  filter(str_detect(crop_type, "Grass")) %>%
  summarise(area = sum(area))
#' CH 2.1 - data on texture (sand, silt, clay) and pH and OM.
#' load texture data, in micro meter
tex_id <- read_csv("data/tex_range_id.csv")
tex_data <- read_csv("data/texture_data.csv")

#'calculate cumsum for each sample on transposed tibble.
tex_cum <- as_tibble(as.matrix(t(tex_data[-(1:2)]))) %>%
  rename_with(~ tex_data$tex_code) %>%
  mutate(across(everything(), ~ cumsum(.)))
tex_d50 <- tex_cum %>%
  summarise(across(everything(), ~ approx(., tex_id$upper, xout = 50)$y))
tex_d90 <- tex_cum %>%
  summarise(across(everything(), ~ approx(., tex_id$upper, xout = 90)$y))
tex_clay <- tex_cum %>%
  summarise(across(everything(), ~ approx(tex_id$upper, .,  xout = 2)$y))
tex_silt <- tex_cum %>%
  summarise(across(everything(), ~ approx(tex_id$upper, .,  xout = 50)$y))
tex_sand <- tex_cum %>%
  summarise(across(everything(), ~ max(.)))
tex_summary <- t(bind_rows(tex_d50, tex_d90, tex_clay, tex_silt, tex_sand)) %>%
  as_tibble() %>%
  setNames(., c("d50", "d90", "clay", "silt", "sand"))  %>%
  mutate(sand = sand - silt,
         silt = silt - clay,
         tex_code = tex_data$tex_code,
         date = tex_data$date) %>%
  filter(date == "2019-08-29" | date == "2020-03-04") # select catchment
# calculate mean, and se for d50, d90, sand, silt, clay
means_tex <- tex_summary %>%
  summarise(across(everything(), ~ mean(.)))
se_tex <- tex_summary %>%
  summarise(across(everything(), ~ sd(.)/sqrt(length(.))))

#' OM and pH overview
catch_samp_an <- read_csv("data/catchment_sample_analysis.csv") %>%
  filter(!is.na(field_nr))
# mean and SE of mean
mean(catch_samp_an$pH)
mean(catch_samp_an$OM)
sd(catch_samp_an$pH)/sqrt(nrow(catch_samp_an))
sd(catch_samp_an$OM)/sqrt(nrow(catch_samp_an))

#' CH 2.1 - long term rainfall and temperature for KNMI.
#' this is secondary data which can be downloaded from:
#' https://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/maandgegevens/mndgeg_380_tg.txt
#' https://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/maandgegevens/mndgeg_380_rh24.txt
temp_knmi <- read_csv("ext_data/knmi_ma_t.csv") %>%
  filter(yyyy <= 2020 & yyyy >= 1991)
mean(temp_knmi$year)/10   # mean annual temperature
sd(temp_knmi$year)/10     # sd of mean temperature

p_knmi <- read_csv("ext_data/knmi_ma_p.csv") %>%
  filter(yyyy <= 2020 & yyyy >= 1991)
mean(p_knmi$year)/10      # mean annual precipitation
sd(p_knmi$year)/10        # sd of mean P

#' CH 2.2 - number of soil samples
count(distinct(filter(soil_samp, !is.na(pest_ID)), pest_ID))

# Event selection ---------------------------------------------------------
ev_samples <- read_csv("data/Event_samples.csv")
#make dataframe with data for each event
ev_dates <- ev_samples %>%
  mutate(date = date(timestamp)) %>%
  group_by(date) %>%
  mutate(start_t = min(timestamp) - hours(1),
         end_t = max(timestamp) + hours(3)) %>%
  select(date, start_t, end_t) %>%
  distinct() %>%
  ungroup() %>%
  mutate(t_diff = if_else(start_t > lag(end_t), 1, 0)) %>%
  filter(!is.na(date)) %>%
  filter(t_diff == 1 | is.na(t_diff)) %>%
  select(-t_diff)
#write_csv(ev_dates, "ev_dates.csv")

#' this are the final used start and end times, which fit the hydrographs for each 
#' individual event in more detail.
ev_dates <- read_csv("sources/ev_dates.csv")

# Calculate discharge -----------------------------------------------------

# load data
wb_data <- read_csv("data/WB_data.csv") # data from waterboard (Ha)
cr6 <- read_csv("data/CR6.csv", col_types = "TddddddddddddddD") # from cr6 (Hb)
crest_nap <- 79.331 # m + NAP as measured by waterboard

# use Q_calc to calculate discharge for each event and store in list
# the dicharge calculations can be found in the 'sources' folder for more details.
Q_tcor <- Q_calc(ev_dates, wb_data, cr6, crest_nap)
df_dat <- Q_tcor[[1]] # data with discharge for each event
df_tcor <- Q_tcor[[2]] # based on synchronization between different time series a time correction is performed

#' CH 2.2 - runoff duration 
#' A count based on the event graphs made (not in this code), gives 20 out of the 25 observed events which
#' are covered for the full peak within 72 minutes. Difficult to catch this easily
#' with a calculated number.
#' CH 2.2 - calculate number of days with runoff events in 2019 and 2020.
runoff_days <- wb_data %>%
  filter(year(timestamp) == 2019 | year(timestamp) == 2020) %>%
  filter(Q >= 0.020) %>%
  mutate(date = date(timestamp)) %>%
  group_by(date) %>%
  mutate(count = n()) %>%
  distinct(date, count) %>%
  filter(count >= 10)
length(runoff_days$date)

# analysis of Qmax and Qtotal for the 39 observed runoff events.
runoff_analysis <- wb_data %>%
  mutate(date = date(timestamp)) %>%
  semi_join(runoff_days, by = "date") %>%
  arrange(timestamp) %>%
  mutate(secs = as.numeric(lead(timestamp) - timestamp),
         Ha = Wh - crest_nap,
         Q = Qparshall(Ha, 0.609),
         Q_m3 = Q * secs) %>%
  group_by(date) %>%
  filter(Q >= 0.015) %>%
  summarise(Q_max = max(Q),
            Q_m3= sum(Q_m3, na.rm = T),
            D = sum(secs) / 60) %>%
  mutate(group = "All observed")
  
# Precipitation data ------------------------------------------------------
#' this data is provided as secondary source, it is clipped from the 5 minute radar
#' of KNMU for 2018 - 2020.
knmi_p <- read_csv("ext_data/KNMI_rain.csv") %>%
  group_by(timestamp) %>%
  summarize(mean(P))
names(knmi_p) <- c("timestamp", "P_m_knmi")
knmi_p <- mutate(knmi_p, P_i_knmi = P_m_knmi * 12)
#add rainfall from knmi to event data, also keep rainfall by cr6
name <- str_remove_all(ev_dates$date, "-")
for (i in seq_along(name)) {
  df_dat[[i]] <- left_join(df_dat[[i]], knmi_p, by = "timestamp") %>%
    mutate(P_i_cr6 = rain_tot * 60)
}

# Erosion data ------------------------------------------------------------
# add tss estimate for 20200812 (made in pesticide calculations.R)
# which was not available due to loss of samples.
tss_est <- read_csv("sources/tss_20200812.csv", lazy = FALSE) %>%
  left_join(ev_samples, by = "pest_ID") %>%
  select(-sed_conc.y) %>%
  rename(sed_conc = sed_conc.x)

# add tss estimate to the event samples data
ev_samples <- ev_samples %>%
  anti_join(tss_est, by = "pest_ID") %>%
  bind_rows(tss_est)

# make sediment runoff data for each event
df_er <- vector("list", length = length(name))
for (i in seq_along(name)) {
  ev_start <- ev_dates$start_t[i]
  ev_end <- ev_dates$end_t[i]
  df_er[[i]] <- ev_samples %>%
    filter(timestamp > ev_start & timestamp < ev_end) %>%
    select(timestamp, sed_conc) %>%
    mutate(timestamp = timestamp + df_tcor[[i]]) %>%
    group_by(timestamp) %>%
    summarise(sed_conc = mean(sed_conc)) %>%
    left_join(df_dat[[i]], by = "timestamp") %>%
    mutate(load = sed_conc * Q_int) %>%
    select(timestamp, sed_conc, load)
}

# Total Qwat, Qsed, rain --------------------------------
## Fig 2 calculate --------------------------
df_totals <- vector("list", length = length(name))
for (i in seq_along(name)) {
Q_tot <- df_dat[[i]] %>%
  arrange(timestamp) %>%
  mutate(secs = as.numeric(lead(timestamp) - timestamp)*60,
         Q_m3 = Q_int * secs) %>%
  full_join(df_er[[i]], by = "timestamp") %>%
  arrange(timestamp) %>%
  mutate(sed_int = na.approx(sed_conc, maxgap = 10),
         Sed_kg = sed_int * Q_m3) %>%
  mutate(Q_m3 = if_else(Q_int > 0.012, Q_m3, 0),
         Q_duration = if_else(Q_int > 0.012, 1, 0)) %>%
  summarise(Q_max = max(Q_int, na.rm = T),
            Q_m3 = sum(Q_m3, na.rm = T),
            Sed_kg = sum(Sed_kg, na.rm = T),
            Sed_conc_max = max(sed_int, na.rm = T),
            P = sum(P_m_knmi, na.rm = T),
            Pi_max = max(P_i_knmi, na.rm = T),
            D = sum(P_i_knmi > 0.5, na.rm = T)*5,
            Temp_cr6 = mean(temp, na.rm = T),
            Q_duration = sum(Q_duration > 0)) %>%
  mutate(undisc_w = 2.5 * Q_m3)
df_totals[[i]] <- Q_tot
}

# combine list into tibble with each event in a row
# CH 3 - mean TSS for events 
total <- bind_rows(df_totals) %>% 
  bind_cols(ev_dates) %>%
  mutate(Q_mm = Q_m3 / 38 / 10,
         QP = Q_mm / P * 100,
         Skg_ha = Sed_kg / 38,
         evi = (Pi_max * P)/D,
         mean_tss= Sed_kg / Q_m3)

# store data for runoff overview visualization
erosion25 <- total

# CH 3.1 - analysis of relation between rain and erosion etc.
fit_total <- total[-17, ]   # remove event 16-08-2020 because this is an outlier
ggplot(fit_total) + geom_point(aes(x = Q_m3, y = P))
fit_water <- lm(Q_m3 ~ P + D, data = fit_total)
summary(fit_water)

# no clear relation - very low r2
fit_sed <- lm(Sed_kg ~ Pi_max, data = fit_total)
summary(fit_sed)
ggplot(fit_total) + geom_point(aes(x = Sed_kg, y = evi))

# Pesticide data ----------------------------------------------------------

#' Load all pesticide concentrations and application data 
#' Pesticide concentration in W and S of all samples, calculations from raw data
#' are performed in pesticide_calculations.R
lc_all_data <- read_csv("sources/lc_all_data.csv", na = c("-999", "-888", "-777", "NA"))
# Pesticide application data
pest_appl <- read_csv("data/Pest_application.csv") %>%
  left_join(fields_large, by = "field_nr")
# description of applied pesticides - link to active ingredients
pest_desc <- read_csv("data/Pest_description.csv") %>%
  mutate(ai = capwords(ai),
         ai_dose = if_else(ai_dose_unit == "%", ai_dose * 10, ai_dose),
         ai = if_else(ai == "Glyphosate-isopropylammonium", "Glyphosate", ai)) %>%
  rename(compound = ai)
# chemical characteristics of active ingredients
compound_char <- read_csv("ext_data/AI_characteristics.csv") %>%
  mutate(koc_class = if_else(koc < 15, "very mobile", "mobile"),
         koc_class = if_else(koc > 75 & koc < 500, "moderately mobile", koc_class),
         koc_class = if_else(koc > 500 & koc < 4000, "slightly mobile", koc_class),
         koc_class = if_else(koc > 4000, "non-mobile", koc_class),
         dt50_class = if_else(dt50 < 30, "non-persistent", "moderateley persistent"),
         dt50_class = if_else(dt50 > 100 & dt50 < 365, "persistent", dt50_class),
         dt50_class = if_else(dt50 > 365, "very persistent", dt50_class),
         Sw_class = if_else(Sw < 50, "Low", "Moderate"),
         Sw_class = if_else(Sw > 500, "High", Sw_class))

# analysed AS
analysed <- str_remove(str_subset(names(lc_all_data), "conc"), "conc_")

# remove terbuthylazine from analysis, is not applied or detected, we added it to analysis
# because it is often detected on fields, however does not have a function in this research
analysed_comp <- tibble(compound = analysed) %>%
  filter(compound != "Terbuthylazine")

#' due to point source pollution, the Glyphosate and AMPA results for 3 events are removed.
point_src_dates <- c("2020-09-23", "2020-10-08", "2020-10-22", "2020-09-26") 
point_src_compounds <- c("Glyphosate", "AMPA")

## Methods and analysis -------------------------
#' AI analysed
ai_analysed <- analysed_comp %>%
  mutate(analysed = rep(1, nrow(analysed_comp))) %>%
  mutate(compound = str_remove(compound, "conc_")) %>%
  left_join(compound_char, by = "compound") %>%
  select(compound, analysed, type)

# CH 2.4 - number of compounds + type
ai_summary <- ai_analysed %>%
  mutate(analysed = sum(analysed)) %>%
  group_by(type) %>%
  summarise(analysed = max(analysed),
            count = n())

ai_char_range <- compound_char %>%
  semi_join(ai_analysed, by = "compound") %>%
  summarize(koc = str_c(min(koc, na.rm = T), " - ", max(koc, na.rm = T)),
            dt50 = str_c(min(dt50, na.rm = T), " - ", max(dt50, na.rm = T)))
## Catchment data --------------------------------
nms <- str_subset(names(lc_all_data), "conc_")

pest_catch <- lc_all_data %>%
  filter(samp_type == "C")%>%
  left_join(soil_samp, by = "pest_ID") %>%
  left_join(fields_large, by = "field_nr") %>%
  select(pest_ID, field_nr, sub_loc, everything()) %>%
  pivot_longer(starts_with("conc_"),
               names_to = "compound", values_to = "conc") %>%
  mutate(compound = str_remove(compound, "conc_"),
         conc = if_else(conc == -999, 0, conc)) %>%
  semi_join(ai_analysed, by = "compound") %>%
  left_join(compound_char, by = "compound")

## Outlet data ------------------------------------
#' pesticide data for each event
pest_ts <- ev_samples %>%
  mutate(date = date(timestamp),
         sed_gr = if_else(is.na(sed_gr), w_d_cup - w_e_cup, sed_gr)) %>%
  select(date, pest_ID, sed_gr, timestamp) %>%
  group_by(timestamp) %>%
  summarise_all(mean, na.rm = T) %>%
  group_by(pest_ID) %>%
  summarise_all(mean, na.rm = T) %>%
  distinct() %>%
  filter(!is.na(pest_ID)) %>%
  select(-timestamp)
# store point source dates for GLY and AMPA in separate file and remove from analysis
lc_all <- lc_all_data %>%
  left_join(pest_ts, by = "pest_ID") %>%
  filter(!is.na(date)) %>%
  select(-an_name, -pest_ID, -week, -conc_Terbuthylazine) %>%
  replace(is.na(.), 0)
pnt_src_lc <- lc_all %>%
  filter(date > "2020-09-10") %>%
  select(samp_type:sed_gr)
lc_all <- lc_all %>%
  mutate(conc_Glyphosate = if_else(date > "2020-09-10", 0, conc_Glyphosate),
         conc_AMPA = if_else(date > "2020-09-10", 0, conc_AMPA))

lc_sed <- lc_all %>%
  filter(samp_type == "S")%>%
  group_by(date, samp_type) %>%
  summarise_at(vars(starts_with('conc_')), ~weighted.mean(., sed_gr, na.rm = T))
lc_all <- lc_all %>%
  select(-sed_gr) %>%
  filter(samp_type == "W") %>%
  group_by(date, samp_type) %>%
  summarise_at(vars(starts_with('conc_')), ~mean(., na.rm = T)) %>%
  bind_rows(lc_sed) %>%
  mutate(date = date(date))

#' add the pesticide data to the totals table. Calculate total 
#' mass transported for all compounds in water and sed.
df_pest <- vector("list", length = length(name))
for (i in seq_along(name)) {
  ev_date <- ev_dates$date[i]
  df_pest[[i]] <- lc_all %>%
    filter(date == ev_date)
  nms <- names(df_pest[[i]])
  oob_lim <-  5 * df_totals[[i]]$undisc_w[1] #limit to end barplots, so the lower values are still vissible
  df_pest[[i]] <- pivot_longer(df_pest[[i]], nms[3:length(nms)],
                                names_to = "compound", values_to = "conc") %>%
     mutate(compound = str_replace(compound, "conc_", ""),
            tot_m = if_else(samp_type == "S", conc * 0.001 * df_totals[[i]]$Sed_kg, conc * df_totals[[i]]$Q_m3),
            val_oob = if_else(tot_m > oob_lim, 1, 0),
            label_oob = signif(if_else(val_oob == 1, tot_m, NaN), digits = 4))
  df_totals[[i]] <- df_totals[[i]] %>%
    mutate(TP_sed = sum(filter(df_pest[[i]], samp_type == "S")$tot_m),
           TP_wat = sum(filter(df_pest[[i]], samp_type == "W")$tot_m),
           P_sed_n = length(filter(df_pest[[i]], samp_type == "S" & tot_m > 0)$tot_m),
           P_wat_n = length(filter(df_pest[[i]], samp_type == "W" & tot_m > 0)$tot_m),
           TP_wat_loq = (P_sed_n - P_wat_n) * undisc_w + TP_wat) 
}
#'remove dates that have wrong data, or must be excluded from overview
dates_rem <- tibble(date = date(c("2019-03-15", "2019-08-18", "2020-09-26")))
total_pest <- bind_rows(df_totals) %>% bind_cols(ev_dates) %>%
  mutate(date = date(date)) %>%
  filter(P_sed_n != 0 | P_wat_n != 0) %>%
  anti_join(dates_rem, by = "date") %>%
  mutate(number = c(1:14)) %>%
  select(date, TP_sed, TP_wat, P_sed_n, P_wat_n, TP_wat_loq, number)
# add pesticide data to total overview table
total <- left_join(total_pest, erosion25) %>%
  mutate(date_lab = format(date, "%d-%m-%Y"),
         date_lab = if_else(date > "2020-09-10", str_c("*", date_lab), date_lab),
         PP_rat = TP_sed / (TP_sed + TP_wat),
         conc_DP = TP_wat / (Q_m3 * 1000),
         conc_PP = TP_sed / Sed_kg,
         fact_phase = conc_PP / conc_DP,
         Tot_P = (TP_wat + TP_sed) / 1000,
         mean_p = ((TP_wat + TP_sed) * 1000) / ((Q_m3 * 1000) + (Sed_kg / 2.65)))

# 3.1 Runoff, erosion and pesticide transport ------------------------

# TSS for the two events from 'total$mean_tss'
mean(total$mean_tss)   # mean tss for all events
# sum TSS and water
sum(total$Q_m3)        # total water discharges
sum(total$Sed_kg)/1000 # total erosion
sum(total$Q_m3)*0.05   # uncertainty in water discharge
sum(total$Sed_kg)/1000 * 0.06   # uncertainty in erosion
5800/79                # ratio Qsed:Qwat
(sum(total$Sed_kg)/1000) / (arable_area[1,1]/10000) # erosion t/ha

# CH 3.1 - check relation between tot PP and erosion
fit_total <- total[-10, ]   # remove event 16-08-2020 because this is an outlier
fit_pest <- lm(TP_sed ~ Sed_kg, data = fit_total)
summary(fit_pest)
ggplot(total) + geom_point(aes(x = mean_tss, y = Pi_max))

#' CH 3.1 - total pesticide and runoff losses over 14 events
pest_total <- tibble(S = sum(total$TP_sed, na.rm = T) / 1000, # g -pesticides in water
                     W = sum(total$TP_wat, na.rm = T) / 1000, # g - pesticides in TSS
                     sum = S + W,
                     tot_qw = sum(total$Q_m3),             # m3 - discharged water
                     tot_qs = sum(total$Sed_kg) / 1000,    # ton - discharged sediment
                     mean_pw = W / tot_qw * 1000,  # ug L - mean load pw
                     mean_ps = S / tot_qs * 1000,  # ug kg - mean load ps
                     mean_p = (S + W) / (tot_qw + (tot_qs/2.65)) * 1000)  # overall mean load

### Fig 2 calculate -------------------------------
# set unit of pesticides to grams
pest_tot <- total %>%
  select(number, TP_sed, TP_wat, TP_wat_loq) %>%
  pivot_longer(TP_sed:TP_wat, names_to = "samp_type", values_to = "P") %>%
  mutate(TP_wat_loq = if_else(samp_type == "TP_sed", 0, TP_wat_loq / 1000),
         P = P / 1000,
         error = if_else(samp_type == "TP_sed", 0.18 * P, 0.17 * P)) #relative uncertainty

### Fig 3 calculate ----------------------------------
#' add pesticide types to each compound, make compounds a factor in order of type and descending 
#' mean per compound - to get the right order in the boxplot.
oob_fig3_s <- 15
oob_fig3_w <- 0.13
comp_box <- bind_rows(df_pest) %>%
  semi_join(total, by = "date") %>%
  left_join(compound_char, by = "compound") %>%
  group_by(compound, samp_type) %>%
  mutate(mean_c = mean(conc, na.rm = T)) %>%
  ungroup() %>%
  mutate(other = if_else(type == "herbicide", str_c(compound, "^~H"), compound),
         other = if_else(type == "insecticide", str_c(other, "^~I"), other),
         other = if_else(type == "fungicide", str_c(other, "^~F"), other),
         other_s = if_else(conc < 101, "AS < 0.1~mg~kg^-1", other),
         other_w = if_else(conc < 5, "AS < 0.005~mg~L^-1", other),
         conc = round(conc / 1000, digits = 2),
         type_s = if_else(other_s == "AS < 0.1~mg~kg^-1", "other", type),
         type_w = if_else(other_w == "AS < 0.005~mg~L^-1", "other", type),
         label_oob_s = if_else(conc > oob_fig3_s, conc, NULL),
         label_oob_w = if_else(conc > oob_fig3_w & samp_type == "W", 
                               as.character(signif(conc), digits = 2), NULL),
         label_oob_s = as.character(round(label_oob_s)),
         date = as.character(date)) %>%
  arrange(type_s, desc(mean_c), other_s)
level_fig3_s <- unique(comp_box$other_s)
level_fig3_w <- unique(comp_box$other_w)
comp_box$other_s <- factor(comp_box$other_s, levels = level_fig3_s)
comp_box$other_w <- factor(comp_box$other_w, levels = level_fig3_w)

# CH 3.1 - number of compounds detected
str_c(min(total$P_sed_n), " - ", max(total$P_sed_n))
str_c(min(total$P_wat_n), " - ", max(total$P_wat_n))

# CH 3.1 - total different compounds detected in TSS and water
comp_overview <- bind_rows(df_pest) %>%
  semi_join(total, by = "date") %>%
  mutate(year = year(date)) %>%
  group_by(compound, samp_type, year) %>%
  summarise(tot_m = sum(tot_m, na.rm = T),
            freq = sum(conc > 0, na.rm = T)) %>%
  left_join(compound_char, by = "compound") %>%
  group_by(compound) %>%
  mutate(mean_c = mean(tot_m, na.rm = T))
comp_tot_n <- comp_overview %>%
  group_by(samp_type, compound) %>%
  summarise(freq = sum(freq > 0)) %>%
  group_by(samp_type) %>%
  filter(freq > 0) %>%
  distinct(compound) %>%
  ungroup()
# LOQ etc not in this code...

max(total$TP_wat/total$Q_m3)    # max mean events concentration
sum((total$TP_wat/total$Q_m3) > 0.5) # events exceeding 0.5 ug/L limit

# 3.2 Chemical characteristics and pesticide transport ------------------------

## Transport over time -----------------
#'first make a table with transport and related data for each AS
AS_trans_total <- total %>%
  select(date, Q_m3)

AS_transport_indiv <- bind_rows(df_pest) %>%
  anti_join(dates_rem, by = "date") %>%
  select(-(val_oob:label_oob)) %>%
  pivot_wider(names_from = samp_type, values_from = c(conc, tot_m)) %>%
  filter(conc_W > 0 | conc_S > 0) %>%
  left_join(AS_trans_total, by = "date") %>%
  mutate(tot_m_W_loq = if_else(conc_W == 0 & conc_S > 0, 0.8 * 2.5 * Q_m3, tot_m_W),
         conc_W = if_else(conc_W == 0, 1.25, conc_W), #adjust 0 values to LOQ/2 for Kp calculation
         conc_S = if_else(conc_S == 0, 5, conc_S),
         Kp = conc_S / conc_W,
         PP_rat = tot_m_S / (tot_m_S + tot_m_W),
         PP_type = if_else(PP_rat < 0.5, 0, 2),
         PP_type = if_else(PP_rat > 0.5, 2, PP_type),
         year = year(date))
AS_transport <- AS_transport_indiv %>%
  group_by(compound) %>%
  summarise(KP_mean = mean(Kp),
            Kp_sd = sd(Kp),
            PP_load = sum(tot_m_S),
            DP_load = sum(tot_m_W),
            DP_LOQ_load = sum(tot_m_W_loq),
            DP_max = sum(PP_type == 0),
            PP_max = sum(PP_type == 2),
            no_max = sum(PP_type == 1),
            n = n(),
            PP_rat_all = PP_load / (PP_load + DP_load)) %>%
  ungroup() %>%
  mutate(PP_type = if_else(PP_rat_all < 0.4, 0, 1),
         PP_type = if_else(PP_rat_all > 0.7, 2, PP_type)) %>%
  arrange(PP_type) %>%
  left_join(compound_char, by = "compound")

#write_csv(AS_transport, "AS_transport.csv")

## Table S.8 ---------------
#' make table with relation between catchment application & concentration and 
#' observed concentrations at the outlet during events
field_appl_all <- fields %>%
  select(field_nr, large, year, crop_type, area) %>%
  rename(area_m2 = area) %>%
  left_join(pest_appl, by = c("field_nr", "large")) %>%
  left_join(pest_desc, by = "pesticide") %>%
  rename(appl_date = date) %>%
  filter(!(is.na(compound))) %>%
  mutate(ai_dose_unit = if_else(dose_unit == "l/ha", "g/l", "g/kg"),
         gram_ha = ai_dose * dose) %>%
  semi_join(df_pest[[1]], by = "compound")

df_relation <- vector("list", length = length(name))
for (i in seq_along(name)) {
  if (nrow(df_pest[[i]]) > 0) {
    relation <- df_pest[[i]] %>%
      select(-(tot_m:label_oob)) %>%
      pivot_wider(names_from = "samp_type", values_from = "conc", 
                  names_prefix = "conc_")
    AMPA_pest <- relation %>%
      filter(compound == "AMPA")
    relation <- relation %>%
      left_join(field_appl_all, by = "compound") %>%
      mutate(daa = as.numeric(days(date - appl_date))/86400,
             daa = ifelse(daa < 0, 9999, daa)) %>%
      group_by(compound, large) %>%
      mutate(daa_min = pmin(daa, na.rm = TRUE)) %>%
      ungroup() %>%
      select(date, conc_W, conc_S, compound, large, daa_min, gram_ha, appl_date) %>%
      #filter(daa_min != 9999) %>%
      distinct(date, compound, daa_min, .keep_all = TRUE) %>%
      arrange(desc(conc_S)) %>%
      left_join(compound_char, by = "compound")
    AMPA_add <- relation %>%
      filter(compound == "Glyphosate") %>%
      select(large:dt50_class) %>%
      bind_cols(AMPA_pest)
    df_relation[[i]] <- bind_rows(relation, AMPA_add)
  }
}

# First calculate for events with field campaign
pos_i <- c(2, 3, 15, 16, 17)
date_catch <- c("2019-05-11", "2019-05-28", "2020-08-20", "2020-08-20", "2020-08-20")
df_outlet_tab <- vector("list", length = length(date_catch))
#view(df_relation[[i]])
for (i in seq_along(date_catch)) { 
  # table with detection on 28-may
  pest_cath_28 <- pest_catch %>%
    filter(date == date_catch[i]) %>%
    filter(!is.na(large)) %>%
    select(date, large, compound, conc, type) %>%
    group_by(date, large, type, compound) %>%
    summarise(conc = signif(mean(conc, na.rm = T), digits = 2)) %>%
    filter(!is.na(conc)) %>%
    ungroup() %>%
    pivot_wider(names_from = large, values_from = conc) %>%
    select((compound:C))
  j = pos_i[i]
  outlet28 <- df_relation[[j]] %>%
    group_by(compound) %>%
    slice_min(order_by = daa_min) %>%
    filter(conc_W != 0 | conc_S != 0) %>%
    left_join(pest_cath_28, by = "compound") %>%
    select(date, compound, conc_W, conc_S, A, B, C, everything()) %>%
    arrange(desc(conc_S)) %>%
    rename(last_applied = large,
           DBE = daa_min) %>%
    select((date:appl_date))
  df_outlet_tab[[i]] <- outlet28
}

outlet_table_all <- bind_rows(df_outlet_tab)

# add events without field campaign data
df_outlet_rest <- vector("list", length = length(name)) 

for (i in seq_along(name)) {
  if (nrow(df_pest[[i]]) > 0) {
    df_outlet_rest[[i]] <- df_relation[[i]] %>%
      group_by(compound) %>%
      slice_min(order_by = daa_min) %>%
      filter(conc_W != 0 | conc_S != 0) %>%
      select(date, compound, conc_W, conc_S, everything()) %>%
      arrange(desc(conc_S)) %>%
      rename(last_applied = large,
             DBE = daa_min) %>%
      select((date:appl_date))
  }
}
outlet_rest <- bind_rows(df_outlet_rest) %>%
  anti_join(outlet_table_all, by = "date")

# combine to 1 large table
outlet_all_dates <- bind_rows(outlet_table_all, outlet_rest) %>%
  arrange(date)

#write_csv(outlet_all_dates, "outlet_all_dates.csv")

#' Add application data to relate transport to field management
large_area <- fields %>%
  select(field_nr, area) %>%
  distinct()

appl_rate_large <- left_join(pest_appl, pest_desc, by = "pesticide") %>%
  semi_join(analysed_comp, by = "compound") %>%
  select(-area) %>%
  mutate(ai_appl_dose = dose * ai_dose,
         year = year(date)) %>%
  left_join(large_area, by = "field_nr") %>%
  mutate(ai_appl_large = as.numeric(ai_appl_dose * (area / 10000))) %>%
  group_by(compound, large, year) %>%
  summarise(appl_rate = sum(ai_appl_large),
            n_appl = n())
appl_rate_AS <- appl_rate_large %>%
  group_by(compound) %>%
  summarise(appl_rate = sum(appl_rate),
            n_appl = sum(n_appl))
AS_transport <- AS_transport %>%
  left_join(appl_rate_AS, by = "compound") %>%
  mutate(trans_DP = (DP_load / 1000) / appl_rate,
         trans_PP = (PP_load / 1000) / appl_rate,
         trans_tot = trans_DP + trans_PP)
appl_tot_large_analysed <- appl_rate_large %>%
  group_by(large) %>%
  summarise(appl_rate = sum(appl_rate),
            n_appl = sum(n_appl))
analysed_appl_all <- sum(appl_tot_large_analysed$appl_rate)


DBE_table <- outlet_all_dates %>%
  select(date, compound, DBE)
AS_transport_rate <- AS_transport_indiv %>%
  left_join(appl_rate_large, by = c("compound", "year")) %>%
  mutate(trans_DP = (tot_m_W / 1000) / appl_rate,
         trans_PP = (tot_m_S / 1000) / appl_rate,
         trans_tot = trans_DP + trans_PP) %>%
  left_join(DBE_table, by = c("date", "compound")) %>%
  distinct(date, compound, .keep_all = T)
# calculate the fraction of total transport within a DBE range
# DP for 5, 10, 50 days
sum(filter(AS_transport_rate, DBE < 5)$tot_m_W) / sum(AS_transport_rate$tot_m_W)
sum(filter(AS_transport_rate, DBE < 10)$tot_m_W) / sum(AS_transport_rate$tot_m_W)
sum(filter(AS_transport_rate, DBE < 50)$tot_m_W) / sum(AS_transport_rate$tot_m_W)
# PP for 5, 10, 50, 100 days
sum(filter(AS_transport_rate, DBE < 5)$tot_m_S) / sum(AS_transport_rate$tot_m_S)
sum(filter(AS_transport_rate, DBE < 10)$tot_m_S) / sum(AS_transport_rate$tot_m_S)
sum(filter(AS_transport_rate, DBE < 50)$tot_m_S) / sum(AS_transport_rate$tot_m_S)
sum(filter(AS_transport_rate, DBE < 100)$tot_m_S) / sum(AS_transport_rate$tot_m_S)

## Table 1 -------------------
# dt50
appl_dt50 <- appl_rate_AS %>%
  left_join(compound_char, by = "compound") %>%
  group_by(dt50_class) %>%
  summarise(m_appl = sum(appl_rate),
            n_appl = n())
mass_dt50 <- AS_transport %>%
  filter(compound != "AMPA") %>%
  group_by(dt50_class) %>%
  summarise(n_trans = n(),
            m_PP = sum(PP_load),
            m_DP = sum(DP_load)) %>%
  ungroup() %>%
  left_join(appl_dt50, by = "dt50_class") %>%
  mutate(trans_rate = (m_PP + m_DP) / (sum(m_PP) + sum(m_DP)),
         appl_rate = m_appl / sum(m_appl)) 

# koc
appl_koc <- appl_rate_AS %>%
  left_join(compound_char, by = "compound") %>%
  group_by(koc_class) %>%
  summarise(m_appl = sum(appl_rate),
            n_appl = n())
mass_koc <- AS_transport %>%
  filter(compound != "AMPA") %>%
  group_by(koc_class) %>%
  summarise(n_trans = n(),
            m_PP = sum(PP_load),
            m_DP = sum(DP_load)) %>%
  ungroup() %>%
  left_join(appl_koc, by = "koc_class") %>%
  mutate(trans_rate = (m_PP + m_DP) / (sum(m_PP) + sum(m_DP)),
         appl_rate = m_appl / sum(m_appl))

# Sw
appl_Sw <- appl_rate_AS %>%
  left_join(compound_char, by = "compound") %>%
  group_by(Sw_class) %>%
  summarise(m_appl = sum(appl_rate),
            n_appl = n())
mass_Sw <- AS_transport %>%
  filter(compound != "AMPA") %>%
  group_by(Sw_class) %>%
  summarise(n_trans = n(),
            m_PP = sum(PP_load),
            m_DP = sum(DP_load)) %>%
  ungroup() %>%
  left_join(appl_Sw, by = "Sw_class") %>%
  mutate(trans_rate = (m_PP + m_DP) / (sum(m_PP) + sum(m_DP)),
         appl_rate = m_appl / sum(m_appl))
# save tables 
# write_csv(mass_dt50, "table_dt50.csv")
# write_csv(mass_koc, "table_koc.csv")
# write_csv(mass_Sw, "table_sw.csv")

#'check if all compounds in runoff are also detected on fields
comp_runoff <- comp_tot_n %>%
  mutate(compound = as.character(compound)) %>%
  select(compound) %>%
  distinct(compound)

#' range koc of AS in DP
DP_AS_range <- AS_transport %>%
  filter(DP_load > 0) %>%
  select(koc) %>%
  summarise(koc = str_c(min(koc), " - ", max(koc)))

#' compound only in PP
PP_no_DP <- outlet_all_dates %>%
  filter(DBE < 11) %>%
  filter(conc_S > 0 & conc_W < 0.1) %>%
  distinct(compound) %>%
  left_join(compound_char, by = "compound")

#' dt50 of AS that are not transported
dt50_no_transport <- anti_join(ai_analysed, comp_runoff, by = "compound") %>%
  left_join(compound_char) %>%
  select(compound, dt50)

## Fig 4 calculate -----------------------
#comp_labels <- read_csv("comp_labels.csv")
dt50_compounds <- c("Mandipropamid", "Metobromuron", "Metribuzin", "Prosulfocarb", "Epoxiconazole", "Bixafen")
dt50_compounds_w <- c("Mandipropamid", "Metobromuron", "Metribuzin")
dt50_compounds_s <- c("Prosulfocarb", "Epoxiconazole", "Bixafen")
dt50_analysis <- vector("list", length = length(dt50_compounds))
dt50_predict <- vector("list", length = length(dt50_compounds))
AS_range <- tibble(compound = dt50_compounds,
                   PP_r = rep(NA, 6),
                   DP_r = rep(NA, 6),
                   xpos = rep(NA, 6),
                   rse_w = rep(NA, 6),
                   rse_s = rep(NA, 6),
                   df_w = rep(NA, 6),
                   df_s = rep(NA, 6))
decay_par <- tibble(decay_w = c(),
                    decay_s = c())
decay <- tibble(decay_w = c(),
                decay_s = c())
for (i in seq_along(dt50_compounds)) {
  dt50_analysis[[i]] <- outlet_all_dates %>%
    left_join(compound_char, by = "compound") %>%
    group_by(compound) %>%
    mutate(conc_W = if_else(conc_W == 0, 0.01, conc_W),
           conc_W = round(conc_W, digits = 0),
           conc_S = round(conc_S, digits = 0),
           rel_w = conc_W / max(conc_W),
           rel_s = conc_S / max(conc_S),
           rel_w = if_else(is.nan(rel_w), 0.01, rel_w),
           rel_w = if_else(rel_w == 0, 0.01, rel_w),
           Kp = conc_S / (conc_W + 0.1)) %>% 
    filter(DBE != 9999) %>%
    group_by(compound) %>%
    arrange(dt50_class) %>%
    filter(compound == dt50_compounds[i])
  
  # fit first order decay for water and sediment
  dt50_predict[[i]] <- expand.grid(DBE = seq(min(dt50_analysis[[i]]$DBE), max(dt50_analysis[[i]]$DBE), 1)) 
  if (dt50_compounds[i] %in% dt50_compounds_w) {
    fit_w <- drm(rel_w ~ DBE, fct = DRC.expoDecay(),
                 data = dt50_analysis[[i]])
    rse_w <- summary(fit_w)$rseMat
    #predict new values for DP
    dt50_predict[[i]]$pred_w <- predict(fit_w, newdata = dt50_predict[[i]])
    #decay[i,1] <- log(2) / fit$coefficients[[2]]
    # DP range labels
    AS_range$DP_r[i] <- str_c(as.character(min(dt50_analysis[[i]]$conc_W)),
                              " - ", as.character(max(dt50_analysis[[i]]$conc_W)))
    AS_range$rse_w[i] <- round(rse_w[1], digits = 2)
    AS_range$df_w[i] <- rse_w[2]
  }
  
  fit_s <- drm(rel_s ~ DBE, fct = DRC.expoDecay(),
               data = dt50_analysis[[i]])
  rse_s <- summary(fit_s)$rseMat
  # predict new values for PP
  dt50_predict[[i]]$pred_s <- predict(fit_s, newdata = dt50_predict[[i]])
  
  dt50_predict[[i]] <- dt50_predict[[i]] %>%
    mutate(compound = dt50_compounds[i])
  
  # PP range labels
  AS_range$PP_r[i] <- str_c(as.character(min(dt50_analysis[[i]]$conc_S)),
                            " - ", as.character(max(dt50_analysis[[i]]$conc_S)))
  AS_range$xpos[i] <- max(dt50_analysis[[i]]$DBE) * 0.6
  AS_range$rse_s[i] <- round(rse_s[1], digits = 2)
  AS_range$df_s[i] <- rse_s[2]
}

rse_range <- AS_range %>%
  select(compound, rse_s, rse_w, df_s, df_w)
dt50_analysis_all <- bind_rows(dt50_analysis)

dt50_predict_all <- bind_rows(dt50_predict)
dt50_r2 <- dt50_analysis_all %>%
  select(compound, DBE, rel_w, rel_s) %>%
  left_join(dt50_predict_all) %>%
  group_by(compound) %>%
  mutate(residuals_w = (rel_w - pred_w)^2,
         squares_w = (rel_w - mean(rel_w, na.rm = T))^2,
         explained_w = (pred_w - mean(rel_w, na.rm = T))^2,
         residuals_s = (rel_s - pred_s)^2,
         squares_s = (rel_s - mean(rel_s, na.rm = T))^2,
         explained_s = (pred_s - mean(rel_s, na.rm = T))^2) %>%
  summarise(SSR_w = sum(residuals_w),
            SST_w = sum(squares_w),
            SSE_w = sum(explained_w),
            SSR_s = sum(residuals_s),
            SST_s = sum(squares_s),
            SSE_s = sum(explained_s),
            r2_w = 1 - SSR_w/SST_w,
            r2_s = 1- SSR_s/SST_s) %>%
  mutate(sum_w = SSR_w + SSE_w,
         sum_s = SSR_s + SSE_s)

dt50_plot_all <- full_join(dt50_analysis_all, dt50_predict_all, by = c("compound", "DBE")) %>%
  left_join(rse_range, by = "compound") %>%
  mutate(rel_w = if_else(compound %in% dt50_compounds_s, NaN, rel_w),
         ypos = if_else(compound == "Bixafen", 0.25, 0.80)) %>%
  group_by(compound) %>%
  mutate(xpos = round(0.6 * max(DBE), digits = 0),
         PP_r = str_c(as.character(min(conc_S, na.rm = T)),
                      " - ", as.character(max(conc_S, na.rm = T))),
         DP_r = str_c(as.character(min(conc_W, na.rm = T)),
                      " - ", as.character(max(conc_W, na.rm = T)))) %>%
  ungroup() %>%
  mutate(xpos = if_else(xpos == DBE, xpos, NaN),
         DP_r = ifelse(is.nan(rel_w), NA, DP_r),
         label = if_else(is.na(DP_r), str_c("PP: ", PP_r, " ppb, rse = ", rse_s), 
                         str_c("PP: ", PP_r, " ppb, rse = ", rse_s, " \n DP: ", DP_r, " ppb, rse = ", rse_w)),
         label = ifelse(is.nan(xpos), NA, label))

dt50_plot_all$compound <- factor(dt50_plot_all$compound, levels = dt50_compounds)

# 3.3 Land use and pesticide transport -----------------------------

# % transport per landuse
source_cereals <- c("Bixafen", "Epoxiconazole", "Isopyrazam")
source_apples <- c("Fluxapyroxad", "Pyraclostrobin", "Chlorantraniliprole", "Difenoconazole") 
source_ap_pot <- c("Glyphosate", "Boscalid") # this two compounds can originate from apples or potato cultivation
# contribution of landuse types to total transport
appl_trans_DP <- sum(filter(AS_transport_rate, compound %in% source_apples)$tot_m_W) / pest_total[1,2] / 1000
appl_trans_PP <- sum(filter(AS_transport_rate, compound %in% source_apples)$tot_m_S) / pest_total[1,1] / 1000
trans_other_DP <- sum(filter(AS_transport_rate, compound %in% source_ap_pot)$tot_m_W) / pest_total[1,2] / 1000
trans_other_PP <- sum(filter(AS_transport_rate, compound %in% source_ap_pot)$tot_m_S) / pest_total[1,1] / 1000
cereal_trans_DP <- sum(filter(AS_transport_rate, compound %in% source_cereals)$tot_m_W) / pest_total[1,2] / 1000
cereal_trans_PP <- sum(filter(AS_transport_rate, compound %in% source_cereals)$tot_m_S) / pest_total[1,1] / 1000

pot_trans_DP <- 1 - appl_trans_DP - cereal_trans_DP - trans_other_DP
pot_trans_PP <- 1 - appl_trans_PP - cereal_trans_PP - trans_other_PP

# occurrence of AS from cereal fields
filter(AS_transport, compound == "Bixafen")$n / 14
filter(AS_transport, compound == "Epoxiconazole")$n / 14
filter(AS_transport, compound == "Isopyrazam")$n / 14

# number of AS associated with potato cultivation
n_AS_potato <- outlet_all_dates %>%
  mutate(pot = if_else(year(date) == 2019 & last_applied == "B", 1, 0),
         pot = if_else(year(date) == 2020 & last_applied == "A", 1, pot)) %>%
  filter(pot == 1) %>%
  distinct(compound)

# Errors and uncertainty -----------------------------------
#' Calculate uncertainty propagation for CH 2.5 and Table S.6

#' Table S.6 - Precipitation
#' in catsop an ARG100 tipping bucket, measures the precipitation. To no the accuracy a test is done.
vol_wat <- 812    # ml water
tip_vol_exp <- 810.4 / 80 #the documented test gives 810.4 ml should result in exactly 80 tips.
tips_exp <- 812 / tip_vol_exp # number of expected tips

cal_tipping <- read_csv("Cal_tipping.dat", skip = 3, locale = nl_tz)
nms <- names(read_csv("Cal_tipping.dat", skip = 1)) %>%
  tolower(.) %>%
  str_replace(., "\\(", "_") %>%
  str_remove(., "\\)")
names(cal_tipping) <- nms

cal_tipping <- cal_tipping %>%
  filter(timestamp > ymd_hms("2021-06-30 11:15:00") & timestamp < ymd_hms("2021-06-30 12:16:00"))

obs_tips <- sum(cal_tipping$rain_tot) / 0.2    # number of tips observed
obs_tip_extra <- 1 / tip_vol_exp    # 1 ml observed in bucket which did not tip
obs_tips_tot <- obs_tips + obs_tip_extra

tip_obs_mm <- signif((tips_exp / obs_tips_tot) * 0.2, digits = 3)

#' Table S.6 - Total suspended sediment
# analyse results of laboratory calibration of sediment inlet boom for Catsop
# load data and calculate other variables
rho_sed <- 2.27
cal <- read_csv("sed_conc_cal_ISCO.csv") %>%
  mutate(vol_source = wat_vol_ml + (sed_source / rho_sed),
         sed_conc_source = sed_source / (vol_source/1000),
         w_sample = f_w_bottle - e_w_bottle,
         sed_sample = f_w_cup - e_w_cup,
         vol_sample = w_sample - sed_sample + (sed_sample / rho_sed),
         sed_conc_sample = sed_sample / (vol_sample/1000),
         rel_source = sed_conc_source / 100,
         rel_sample = sed_conc_sample / 100)

# fit linear model with intercept 0
fit <- lm(cal$rel_sample ~ 0 + cal$rel_source)
SDtss <- summary(fit)$sigma * 100

new <- data.frame(cal$rel_source)
cal <- cal %>%
  mutate(pred = predict.lm(fit, new),
         err2 = (rel_source - pred)^2)
RMSE <- sqrt(sum(cal$err2) / nrow(cal))
#Variance explained (Jin Li, 2017)
#VEcv = (1 - (sum((y_obs_i - y_pred_i)^2))/sum((y_obs_i - y_obs_mean)^2)) * 100

#' Table S.6 - sediment discharge
SD_Qtss <- sqrt((RMSE*100)^2 + 5^2)

#' vwc_corrections
gly_data <- read_csv("data/LC_GLY.csv")
diff_fact_vwc <- tibble(an_name = c("20200627_S_1_5-GLY", "20190315_S_1-GLY", "20190606_S_1_2-GLY",
                                "20200928_S_1_12-GLY", "20200924_S_1_9-GLY", "20190510_S_1_9-GLY"),
                    fact = 2 * 10)
# vwc_gly - correct for vwc
KOH_w <- 0.0395 # g per ml solution.
vwc_gly <- gly_data %>%
  select(-(all_of(str_subset(names(gly_data), "qn|ql|IS"))), - an_code, -number) %>%
  left_join(diff_fact_vwc, by = "an_name") %>%
  mutate(KOH_ml = aimed_w_sample * 5,
         KOH_ml = if_else(!is.na(fact), KOH_ml * 2, KOH_ml),
         samp_w = (d_w_tube - e_w_tube) - (KOH_w * KOH_ml),
         vwc_gly = 1 - (samp_w / aimed_w_sample),
         an_name = str_replace(an_name, "GLY", "LC")) %>%
  select(an_name, vwc_gly) %>%
  filter(!is.na(vwc_gly))
## Check vwc lc-multi samples -----------------------
lc_data <- read_csv("data/LC_multi.csv") %>%
  mutate(ACN_ml = aimed_w_sample * 2,
         an_name = str_replace(an_name, "-B", "-LC")) %>%
  arrange(week, an_code)
vwc_lc <- lc_data %>%
  select(-(all_of(str_subset(names(.), "RT_|Area"))), - number, - an_code, - dil_ratio) %>%
  mutate(samp_w_lc = dry_w - empty_w - aimed_w_sample,
         samp_w_lc = if_else(an_name == "20190819_C_L-LC", samp_w_lc - 2.5, samp_w_lc),
         vwc_lc = 1 - (samp_w_lc / aimed_w_sample)) %>% #adjust for double salt weight
  left_join(vwc_gly, by = "an_name") %>%
  mutate(salt_eff = (vwc_lc - vwc_gly) / aimed_w_sample,
         diff_vwc = (vwc_lc - vwc_gly),
         diff_w_lc = diff_vwc * aimed_w_sample) %>%
  filter(!is.na(vwc_lc)) %>%
  group_by(aimed_w_sample) %>%
  mutate(salt_eff = mean(salt_eff))

fit1 <- lm(vwc_lc$vwc_gly ~ vwc_lc$vwc_lc)
vwc_rse <- summary(fit1)$sigma * 100

new <- data.frame(vwc_lc$vwc_lc)
vwc_lc <- vwc_lc %>%
  ungroup() %>%
  mutate(pred = predict(fit1, new),
         err2 = (vwc_gly - pred)^2)
RMSE_vwc <- sqrt(sum(vwc_lc$err2, na.rm = T) / nrow(vwc_lc))

#' Table S.6 - Pesticide concentrations
recov_stats <- read_csv("recovery_overview.csv")
Pest_sd <- max(recov_stats$sd) * 100

#' Table S.6 - Pesticide load water
SD_pest_w <- sqrt(Pest_sd^2 + 5^2)

#' Table S.6 - Pesticide load TSS
SD_pest_tss <- sqrt(Pest_sd^2 + SD_Qtss^2 + vwc_rse^2)

# Supplementary materials ------------------
## Fig S.2 calculate ------------
# analyse the differences in runoff and erosion between the analysed and 
# not analysed events.
erosion25 <- erosion25 %>%
  select(date, Sed_kg, Sed_conc_max, mean_tss)
runoff_analysis <- left_join(runoff_analysis, erosion25, by = "date")
events14 <- total %>%
  select(date, Sed_kg, Sed_conc_max, mean_tss, D, Q_m3, Q_max) %>%
  mutate(group = "Analyzed")
not_analysed <- anti_join(runoff_analysis, events14, by = "date") %>%
  mutate(group = "Not Analyzed")%>%
  select(date, group, Sed_kg, Sed_conc_max, mean_tss, D, Q_m3, Q_max)

analysed_w <- anti_join(runoff_analysis, not_analysed, by = "date") %>%
  mutate(group = "Analyzed")
#' events on 2019-08-18 is removed because the water level readings stay
#' high for the rest of the day, this is an extreme outlier most likely due
#' to sensor failure or obstruction.
df_runoff_w <- bind_rows(analysed_w, not_analysed) %>%
  filter(date != "2019-08-18") %>%
  select(group, D, Q_m3, Q_max) %>%
  pivot_longer(cols = (D:Q_max), names_to = "var", values_to = "value") %>%
  mutate(value = if_else(value == -Inf, 0, value))

df_runoff_s <- bind_rows(events14, not_analysed) %>%
  filter(date != "2019-08-18") %>%
  select(group, Sed_kg, Sed_conc_max, mean_tss) %>%
  pivot_longer(cols = (Sed_kg:mean_tss), names_to = "var", values_to = "value") %>%
  mutate(value = if_else(value == -Inf, 0, value)) %>%
  filter(!is.na(value))

df_runoff <- bind_rows(df_runoff_s, df_runoff_w)
df_runoff$var <- factor(df_runoff$var, levels = c("Q_m3", "Q_max", "D",
                                                  "Sed_kg", "Sed_conc_max", "mean_tss"),
                          labels = c("Total discharge (m3)",
                                     "Maximum discharge (m3/sec)",
                                     "Duration (min)",
                                     "Total erosion (kg)",
                                     "Maximum TSS (g/L)",
                                     "Mean TSS (g/L)"))
                  
# t-test difference between analysed and not analysed
t_test_w <- df_runoff_w %>%
  group_by(var) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "group")
t_test_s <- df_runoff_s %>%
  group_by(var) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "group")

# compare means between total erosion in analysed and not analysed events
mean(filter(df_runoff_s, var == "Sed_kg" & group == "Not Analyzed")$value, na.rm = T)
mean(filter(df_runoff_s, var == "Sed_kg" & group == "Analyzed")$value, na.rm = T)

## Point source calculations -------------------
pnt_scr_sed <- pnt_src_lc %>%
  filter(samp_type == "S")%>%
  group_by(date, samp_type) %>%
  summarise_at(vars(starts_with('conc_')), ~weighted.mean(., sed_gr, na.rm = T))
pnt_src_conc <- pnt_src_lc %>%
  select(-sed_gr) %>%
  filter(samp_type == "W") %>%
  group_by(date, samp_type) %>%
  summarise_at(vars(starts_with('conc_')), ~mean(., na.rm = T)) %>%
  bind_rows(pnt_scr_sed) %>%
  mutate(date = date(date)) %>%
  pivot_wider(names_from = samp_type, values_from = (conc_Glyphosate:conc_AMPA)) %>%
  left_join(total, by = "date") %>%
  select(-(TP_sed:TP_wat_loq)) %>%
  filter(date != "2020-09-26") %>%
  mutate(TP_sed_GLY = conc_Glyphosate_S * Sed_kg,
         TP_sed_AMPA = conc_AMPA_S * Sed_kg,
         TP_wat_GLY = conc_Glyphosate_W * Q_m3 * 1000,
         TP_wat_AMPA = conc_AMPA_W * Q_m3 * 1000,
         PP_rat_GLY = TP_sed_GLY / (TP_sed_GLY + TP_wat_GLY),
         PP_rat_AMPA = TP_sed_AMPA / (TP_sed_AMPA + TP_wat_AMPA),
         fact_phase_GLY = conc_Glyphosate_S / conc_Glyphosate_W,
         fact_phase_AMPA = conc_AMPA_S / conc_AMPA_W)

pnt_table <- pnt_src_conc %>%
  select((date:conc_AMPA_S), starts_with("TP_"))
write_csv(pnt_table, "pnt_table.csv")

## Table S.1 ---------------------------
# Annex A - list of applied and analysed pesticides

all_applied <- pest_desc %>%
  semi_join(pest_appl, by = "pesticide") %>%
  select(compound) %>%
  distinct() %>%
  left_join(compound_char, by = "compound") %>%
  mutate(analysed = if_else(compound %in% analysed, "Yes", "No"),
         kfd = as.numeric(str_remove(kd_note, "f")),
         kfoc = as.numeric(str_remove(koc_note, "f"))) %>%
  select(compound, dt50, kfd, kfoc, Sw, type, analysed) %>%
  rename(AS = compound) %>%
  filter(type %in% c("herbicide", "fungicide", "insecticide"))

#write_csv(all_applied, "compound_overview.csv")
## Table S.7 ------------------------------

#' Table S.7 - number of application days per year for each field
appl_days_pest <- left_join(pest_appl, pest_desc, by = "pesticide") %>%
  left_join(compound_char, by = "compound") %>%
  filter(!is.na(type)) %>%
  filter(type != "plant growth regulator") %>%
  group_by(large, date) %>%
  distinct(date, large) %>%
  mutate(year = year(date)) %>%
  group_by(large, year) %>%
  summarise(days = n())

#' Table S.7 - kg/ha AI per crop type per year.
appl_rate <- left_join(pest_appl, pest_desc, by = "pesticide") %>%
  left_join(compound_char, by = "compound") %>%
  filter(!is.na(type)) %>%
  filter(type != "plant growth regulator") %>%
  select(-area) %>%
  mutate(ai_appl_dose = dose * ai_dose,
         year = year(date)) %>%
  left_join(large_area, by = "field_nr") %>%
  group_by(field_nr, large, year) %>%
  summarise(kg_ai_ha_f = sum(ai_appl_dose, na.rm = T) / 1000,
            area = as.numeric(mean(area))) %>%
  group_by(large, year) %>%
  summarise(kg_ai_ha = weighted.mean(kg_ai_ha_f, area, na.rm = T))

#' Table S.7 - number of different AI applied
number_ai <- left_join(pest_appl, pest_desc, by = "pesticide") %>%
  left_join(compound_char, by = "compound") %>%
  filter(!is.na(type)) %>%
  filter(type != "plant growth regulator") %>%
  mutate(year = year(date)) %>%
  group_by(large, year) %>%
  distinct(large, compound) %>%
  summarise(n = n())

#' Table S.7 - number of AI analysed in this study
ai_analysed_number <- left_join(pest_appl, pest_desc, by = "pesticide") %>%
  left_join(compound_char, by = "compound") %>%
  filter(!is.na(type)) %>%
  filter(type != "plant growth regulator") %>%
  mutate(year = year(date)) %>%
  group_by(large, year) %>%
  distinct(large, compound) %>%
  left_join(ai_analysed, by = "compound") %>%
  summarise(analysed = sum(analysed, na.rm = T))

#' Table S.7 - number of ai detected on per field per year
ai_detected <- pest_catch %>%
  mutate(year = year(date)) %>%
  filter(!is.na(large) & conc > 0,
         compound != "Tebuconazole") %>%
  group_by(year, large) %>%
  distinct(compound) %>%
  summarise(detected = n())

#' Table S.7 - number of samples collected per field per year
sample_large <- soil_samp %>%
  left_join(fields_large, by = "field_nr") %>%
  mutate(year = year(date)) %>%
  filter(!is.na(large)) %>%
  group_by(year, large) %>%
  summarise(samples = n())

# Visualization -----------------------------------------------------------
#' List of figures in paper:
#' fig 1: map of catchment - QGIS
#' fig 2: overview of rain, runoff, erosion and DP & PP transport per event
#' fig 3: concentrations at outlet per event
#' fig 4: decrease in concentration over time after application for 6 AS

## graph settings ----------------------
coeff <- 0.25
Q_scale_min <- 200
P_scale_min <- 30
margin <- c(0,0)

# Define different color pallets for figures
# 2 colors for water and sediment OR DP and PP
colors1 <- c("#E6BB00", "#0072B2")

colors2 <- c("#E6BB00", "#0072B2", "#CC0000", "#8df55b")

# 12 color version
colors3 <- c("#7465AC", "#E6BB00", "#56B4E9",
          "#50ECCC", "#009E73", "#F0E442", "#885522",
          "#0072B2", "#CC0000", "#8df55b", "#BFC816", "#999999")

my_theme = theme(
  text = element_text(size = 12, family = "Times New Roman"),
  axis.title.x = element_text(size = 12, family = "Times New Roman"),
  axis.text.x = element_text(size = 12, family = "Times New Roman"),
  axis.title.y = element_text(size = 12, family = "Times New Roman"),
  plot.margin = margin(6, 6, 6, 0)
)

## figure 2 -----------
#' fig 2: overview runoff and precipitation events
Q_scale <- 10000
P_scale <- 30

Q_plot <- ggplot(total)+
  geom_bar(aes(x = number, y = Q_m3), stat = "identity", width = 0.5, fill = colors1[2]) +
  theme_classic() +
  #geom_vline(xintercept = 4.5, linetype = "dashed") +
  labs(x = "", y = bquote(Runoff~(m^"3")), parse = T) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(position = "bottom", breaks = NULL) +
  my_theme

#Q_plot
Sed_plot <- ggplot(total)+
  geom_bar(aes(x = number, y = Sed_kg), stat = "identity", width = 0.5, fill = colors1[1]) +
  theme_classic() +
  #geom_vline(xintercept = 4.5, linetype = "dashed") +
  labs(x = "", y = "Erosion (kg)", parse = T) +
  scale_x_continuous(position = "bottom", breaks = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1)) +
  my_theme

#Sed_plot
rain_plot <- ggplot(total)+
  geom_bar(data = total, aes(x = number, y = P), stat = "identity", width = 0.5, alpha = 0.5) +
  geom_point(data = total, aes(x = number, y = evi)) +
  scale_y_reverse(limits = c(P_scale, 0), expand = margin) +
  #geom_vline(xintercept = 4.5, linetype = "dashed") +
  theme_classic() +
  scale_x_continuous(position = "top", breaks = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  labs(x = "", y = "P (mm) & \n EVI [-]", title = "", parse = T) +
  my_theme
#rain_plot

pest_plot <- ggplot(pest_tot) +
  geom_bar(aes(x = number, y = P, fill = samp_type), stat = "identity", position = "dodge", width = 0.5) +
  geom_errorbar(aes(x = number, y = P, ymin = P - error, ymax = P + error, fill = samp_type),
                position = position_dodge(width = 0.5), width = 0.5) +
  theme_classic() +
  labs(x = "", y = "Discharged \n load (g)", title = "",
       parse = T) +
  scale_fill_manual(name = "", labels = c("PP", "DP"), values = colors1) +
  scale_x_continuous(breaks = 1:14, labels = total$date_lab) +
  theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1)) +
  my_theme +
  scale_y_continuous(expand = c(0,0)) 
#pest_plot

# pest_plot without a legend
pest_plot_ng <- pest_plot + 
  theme(legend.position = "none")

legend_5 <- get_legend(
  pest_plot + theme(legend.text.align = 0,
                    legend.box.margin = margin(0,0,0,0),
                    legend.position = c(0.5, 0.27)))

plot_pest <- plot_grid(rain_plot, Q_plot, Sed_plot, pest_plot_ng, 
                       nrow = 4, align = "v", 
                       axis = "rl", rel_heights = c(4,4,4,7), labels = c("A", "B", "C", "D"),
                       label_fontfamily = "Times New Roman", hjust = -1.5)
#plot_pest
plot_pest2 <- plot_grid(plot_pest, legend_5, ncol = 2, rel_widths = c(2, .2))
#plot_pest2

ggsave("images/figure2.tiff", width = 180, height = 150, units = "mm",
       dpi = 450, device = "tiff")
## figure 3 ------------
#' fig 3: compound specific transport at outlet conc

#colors for fig 3
colors3_s <- c("#7465AC", "#E6BB00", "#56B4E9",
             "#50ECCC", "#009E73", "#F0E442", "#885522",
             "#0072B2", "#CC0000", "#8df55b", "#999999")

colors3_w <- c("#7465AC", "#56B4E9",
               "#F0E442", "#0072B2", "#CC0000", "#8df55b", "#999999")
# figure PP
event_bar_s <- ggplot()+
  geom_bar(data = filter(comp_box, samp_type == "S"), 
           aes(x = date, y = conc, fill = other_s),
           position = "stack", stat = "identity",
           width = 0.5) +
  theme_classic() +
  coord_cartesian(ylim = c(0, oob_fig3_s)) +
  labs(x = element_blank(), y = bquote(mg~kg^"-1"), 
       title = element_blank(),
       parse = T) +
  geom_text(data = comp_box, aes(x = date, y = oob_fig3_s, label = label_oob_s),
            angle = 90, family = "Times New Roman") +
  scale_y_continuous(expand = expansion(c(0, .1))) +
  scale_x_discrete(labels = element_blank()) +
  scale_fill_manual(values = colors3_s, name = "Active\nSubstance",
                    labels = function(x) str2expression(paste0(x))) +
  my_theme
#event_bar_s
#save plot without legend
event_bar_sp <- event_bar_s + 
  theme(legend.position = "none")
# figure DP
event_bar_w <- ggplot()+
  geom_bar(data = filter(comp_box, samp_type == "W"), 
           aes(x = date, y = conc, fill = other_w),
           position = "stack", stat = "identity",
           width = 0.5) +
  theme_classic() +
  coord_cartesian(ylim = c(0, oob_fig3_w)) +
  scale_x_discrete(labels = total$date_lab) +
  theme(legend.text.align = 0,
        axis.text.x = element_text(angle = 45, hjust=1,vjust=1)) +
  labs(x = element_blank(), y = bquote(mg~L^"-1"), 
       title = element_blank(),
       parse = T) +
  geom_text(data = comp_box, aes(x = date, y = oob_fig3_w, label = label_oob_w),
            angle = 90, family = "Times New Roman") +
  scale_y_continuous(expand = expansion(c(0, .1))) +
  scale_fill_manual(values = colors3_w, name = "Active\nSubstance",
                    labels = function(x) str2expression(paste0(x))) +
  guides(fill = "none") +
  my_theme
#event_bar_w

# get legend form plot and add to plot grid.
legend_4 <- get_legend(
  event_bar_s + theme(legend.text.align = 0,
                      legend.box.margin = margin(0,0,0,0)))

plot <- plot_grid(event_bar_sp, event_bar_w, nrow = 2, align = "v", axis = "rl",
                  labels = c("PP", "DP"), label_fontfamily = "Times New Roman", 
                  label_size = 12, rel_heights = c(4,4))
plot2 <- plot_grid(plot, legend_4, ncol = 2, rel_widths = c(2, .8))

ggsave("images/figure3.tiff", width = 180, height = 120, units = "mm", dpi = 300, device = "tiff")

## figure 4 -------------------------
# PP + DP
vars <- c("PP" = "#E6BB00", "DP" = "#0072B2")
dt50_plot <- ggplot(dt50_plot_all) +
  geom_point(aes(x = DBE, y = rel_s, color = "PP")) +
  geom_line( aes(x = DBE, y = pred_s, color = "PP")) +
  geom_point(aes(x = DBE, y = rel_w, color = "DP")) +
  geom_line(aes(x = DBE, y = pred_w, color = "DP")) +
  scale_colour_manual(name = "", values = vars) +
  facet_wrap(~compound, scales = "free_x", nrow = 2, strip.position = "top") +
  ylim(c(0,1)) +
  expand_limits(x = 0) +
  geom_text(aes(x = xpos, y = ypos, label = label), family = "Times New Roman", size = 2.5) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside",
        panel.grid = element_blank(), plot.margin = margin(10, 10, 10, 10)) +
  labs(x = "Days after application", y = "Relative concentration") +
  my_theme
#dt50_plot  
ggsave("images/figure4.tiff", width = 180, height = 120, units = "mm", dpi = 300, device = "tiff")

## figure S.2 -------------------
# figure to compare runoff and erosion of 14 events with all 39 occurred.
runoff_plot <- ggplot(df_runoff) +
  geom_boxplot(aes(x = var, y = value, fill = group, alpha = group)) +
  facet_wrap(~var, scale = "free", strip.position = "top") +
  theme_bw() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("#0072B2", "#0072B2")) +
  scale_alpha_manual(values = c(1,0.3)) +
  scale_x_discrete(position = "top") +
  theme(strip.background = element_blank(), strip.placement = "outside",
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank(), legend.position = "bottom",
        legend.title = element_blank()) +
  my_theme + 
  theme(plot.margin = margin(5, 5, 5, 5))
runoff_plot

ggsave("images/figureS2.tiff", width = 180, height = 180, units = "mm", dpi = 300, device = "tiff")
