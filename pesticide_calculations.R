# Load and analyze the raw data LC-MS/MS results exported from MassLynx

# Initialization ------------------
library(tidyverse)

# GLY data ----------------------------------------------------------------
### load data ----------------
gly_data <- read_csv("data/LC_GLY.csv") %>%
  mutate(KOH_ml = aimed_w_sample * 5)
# correction for dilutions in preparation process.
# first tibble is for 6 samples with 1 gram sediment but 10 ml KOH, instead of 5 mL.
diff_fact <- tibble(an_name = c("20200627_S_1_5-GLY", "20190315_S_1-GLY", "20190606_S_1_2-GLY",
                                "20200928_S_1_12-GLY", "20200924_S_1_9-GLY", "20190510_S_1_9-GLY"),
                    fact = 2 * 10)
# g/ml (salt weight) of 0.6M KOH solution (purity = 85%)
# prepared 39.5294 g in 1000 ml
KOH_w <- 0.0395 # g per ml solution.
diff_fact_vwc <- diff_fact
### loop over compounds --------------------
#' calculations; loop over each batch and compound to calculate concentrations and evaluation criteria.
#' concentration is calculated with a linear model fit on the combined calibration curve data.
compounds <- c("GLY", "AMPA")
batch <- c(49, 50)
df_list <- vector("list", length = length(compounds))
df_recov <- vector("list", length = length(compounds))
batch_list <- vector("list", length = length(batch))
batch_recov <- vector("list", length = length(batch))
models <- vector("list", length = length(batch)*length(compounds))
models_stats <- tibble(name = c(), batch = c(), intercept = c(), b = c(), r_sqrd = c())
nms_basic <- c(names(gly_data)[1:6], names(gly_data)[19:24])
for (j in seq_along(compounds)) {
nms <- str_subset(names(gly_data), compounds[j])
nms <- c(nms, nms_basic)
new_nms <- str_replace(nms, compounds[j], "")
temp_c <- gly_data %>%
  select(all_of(nms)) %>%
  filter(!str_detect(an_name, "^Blank"))
names(temp_c) <- new_nms
temp_c <- temp_c %>%
  mutate(RT_m = if_else(an_name == "Blank", 0, (RT__qn + RT__ql +RT__IS)/3),
         ir = if_else(an_name == "Blank", 0, Area__ql / Area__qn),
         resp = if_else(an_name == "Blank", 0, Area__qn / Area__IS))
results <- temp_c %>%
  filter(!is.na(samp_type)) %>%
  select(an_name, pest_ID, samp_type)
recov <- temp_c %>%
  filter(an_type == "QC") %>%
  select(an_name, week)
### loop over batch -----------------------------------
#' calculate calibration curve
for (i in seq_along(batch)) {
  #### fit model ---------------------------
  cal <- temp_c %>%
  filter(an_type == "curve") %>%
  filter(week == batch[i]) %>%
  select(an_name, resp, ir) %>%
  mutate(cal_conc = str_extract(an_name, "(\\d\\.\\d*)|(\\d)"),
         cal_conc = as.numeric(if_else(is.na(cal_conc), "0", cal_conc)))
#fit model on calibration curve
x <- (j-1)*length(batch)+i
models[[x]] <- lm(resp ~ cal_conc, data = cal)
coeff <- models[[x]]$coefficients
models[[x]]$name <- str_c(compounds[j], "_", batch[i])
model_x_stats <- tibble(name = str_c(compounds[j]),
                        batch = str_c(batch[i]),
                        intercept = coeff[1],
                        b = coeff[2],
                        r_sqrd = summary(models[[x]])$r.squared)
models_stats <- bind_rows(models_stats, model_x_stats)
#' calculate final concentrations, and adjust for water in initial samples
#' the sediment samples are often wet, so after analysis they are dried, and the weight is corrected for
#' the KOH salts. The difference between aimed weight (2gr) and real weight is water. This does also contain
#' pesticides, so the water concentration is used to adjust for that.
temp <- temp_c %>%
  filter(week == batch[i]) %>%
  left_join(diff_fact, by = "an_name") %>%
  mutate(detec_conc = (resp - coeff[1]) / coeff[2],
         samp_conc = if_else(an_type == "curve" | 
                                an_type == "point" | 
                                samp_type == "W", detec_conc * 2, detec_conc * 10),
         samp_conc = if_else(an_type == "QC", detec_conc * 10, samp_conc),
         samp_conc = if_else(is.na(fact), samp_conc, detec_conc * fact),
         LOQ = if_else(detec_conc < 0.005, "<LOQ", ">"),
         KOH_ml = if_else(is.na(fact), KOH_ml, 10),
         samp_w = (d_w_tube - e_w_tube) - (KOH_w * KOH_ml))
temp_ir <- cal %>%
  filter(str_detect(an_name, "Blank", negate = T))
m_ir <- mean(temp_ir$ir)
min_ir <- m_ir * 0.7
max_ir <- m_ir * 1.3

temp <- temp %>%
  mutate(IR_error = if_else(ir > max_ir | ir < min_ir, 1, 0))
### calculate recovery -----------------------
#' recovery for calibration curve is compared with the detection conc. For the QC it is with the 
#' concentration in the sample tubes (not corrected for water). The QC samples are corrected for already
#' available concentration in the mother material.
#' forgot to add IS to QC1-0.2, week 50 which gives strange values for both GLY and AMPA (see also peaks in MassLynx) -> remove
origin_s <- temp %>%
  select(an_name, samp_conc) %>%
  rename(origin = an_name, samp_conc_o = samp_conc) %>%
  mutate(samp_conc_o = if_else(samp_conc_o < 0, 0, samp_conc_o))
QC_sample <- temp %>%
  filter(an_type == "QC") %>%
  mutate(origin = if_else(an_name == "XY-0.5 Spiked", "XY-BLANK", str_c(str_extract(notes, "\\d\\d.*$"), "-GLY"))) %>%
  left_join(origin_s, by = "origin") %>%
  select(an_name, origin, samp_conc_o) %>%
  mutate(samp_conc_o = if_else(str_detect(origin, "20200304"), 0.0123, samp_conc_o))
#'calculate recovery based on the real concentration: the concentration of the sample before the preparation
#' process.
recov_temp <- temp %>%
  filter(an_type == "curve" | an_type == "QC") %>%
  left_join(QC_sample, by = "an_name") %>%
  mutate(real_conc = if_else(an_type == "curve" | an_type == "point",
                             as.numeric(str_extract(an_name, "(\\d+\\.\\d*)|(\\d+)")) * 2, 
                             as.numeric(str_extract(an_name, "(\\d\\.\\d*$)|(\\d$)"))),
         real_conc = if_else(is.na(real_conc), 0, real_conc),
         real_conc = if_else(an_name == "XY-0.5 Spiked", 0.5, real_conc),
         real_conc = if_else(str_detect(an_name, "RS"), real_conc / 2, real_conc),
         recov = if_else(an_type == "curve" | an_type == "point", samp_conc / real_conc, 
                         (samp_conc - samp_conc_o) / real_conc),
         orig_eff = samp_conc_o / real_conc,
         recov = if_else(orig_eff > 0.3, NaN, recov)) # when original conc is >30% QC is not valid
batch_recov[[i]] <- recov_temp %>%
  filter(an_type == "QC") %>%
  mutate(recov = if_else(IR_error == 1, NaN, recov)) %>%
  select(recov, orig_eff) %>%
  rename_with(.cols = everything(), ~ str_c(. , "_", compounds[j]))


## VWC adjustment --------------------------------
water_s <- temp %>%
  filter(samp_type == "W") %>%
  select(pest_ID, samp_conc, LOQ) %>%
  rename(samp_conc_w = samp_conc,
         w_LOQ = LOQ)

temp <- temp %>%
  left_join(water_s, by = "pest_ID") %>%
  mutate(cor_w_part = if_else(w_LOQ == "<LOQ"| is.na(w_LOQ), 0, 
                              (((aimed_w_sample - samp_w) / 0.998) * samp_conc_w)/ aimed_w_sample),
         final_conc = signif(if_else(is.na(aimed_w_sample), samp_conc, 
                                     ((aimed_w_sample / samp_w) * samp_conc)-cor_w_part),
                             digits = 3))
         
#' prepare one column with final concentrations, give -999 for all samples that do not pass evaluation criteria
batch_list[[i]] = temp %>%
  filter(!is.na(samp_type)) %>%
  mutate(final_conc = signif(if_else(LOQ == "<LOQ" | IR_error == 1, -999, final_conc), digits = 3)) %>%
  select(final_conc) # | IR_error == 1
names(batch_list[[i]]) <- str_c("conc_", compounds[j])
}
#' combine batches
df_list[[j]] <- bind_rows(batch_list)
df_recov[[j]] <- bind_rows(batch_recov)
}
### combine all data to 1 overview ---------------------
#' combine compounds and add sample names etc.
pest_result <- bind_cols(df_list)
names(pest_result) <- c("conc_Glyphosate", "conc_AMPA")
pest_recovery <- bind_cols(recov, df_recov) %>%
  filter(!(an_name == "GLY-QC1 - 0.2" & week == 50))
pest_result <- bind_cols(results, pest_result) %>%
  arrange(pest_ID)
# calculate recovery statistics
recov_temp <- pest_recovery %>%
  select(all_of(starts_with("recov")))
recov_stats <- tibble(compound = compounds,
                     mean = round(sapply(recov_temp, mean, na.rm = T), digits = 1),
                     sd = round(sapply(recov_temp, sd, na.rm = T), digits = 2),
                     n = sapply(recov_temp, function(x){sum(!is.na(x))}))
# change unit to ng/ml so it is the same as the LC_multi results
# adjust for recovery
d1 <- as.matrix(pest_result[4:5])
v1 <- as.matrix(recov_stats[2])
lc_gly <- as_tibble(sweep(d1, 2, v1, FUN="/")) %>%
  bind_cols(pest_result[1:3]) %>%
  mutate(conc_Glyphosate = if_else(conc_Glyphosate > 0, conc_Glyphosate * 1000, conc_Glyphosate),
         conc_AMPA = if_else(conc_AMPA > 0, conc_AMPA * 1000, conc_AMPA)) %>%
  select(-an_name) %>%
  mutate(across(starts_with("conc_"), ~ if_else(. < -500, -999, .)))
recov_gly <- recov_stats
model_gly <- models_stats

# LC-multi ----------------------------------------------------------------

### load data --------------------
lc_data <- read_csv("data/LC_multi.csv") %>%
  mutate(ACN_ml = aimed_w_sample * 2,
         an_name = str_replace(an_name, "-B", "-LC")) %>%
  arrange(week, an_code)
# batch number of weeks with diluted samples form earlier batches
dil_wks <- c("12", "18")
#' Drop 2nd 'B4 Sta 1.25' from batch 12 - this has a bad overall response, and can not be included.
lc_data <- lc_data %>%
  filter(an_code != "TQS3_210326_075")
#' add old batch number for dilutions batch, to be able to analyse in two groups.
old_batch <- lc_data %>%
  filter(str_detect(an_name, "B3"))
old_batch_n12 <- max(filter(old_batch, week == "12")$number)
old_batch_n18 <- max(filter(old_batch, week == "18")$number)
old_batch <- lc_data %>%
  filter(week %in% dil_wks) %>%
  mutate(old_batch = if_else(week == "12" & number <= old_batch_n12, "03", "04"),
         old_batch = if_else(week == "18" & number <= old_batch_n18, "03", old_batch))

lc_data <- lc_data %>%
  filter(!(week %in% dil_wks)) %>%
  bind_rows(old_batch) %>%
  arrange(week, an_code)
#' add a 'batch' variable, with week and specific old batch numbers for the dilution batch; 3, 4, 123, 124.
lc_data <- lc_data %>%
  mutate(batch = if_else(week %in% dil_wks, str_c(week, old_batch), week))

#' correction for dilutions in preparation process.
#' first tibble is for samples with different sediment to ACN ratios.
diff_nms <- tibble(an_name = c("20190819_C_L-LC", "20190315_S_1_2-LC", "20190606_S_1_2-LC", "20200924_S_1_9-LC"),
                   ACN_ml = c(10, 5, 2.5, 2.5),
                   aimed_w_sample = c(2.5, 1.7, 1.35, 1.25),
                   pest_ID = c(10081, 10001, 10011, 10044))
diff_fact <- lc_data %>%
  select(-ACN_ml, -aimed_w_sample, -pest_ID) %>%
  semi_join(diff_nms, by = "an_name") %>%
  left_join(diff_nms, by = "an_name")
lc_data <- anti_join(lc_data, diff_fact, by = "an_name") %>%
  bind_rows(diff_fact) %>%
  mutate(samp_type = if_else(str_detect(an_name, "(Std |Sta )"), "std", samp_type),
         samp_type = if_else(an_type == "QC", "C", samp_type),
         samp_type = if_else(str_detect(an_name, "W-QC1"), "W", samp_type),
         an_name = str_replace(an_name, "-B", "-LC"),
         aimed_w_sample = if_else(samp_type == "C" & an_type == "QC", 5, aimed_w_sample),
         ACN_ml = aimed_w_sample * 2) %>%
  arrange(week, an_code)

# vwc_gly
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

vwc_rse <- summary(lm(vwc_lc$vwc_lc ~ vwc_lc$vwc_gly))$sigma

#' The salts in the lc-multi have attracted some water after drying, this causes negative vwc values
#' the solution for this problem is to use the vwc values of the gly samples.

lc_data <- lc_data %>%
  left_join(vwc_gly, by = "an_name")

### calculate TSS for 2020-08-12 --------------------
#' this event only had 1 bottle per sample moment. So either pesticides or TSS could be analyzed
#' We decided to analyze pesticides, and afterwards estimate TSS from the pesticide samples.
tss_est_gly <- gly_data %>%
  select(pest_ID, aimed_w_sample) %>%
  rename(aimed_w_gly = aimed_w_sample)

cent_tube_w <- read_csv("sources/cent_tube_w.csv") %>%
  summarise(e_w = mean(e_w))
cent_tube_w <- cent_tube_w$e_w
est_vol_samp <- 645 # based on the mean + sd of other samples in 2020.
tss_est <- lc_data %>%
  filter(!is.na(B_w)) %>%
  select(an_name, B_w, vwc_gly, pest_ID, aimed_w_sample) %>%
  rename(aimed_w_lc = aimed_w_sample) %>%
  left_join(tss_est_gly, by = "pest_ID") %>%
  filter(!is.na(aimed_w_gly)) %>%
  mutate(B_w = B_w - cent_tube_w, 
         est_mass = (aimed_w_gly + aimed_w_lc + B_w) * (1 - vwc_gly),
         n_bot = c(2,1,1,1,1,1,3),
         vol_wat = n_bot * est_vol_samp,
         vol_sed = est_mass / 2.65,
         sed_conc = (est_mass * 1000) / (vol_sed + vol_wat),
         sed_gr = est_mass / n_bot) %>%
  select(pest_ID, sed_conc, sed_gr)
write_csv(tss_est, "tss_20200812.csv")

### combine A + B peaks -----------------
#' combine peak areas of A and B compounds and update lc_data table
compounds <- unique(str_replace(unique(str_replace(str_subset(names(lc_data), "(_A|_B)$"), "(_A|_B)$", "")), "RT_|Area_|Area1_", ""))
nms <- str_replace(str_subset(names(lc_data), compounds[1]), compounds[1], "")
df_com <- vector("list", length = length(compounds))
for (i in seq_along(compounds)) {
combin <- lc_data[str_subset(names(lc_data), compounds[i])]
names(combin) <- nms
combin <- combin %>%
  replace(is.na(.), 0) %>%
  mutate(RT_ = (RT__A + RT__B) / 2,
         Area_ = (Area__A + Area__B),
         Area1_ = (Area1__A + Area1__B)) %>%
  select(-all_of(nms))
new_nms <- str_c(names(combin), compounds[i])
names(combin) <- new_nms
df_com[[i]] <- combin
}
com <- bind_cols(df_com)
#' drop all columns with A and B data
drop_nms <- str_subset(names(lc_data), "(_A|_B)$")
lc_data <- lc_data %>%
  select(-all_of(drop_nms)) %>%
  bind_cols(com)

#' drop compounds that have a bad response and will not be included in the analysis.
#' Or are not relevant for further analysis
drop_nms <- str_subset(names(lc_data), "Mesosulfuron-methyl|Pendimethalin|Trinexapac|Tebuconazole")
lc_data <- lc_data %>%
  select(-all_of(drop_nms))

#' calculations; loop over each batch and compound to calculate concentrations and evaluation criteria.
#' concentration is calculated with a linear model fit on the combined calibration curve data.
compounds <- str_replace(str_subset(names(lc_data), "RT_"), "RT_", "")
#' remove Caffeine from calculations
compounds <- compounds[-1]

btch <- unique(lc_data$batch)
df_list <- vector("list", length = length(compounds))
df_qual <- vector("list", length = length(compounds))
df_recov <- vector("list", length = length(compounds))
batch_list <- vector("list", length = length(btch))
qual_list <- vector("list", length = length(btch))
batch_recov <- vector("list", length = length(btch))
models <- vector("list", length = length(btch)*length(compounds))
models_stats <- tibble(name = c(), batch = c(), intercept = c(), b = c(), r_sqrd = c())
curve_max <- vector("list", length = length(btch)*length(compounds))
curve_min <- vector("list", length = length(btch)*length(compounds))
limit <- tibble(compound = c(),
                batch = c(),
                loq = c(),
                limit = c())
nms_basic <- subset(names(lc_data), !str_detect(names(lc_data), "(RT|Ar)"))
### loop over compounds -----------------------
for (j in seq_along(compounds)) {
  nms <- str_subset(names(lc_data), compounds[j])
  nms <- c(nms, nms_basic)
  new_nms <- str_replace(nms, compounds[j], "")
  temp_c <- lc_data %>%
    select(all_of(nms)) %>%
    filter(!str_detect(an_name, "^Blank"))
  names(temp_c) <- new_nms
  temp_c <- temp_c %>%
    mutate(Area_ = if_else(is.na(Area_), 0, Area_),
           Area1_ = if_else(is.na(Area1_), 0, Area1_),
           ir = Area_ / Area1_,
           ir = if_else(is.infinite(ir), NaN, ir))
  results <- temp_c %>%
    filter(an_type == "sample") %>%
    select(an_name, pest_ID, samp_type, week)
  recov <- temp_c %>%
    filter(an_type == "QC") %>%
    select(an_name, batch)
  ### loop over batch --------------------
  for (i in seq_along(btch)) {
    #### fit model -------------------------
    cal <- temp_c %>%
      filter(an_type == "curve") %>%
      filter(batch == btch[i]) %>%
      select(an_name, Area_, ir) %>%
      mutate(cal_conc = as.numeric(str_extract(an_name, "( \\d+\\.\\d*)|( \\d+)"))) %>%
      group_by(cal_conc) %>%
      mutate(mean = mean(Area_)) %>%
      ungroup()
    pnt <- 3.125
    pnt_mean <- filter(cal, cal_conc == pnt)
    pnt_mean <- mean(pnt_mean$mean)
    #' only calculate model etc if data is available (in batch 12 several compounds are not measured)
    if (pnt_mean > 0) {
    #' select all cal points between highest en lowest linear points
    cal <- cal %>%
      mutate(linearity = ((Area_ / cal_conc) * pnt) / pnt_mean)
    cal_max <- max(filter(cal, linearity <= 1.2 & linearity >= 0.8)[ ,4])
    cal_min <- min(filter(cal, linearity <= 1.2 & linearity >= 0.8)[ ,4])
    cal <- cal %>%
      filter(cal_conc >= cal_min & cal_conc <= cal_max)
    #'fit model on calibration curve
    x <- (j-1)*length(btch)+i
    models[[x]] <- lm(Area_ ~ cal_conc, data = cal)
    coeff <- models[[x]]$coefficients
    models[[x]]$name <- str_c(compounds[j], "_", btch[i]) 
    model_x_stats <- tibble(name = str_c(compounds[j]),
                            batch = str_c(btch[i]),
                            intercept = coeff[1],
                            b = coeff[2],
                            r_sqrd = summary(models[[x]])$r.squared)
    models_stats <- bind_rows(models_stats, model_x_stats)
    curve_max[x] <- max(cal$cal_conc)
    curve_min[x] <- min(cal$cal_conc)
    curvelim <- tibble(compound = compounds[j],
                       batch = btch[[i]],
                       loq = curve_min[[x]],
                       limit = curve_max[[x]])
    limit <- bind_rows(limit, curvelim)
        #'calculate final concentrations
    #' basic formula for dilution during preparation process: (weight sample / (ACN_ml * 1.1)) / 2
    #' an other dilution ratio is included for the dilution batch, this is 1 for all not diluted samples
    temp <- temp_c %>%
      filter(batch == btch[i]) %>%
      mutate(detec_conc = (Area_ - coeff[1]) / coeff[2],
             samp_conc = if_else(samp_type == "std", detec_conc * 8, detec_conc * 2),
             samp_conc = if_else(samp_type == "S" | samp_type == "C", detec_conc / (aimed_w_sample/(ACN_ml*1.1)/2), samp_conc),
             samp_conc = if_else(week %in% dil_wks, samp_conc * dil_ratio, samp_conc),
             LOQ = if_else(detec_conc < curve_min[x], "<LOQ", ">"),
             dilute = if_else(detec_conc > curve_max[x], 1, 0),
             dil_factor = detec_conc / as.numeric(curve_max[x]))
    } else {
      temp <- temp_c %>%
        filter(batch == btch[i]) %>%
        mutate(samp_conc = NaN,
               LOQ = NA,
               dilute = NA,
               dil_factor = NaN, 
               detec_conc = NaN)
    }
    #'calculate Ion Ratio
    temp_ir <- cal %>%
      filter(str_detect(an_name, "Std 0", negate = T))
    m_ir <- mean(temp_ir$ir, na.rm = T)
    min_ir <- m_ir * 0.7
    max_ir <- m_ir * 1.3

    temp <- temp %>%
      mutate(IR_error = if_else(ir > max_ir | ir < min_ir, 1, 0))
        ### calculate recovery ---------------
    #' recovery for calibration curve is compared with the detection conc. For the QC it is with the 
    #' concentration in the sample tubes (not corrected for water). The QC samples are corrected for already
    #' available concentration in the mother material.
    origin_s <- temp %>%
      select(an_name, samp_conc) %>%
      rename(origin = an_name, samp_conc_o = samp_conc) %>%
      mutate(samp_conc_o = if_else(samp_conc_o < 0, 0, samp_conc_o))
    QC_sample <- temp %>%
      filter(an_type == "QC") %>%
      mutate(origin = if_else(str_detect(an_name, "\\+ IS"), "BLANK - IS", str_extract(notes, "\\d\\d.*$")),
             origin = if_else(str_detect(an_name, "IS-std"), "XY-blank-IS", origin)) %>%
      left_join(origin_s, by = "origin") %>%
      select(an_name, origin, samp_conc_o)
   #all QC samples have a real_conc of 50 ng/ml
    recov_temp <- temp %>%
      filter(an_type == "curve" | an_type == "QC" | an_type == "point") %>%
      left_join(QC_sample, by = "an_name") %>%
      distinct(an_code, .keep_all = T) %>%
      mutate(real_conc = if_else(an_type == "curve" | an_type == "point",
                                 as.numeric(str_extract(an_name, "( \\d+\\.\\d*)|( \\d+)")) * 8, 50),
             real_conc = if_else(samp_type == "W", 25, real_conc),
             recov = if_else(an_type == "curve" | an_type == "point", samp_conc / real_conc, 
                             (samp_conc - samp_conc_o) / real_conc),
             orig_eff = samp_conc_o / real_conc,
             recov = if_else(orig_eff > 0.3, NaN, recov)) # when original conc is >30% QC is not valid
    batch_recov[[i]] <- recov_temp %>%
      filter(an_type == "QC") %>%
      mutate(recov = if_else(IR_error == 1, NaN, recov)) %>%
      select(recov, orig_eff) %>%
      rename_with(.cols = everything(), ~ str_c(. , "_", compounds[j]))

  ## VWC adjustment -----------------------------------------
    #' adjust concentration for moisture content in sample. 
    water_s <- temp %>%
      filter(samp_type == "W") %>%
      select(pest_ID, samp_conc, LOQ) %>%
      rename(samp_conc_w = samp_conc,
             w_LOQ = LOQ)
    
    temp <- temp %>%
      left_join(water_s, by = "pest_ID") %>%
      mutate(samp_w = dry_w - empty_w - aimed_w_sample,
             samp_w = if_else(an_name == "20190819_C_L-LC", samp_w - 2.5, samp_w),
             cor_w_part = if_else(w_LOQ == "<LOQ"| is.na(w_LOQ), 0, 
                                  (((vwc_gly * aimed_w_sample)) / 0.998 * samp_conc_w)/ aimed_w_sample),
             samp_conc = signif(if_else(is.na(samp_w) | is.na(cor_w_part), samp_conc, 
                                        ((aimed_w_sample / samp_w) * samp_conc) - cor_w_part),
                                digits = 3))
    
    
    ### Store concentration and quality ---------------------------------
    #' prepare one column with final concentrations, give -999 for all samples that are <LOQ,
    #' and -777 for samples with IR error,
    #' give -888 to samples that need to be diluted    
    qual_list[[i]] <- temp %>%
      filter(an_type == "sample") %>%
      mutate(Q_code = if_else(IR_error == 1, -777, 1),
             Q_code = if_else(LOQ == "<LOQ", -999, Q_code),
             Q_code = if_else(dilute == 1, -888, Q_code),
             dilution = if_else(dilute == 1, detec_conc / curvelim$limit[1], 0),
             final_conc = if_else(IR_error == 1, -777, samp_conc),
             final_conc = signif(if_else(LOQ == "<LOQ", -999, final_conc), digits = 3),
             final_conc = if_else(dilute == 1, -888, final_conc)) %>%
      select(Q_code, detec_conc, samp_conc, dilution, final_conc)
    batch_list[[i]] <- temp %>%
      filter(an_type == "sample") %>%
      mutate(final_conc = if_else(IR_error == 1, -777, samp_conc),
             final_conc = signif(if_else(LOQ == "<LOQ" & !(week %in% dil_wks), -999, final_conc), digits = 3),
             final_conc = if_else(dilute == 1 & !(week %in% dil_wks), -888, final_conc),
             final_conc = if_else(dilute == 1 & week %in% dil_wks, -888, final_conc)) %>%
      select(final_conc) 
    names(batch_list[[i]]) <- str_c("conc_", compounds[j])
    names(qual_list[[i]]) <- c(str_c("qual_", compounds[j]), str_c("detec_", compounds[j]), 
                               str_c("conc_", compounds[j]), str_c("dilution_", compounds[j]),
                               str_c("final_", compounds[j]))
  } # end batch loop
  df_list[[j]] <- bind_rows(batch_list)
  df_recov[[j]] <- bind_rows(batch_recov)
  df_qual[[j]] <- bind_rows(qual_list)
} # end compound loop
pest_result <- bind_cols(results, df_list)
qual_result <- bind_cols(results, df_qual)
pest_recovery <- bind_cols(recov, df_recov)

# Recovery and uncertainty QC -------------
#'make overview of recovery for point data lc_multi and gly ampa
# calculate recovery statistics
recov_w <- pest_recovery %>%
  filter(str_detect(an_name, "-W-")) %>%
  select(all_of(starts_with("recov"))) %>%
  pivot_longer(cols = everything()) %>%
  mutate(recov = mean(value),
         sd = sd(value),
         n_out = sum(value < 0.8 | value > 1.2))
pest_recov_s <- pest_recovery %>%
  filter(!str_detect(an_name, "-W-"),
         !str_detect(an_name, "std_mix|std MIX"),
         !str_detect(an_name, "BLANK")) %>%
  select(all_of(starts_with("recov")))
recov_stats <- tibble(compound = compounds,
                      mean = round(sapply(pest_recov_s, mean, na.rm = T), digits =  1),
                      sd = round(sapply(pest_recov_s, sd, na.rm = T), digits = 4),
                      n = sapply(pest_recov_s, function(x){sum(!is.na(x))}))
Pest_sd <- max(recov_stats$sd) * 100

models_stats <- bind_rows(models_stats, model_gly)
model_mean_stats <- models_stats %>%
  group_by(name) %>%
  mutate(r_sqrd_m = mean(r_sqrd), n = n())

# analysis of different recovery for different QC classes
recov_analysis <- pest_recovery %>%
  mutate(class = if_else(str_detect(an_name, "QC-"), "QC", "BLANK"),
         class = if_else(str_detect(an_name, "-RS-"), "RS", class),
         class = if_else(str_detect(an_name, "std_mix|std MIX"), "STD", class),
         class = if_else(str_detect(an_name, "-W-"), "W", class)) %>%
  group_by(class) %>%
  summarise(across(starts_with("recov"), ~ mean(.x, na.rm = T))) %>%
  pivot_longer(cols = all_of(str_subset(names(pest_recov_s), "recov")), names_to = "compound", 
               values_to = "recovery") %>%
  pivot_wider(names_from = class, values_from = recovery) %>%
  mutate(compound = str_remove(compound, "recov_"))

# Calculate output tables -----------------

#' save results of batch 3 and 4 to calculate dilutions
batch_3_4 <- filter(qual_result, !(week %in% dil_wks)) %>%
  arrange(pest_ID)

### calculate dilutions ------------------------
#' overview dilution of samples
#' samples that have a response higher than the linear part of the calibration curve can not be measured accurately
#' these samples have to be diluted and measured again. Below an overview is made of the amount of samples per compound
#' and of the amount of limit exceedances per sample.
#' Due to an error in the first version of dilution calculations, some wrong samples were diluted. 
#' The corrections are calculated below.

# overview of dilutions per samples
dilute1 <- function(x) (if_else(x == -888, 1, 0))
temp_dil <- batch_3_4 %>%
  mutate(across(matches("qual_"), dilute1)) %>%
  replace(is.na(.), 0) %>%
  select(-(an_name:week), -(all_of(str_subset(names(batch_3_4), "conc_|detec_|dilution_|final_"))))
names(temp_dil) <- compounds
nm_cols <- batch_3_4[1:4]
lc_dilute <- temp_dil %>%
  mutate(sum = rowSums(.)) %>%
  bind_cols(nm_cols) %>%
  select(an_name, sum, everything())
# overview of dilutions per compound
comp_dil <-  temp_dil %>%
  summarise_all(sum) %>%
  pivot_longer(cols = all_of(compounds), names_to = "compound", values_to = "count") %>%
  mutate(compound = str_replace(compound, "conc_", ""))
limit$limit <- as.numeric(limit$limit)
dilute2 <- limit %>%
  group_by(compound) %>%
  mutate(loq = max(loq)) %>%
  pivot_wider(names_from = batch, values_from = limit, names_prefix = "batch_") %>%
  left_join(comp_dil, by = "compound") %>%
  mutate(batch_03 = as.numeric(batch_03),
         batch_04 = as.numeric(batch_04),
         batch_1203 = as.numeric(batch_1203),
         batch_1204 = as.numeric(batch_1204))

#### find wrong dilutions batch 12----------------
#' dilute_samples.csv is the old version, on which the actual preparation of batch12 was based. 
#' See code version before 26-03-2021 to recreate result.
dilute_old <- read_csv("sources/dilute_samples.csv") %>%
  arrange(an_name)
dilute_new <- select(lc_dilute, -week) %>%
  arrange(an_name)
dilute_ns <- select(dilute_new, an_name, sum)
dilute_os <- select(dilute_old, an_name, sum)
names(dilute_os) <- c("an_name", "sum_o")
dilute_compare <- left_join(dilute_ns, dilute_os, by = "an_name")
dilute_error <- filter(dilute_compare, sum != sum_o)

# 0 = not dilute, 1 = must be diluted not done, 2 = diluted not needed, 3 = diluted and done
dilute_diff <- (select(dilute_old, -an_name, -pest_ID, -sum, -samp_type, -'Trinexapac-ethyl', -Tebuconazole) * 2) + select(dilute_new, -an_name, -pest_ID, -sum, -samp_type)
dilutions_error_count <- table(unlist(dilute_diff))
dilute_diff <- bind_cols(arrange(nm_cols, an_name), dilute_diff) %>%
  left_join(dilute_ns, by = "an_name")

### update diluted samples ---------------
#' update the final results with new diluted values, check if results are correct
#' select(-(all_of(str_subset(names(qual_result), "conc_|detec_")))) %>%
result_batches <- c("12", "18")
dil_batches <- vector("list", length = 2)
dil_batches[[1]] <- filter(qual_result, week == "12") %>%
  arrange(pest_ID, samp_type, week)
diluted <- qual_result %>%
  filter(!(week %in% dil_wks)) %>%
  semi_join(dil_batches[[1]], by = "an_name") %>%
  arrange(pest_ID, samp_type, week)
not_diluted <- anti_join(qual_result, dil_batches[[1]], by = "an_name") %>%
  filter(week != "18") %>%
  arrange(pest_ID, samp_type, week)

#'the solution of updating the final results with the diluted values works,
#'but is very slow

samples <- dil_batches[[1]]$an_name
comp_list <- vector("list", length = length(compounds))
samp_list <- vector("list", length = length(samples))
for (i in seq_along(compounds)) {
  dil_comp <- select(diluted, an_name, str_subset(names(dil_batches[[1]]), str_c(compounds[i], "$")))
  batch_comp <- select(dil_batches[[1]], an_name, str_subset(names(dil_batches[[1]]), str_c(compounds[i], "$")))
  for (j in seq_along(samples)) {
    dilute <- filter(dil_comp, an_name == samples[j])
    batch_samp<- batch_comp %>%
      filter(an_name == samples[j])
    if (dilute[1,2] == -888 & !is.na(dilute[1,2])) {
      dilute = batch_samp}
    samp_list[[j]] <- dilute[-1]
  }
  comp_list[[i]] <- bind_rows(samp_list)
}


batch_cols <- dil_batches[[1]][1:4]
diluted_upd <- bind_cols(comp_list) %>%
  mutate(an_name = samples) %>%
  left_join(batch_cols, by = "an_name")

qual_result2 <- bind_rows(not_diluted, diluted_upd) %>%
  arrange(pest_ID, samp_type, week)

### samples still to dilute after batch12 ------------------
#overview dilution per sample
dil_ratio_12 <- lc_data %>%
  select(an_name, week, dil_ratio)
temp_dil <- qual_result2 %>%
  mutate(across(matches("qual_"), dilute1)) %>%
  replace(is.na(.), 0) %>%
  left_join(dil_ratio_12, by = c("an_name", "week"))
lc_dilute_v2 <- temp_dil %>%
  rowwise() %>%
  mutate(sum = sum(c_across(all_of(str_subset(names(temp_dil), "qual_")))))
# overview dilutions per compound
comp_dil <-  temp_dil %>%
  select(all_of(str_subset(names(temp_dil), "qual_"))) %>%
  summarise_all(sum) %>%
  pivot_longer(cols = all_of(str_subset(names(temp_dil), "qual_")), names_to = "compound", values_to = "count") %>%
  mutate(compound = str_replace(compound, "qual_", ""))
limit$limit <- as.numeric(limit$limit)
dilute2_v2 <- limit %>%
  group_by(compound) %>%
  mutate(loq = max(loq)) %>%
  pivot_wider(names_from = batch, values_from = "limit", names_prefix = "batch_") %>%
  left_join(comp_dil, by = "compound") %>%
  mutate(batch_03 = as.numeric(batch_03),
         batch_04 = as.numeric(batch_04),
         batch_1203 = as.numeric(batch_1203),
         batch_1204 = as.numeric(batch_1204))

## LOQ in diluted samples ----------
#'(which is not logical because previously they where too high)
dil_new <- filter(lc_dilute, sum > 0)
b12_l <- semi_join(dil_batches[[1]], dil_new, by = "an_name") %>% arrange(an_name) %>% 
  select(all_of(str_subset(names(.), "qual_")))
dv2_l <- semi_join(dil_new, dil_batches[[1]], by = "an_name") %>% arrange(an_name) %>% 
  select(c(3:33))
loq_check <- b12_l * dv2_l

# check division LOQ, dilute, 0 and 1.
table(unlist(loq_check))
sum(unlist(loq_check) > 0, na.rm = T)

#select samples that are below LOQ and need attention
loq_check <-  loq_check %>%
  rowwise() %>%
  mutate(Min = min(c_across(), na.rm = T)) %>%
  bind_cols(select(arrange(semi_join(dil_new, dil_batches[[1]], by = "an_name"), an_name), an_name)) %>%
  filter(Min == -999)

### update diluted samples batch 18 ---------------
#' update the final results with new diluted values, check if results are correct
#' select(-(all_of(str_subset(names(qual_result), "conc_|detec_")))) %>%
dil_batches[[2]] <- filter(qual_result, week == "18") %>%
  arrange(pest_ID, samp_type, week)
diluted <- qual_result2 %>%
  semi_join(dil_batches[[2]], by = "an_name") %>%
  arrange(pest_ID, samp_type, week)
not_diluted <- anti_join(qual_result2, dil_batches[[2]], by = "an_name") %>%
  arrange(pest_ID, samp_type, week)
batch_cols <- unique(dil_batches[[2]][1:4])
samples <- batch_cols$an_name
comp_list <- vector("list", length = length(compounds))
samp_list <- vector("list", length = length(samples))
for (i in seq_along(compounds)) {
  dil_comp <- select(diluted, an_name, str_subset(names(dil_batches[[2]]), str_c(compounds[i], "$")))
  batch_comp <- select(dil_batches[[2]], an_name, str_subset(names(dil_batches[[2]]), str_c(compounds[i], "$")))
  for (j in seq_along(samples)) {
    dilute <- filter(dil_comp, an_name == samples[j])
    batch_samp<- batch_comp %>%
      filter(an_name == samples[j]) %>%
      filter_at(2, all_vars(. == max(.))) %>%
      filter_at(3, all_vars(. == max(.)))
    if (dilute[1,2] == -888 & !is.na(dilute[1,2])) {
      dilute = batch_samp}
    samp_list[[j]] <- dilute[-1]
  }
  comp_list[[i]] <- bind_rows(samp_list)
}

diluted_upd <- bind_cols(comp_list) %>%
  mutate(an_name = samples) %>%
  left_join(batch_cols, by = "an_name")

qual_result3 <- bind_rows(not_diluted, diluted_upd) %>%
  arrange(pest_ID, samp_type, week)

### samples still to dilute after batch18 ------------------
# no dilutions are needed!
sum(unlist(qual_result3) == -888)

## IR error ------------
sum(unlist(qual_result3) == -777)

#save table to prepare dilutions batch 18
write_csv(lc_dilute_v2, "dilute_samples_batch18.csv")
write_csv(dilute2_v2, "dilute_compounds_batch18.csv")
write_csv(loq_check, "dilute_loq_batch18.csv")

# Final values --------------

#' correct prothioconazole
#merge GLY and multi data and correct for recovery
lc_multi <- qual_result3 %>%
  select(an_name:week, all_of(str_subset(names(.), "final_"))) %>%
  rename_with(~str_replace(., "final_", "conc_")) 

d1 <- as.matrix(lc_multi[5:35])
v1 <- as.matrix(recov_stats[2])

lc_matrix <- as_tibble(sweep(d1, 2, v1, FUN = "/")) %>%
  bind_cols(lc_multi[1:4]) %>%
  mutate(across(starts_with("conc_"), ~ if_else(. < -500, -999, .))) %>%
  rename_with(~str_replace(., "conc_Prothioconazole_desthio-2", "conc_Prothio2")) %>%
  mutate(conc_Prothioconazole = if_else(conc_Prothio2 > 0 & conc_Prothioconazole < 0, 
                                        ((344.3/312.2)* conc_Prothio2), conc_Prothioconazole),
         conc_Prothioconazole = if_else(conc_Prothio2 > 0 & conc_Prothioconazole > 0,
                                        conc_Prothioconazole + ((344.3/312.2)* conc_Prothio2),
                                        conc_Prothioconazole)) %>%
  select(-conc_Prothio2)

lc_all <- left_join(lc_matrix, lc_gly, by = c("pest_ID", "samp_type"))
recov_stats <- bind_rows(recov_stats, recov_gly)

# Save data -------------------
#' dilute_samples_v2.csv = the selection of samples which need dilution
#' dilute_compounds_v2.csv = overview of dilutions per compound
#' dilute_diff_v1.csv = the difference between the wrong dilution calculations and the new one (v2)
#' dilute_samples_batch12.csv = overview of all samples that still need dilution after batch12.
#' dilution_ratio_batch12.csv = overview of amount of dilution needed per sample
#' lc_all_data.csv = values for each sample and compound in ppb
#' lc_all_quality.csv = quality indicator for each sample and compound 
#'           (<LOQ = -999, IR_error = -777, dilute = -888, Good = 1)
write_csv(lc_dilute, "dilute_samples_v2.csv")
write_csv(dilute2, "dilute_compounds_v2.csv")
write_csv(dilute_diff, "dilute_diff_v1.csv")
write_csv(limit, "lc_multi_limits.csv")
write_csv(lc_all, "lc_all_data.csv")
write_csv(recov_stats, "recovery_overview.csv")
