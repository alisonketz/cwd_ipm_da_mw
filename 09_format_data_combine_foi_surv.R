
################################################################
###
### separating survival data into likelihood specific classes
###
################################################################

# all deer
low_sus <- d_surv$lowtag

# uninfected at recapture, antemortem negative
# low_recap_neg 

# infected at recapture,
# low_recap_pos

# infected at capture, antemortem
low_icap <- d_surv$lowtag[d_surv$cwd_cap == 1]

# test + infected at mortality but not at capture or recapture
# remove infected at mortality from infected antemort
# remove deer that died and weren't tested at mortality
low_idead <- d_surv$lowtag[d_surv$cwd_mort == 1 & !is.na(d_surv$cwd_mort) & !is.na(d_surv$right_period_s)]
low_idead <- low_idead[!(low_idead %in% c(low_recap, low_icap))]

### all individuals that are alive at the end of the study that weren't recaptured
# low_endlive <- d_surv$lowtag[!(d_surv$lowtag %in% 
#          as.integer(unique(c(d_cens$lowtag, d_mort$lowtag))))]
# low_endlive 

# never tested positive, and were tested at
low_sus <- low_sus[!(low_sus %in% unique(c(low_icap,
                                           low_recap_pos,
                                           low_recap_neg,
                                           low_idead,
                                           low_endlive)))]

d_fit_sus <- d_surv[d_surv$lowtag %in% low_sus, ]
d_fit_icap <- d_surv[d_surv$lowtag %in% low_icap, ]
d_fit_idead <- d_surv[d_surv$lowtag %in% low_idead, ]
d_fit_rec_pos <- d_surv[d_surv$lowtag %in% low_recap_pos, ]
d_fit_rec_neg <- d_surv[d_surv$lowtag %in% low_recap_neg, ]
d_fit_endlive <- d_surv[d_surv$lowtag %in% low_endlive, ]

# there is more than one record for deer that die during the study
# so the number of individuals is less than the number of records for each
# class of animals

n_fit_sus <- nrow(d_fit_sus)
n_fit_icap <- nrow(d_fit_icap)
n_fit_idead <- nrow(d_fit_idead)
n_fit_rec_neg <- nrow(d_fit_rec_neg)
n_fit_rec_pos <- nrow(d_fit_rec_pos)
n_fit_endlive <- nrow(d_fit_endlive)

######################
### Recaptures
######################

# separating the censor part from the mortality part
# so to properly structure the data for the likelihood
# keeping sus and inf parts separated

# there is 1 recaptured deer that was right censored and not a mortality and did not live to end of study
# due to a dead battery
# this deer was test negative at capture, test positive at recapture, then interval censored

d_fit_rec_pos_cens <- d_fit_rec_pos[d_fit_rec_pos$lowtag == 6065, ]
d_fit_rec_pos <- d_fit_rec_pos[d_fit_rec_pos$lowtag != 6065, ]
n_fit_rec_pos <- nrow(d_fit_rec_pos)
n_fit_rec_pos_cens <- nrow(d_fit_rec_pos_cens)

# for the cwd test (-) deer at capture
# these deer were interval censored and tested neg at capture and at censoring (or not tested at censoring)
d_fit_rec_neg_cens <- d_fit_rec_neg[d_fit_rec_neg$censor == 1, ]
n_fit_rec_neg_cens <- nrow(d_fit_rec_neg_cens)

# these deer were test negative at capture, and recapture, and test negative at mort,
# all of them had post mortem tests
d_fit_rec_neg_mort <- d_fit_rec_neg[d_fit_rec_neg$censor == 0 & d_fit_rec_neg$cwd_mort == 0, ]
n_fit_rec_neg_mort <- nrow(d_fit_rec_neg_mort)

# these deer were test negative at capture, and recapture, and test positive at censor
d_fit_rec_pos_mort <- d_fit_rec_neg[d_fit_rec_neg$censor == 0, ]
n_fit_rec_pos_mort <- nrow(d_fit_rec_pos_mort)


############################################
### separating susceptible deer that were 
### test negative at capture into whether 
### tested after mortality or censoring
############################################

d_fit_sus_cens <- d_fit_sus[d_fit_sus$censored == 1, ]
d_fit_sus_mort <- d_fit_sus[d_fit_sus$censored == 0, ]

d_fit_sus_cens_posttest <- d_fit_sus_cens[!is.na(d_fit_sus_cens$cwd_mort), ]
d_fit_sus_cens_postno <- d_fit_sus_cens[is.na(d_fit_sus_cens$cwd_mort), ]

d_fit_sus_mort_posttest <- d_fit_sus_mort[!is.na(d_fit_sus_mort$cwd_mort), ]
d_fit_sus_mort_postno <- d_fit_sus_mort[is.na(d_fit_sus_mort$cwd_mort), ]

d_fit_rec_neg_cens_posttest <- d_fit_rec_neg_cens[!is.na(d_fit_rec_neg_cens$cwd_mort), ]
d_fit_rec_neg_cens_postno <- d_fit_rec_neg_cens[is.na(d_fit_rec_neg_cens$cwd_mort), ]


############################################
### adjusting data frames for the 
### deer that enter the study as fawns and 
### die within the first month are assumed to 
### be test negative at mortality
############################################

adj_indx <- which(d_fit_sus_mort_postno$right_age_s - 1 < 5)

d_fit_sus_mort_posttest <- rbind(d_fit_sus_mort_posttest,
                                 d_fit_sus_mort_postno[adj_indx, ])
d_fit_sus_mort_postno <- d_fit_sus_mort_postno[-adj_indx, ]

####################################################
### separating infected-at-capture 
### censored deervs mortalities
#####################################################

d_fit_icap_cens <- d_fit_icap[d_fit_icap$censored == 1, ]
d_fit_icap_mort <- d_fit_icap[d_fit_icap$censored == 0, ]

####################################################
### sample sizes each likelihood type
#####################################################

n_fit_sus_cens_posttest <- nrow(d_fit_sus_cens_posttest)
n_fit_sus_cens_postno <- nrow(d_fit_sus_cens_postno)
n_fit_sus_mort_posttest <- nrow(d_fit_sus_mort_posttest)
n_fit_sus_mort_postno <- nrow(d_fit_sus_mort_postno)
n_fit_icap_cens <- nrow(d_fit_icap_cens)
n_fit_icap_mort <- nrow(d_fit_icap_mort)
n_fit_rec_neg_cens_posttest <- nrow(d_fit_rec_neg_cens_posttest)
n_fit_rec_neg_cens_postno <- nrow(d_fit_rec_neg_cens_postno)

#############################################################################################
# separating Hunter Harvest Data for likelihoods into positive or negative cases
#############################################################################################

###############################################################
###
### Setting up aggregation by study
### Separating hunter harvested
### CWD test positive deer from CWD test negative
###
### d_fit_hunt_pos
### d_fit_hunt_neg
###
##############################################################

# cwd_df_agg <- cwd_df %>% 
#               group_by(teststatus, ageweeks, birthweek, sex, ew) %>%
#               summarise(n_cases = n(), .groups = 'drop')
# d_fit_hunt_neg <- cwd_df_agg[cwd_df_agg$teststatus == 0, ]
# d_fit_hunt_pos <- cwd_df_agg[cwd_df_agg$teststatus == 1, ]

# n_fit_hunt_pos <- nrow(d_fit_hunt_pos)
# n_fit_hunt_neg <- nrow(d_fit_hunt_neg)

# d_fit_hunt <- rbind(d_fit_hunt_neg,d_fit_hunt_pos)
# n_fit_hunt <- nrow(d_fit_hunt)

cwd_df_agg <- cwd_df %>% 
              group_by(teststatus, agemonths, birthmonth, sex, ew) %>%
              summarise(n_cases = n(), .groups = 'drop')
d_fit_hunt_neg <- cwd_df_agg[cwd_df_agg$teststatus == 0, ]
d_fit_hunt_pos <- cwd_df_agg[cwd_df_agg$teststatus == 1, ]

n_fit_hunt_pos <- nrow(d_fit_hunt_pos)
n_fit_hunt_neg <- nrow(d_fit_hunt_neg)

d_fit_hunt <- rbind(d_fit_hunt_neg,d_fit_hunt_pos)
n_fit_hunt <- nrow(d_fit_hunt)



####################################################
### For the deer that are fast mortalities, 
### setting up an indicator to change 
### which likelihood calculation to prevent bad indexing
### sum(d_fit_icap_mort$left_age_e == d_fit_icap_mort$right_age_r)
#####################################################

# d_fit_sus_mort_posttest$fast <- ifelse(d_fit_sus_mort_posttest$left_age_e ==
#                                        d_fit_sus_mort_posttest$right_age_r,1,0)


# d_fit_icap_mort$fast <- ifelse(d_fit_icap_mort$left_age_e ==
#                                        d_fit_icap_mort$right_age_r,1,0)


###############################################################
###
### all data cases
###
###############################################################

# num_observations <- rbind(n_fit_hunt_neg,
#         n_fit_hunt_pos,
#         n_fit_sus_cens_posttest,
#         n_fit_sus_cens_postno,
#         n_fit_sus_mort_posttest,
#         n_fit_sus_mort_postno,
#         n_fit_icap_cens,
#         n_fit_icap_mort,
#         n_fit_rec_neg_cens_posttest,
#         n_fit_rec_neg_cens_postno,
#         n_fit_rec_neg_mort,
#         n_fit_rec_pos_cens,
#         n_fit_rec_pos_mort,
#         n_fit_idead,
#         n_fit_endlive)


# data_cases <- rownames(rbind(n_fit_hunt_neg,
#         n_fit_hunt_pos,
#         n_fit_sus_cens_posttest,
#         n_fit_sus_cens_postno,
#         n_fit_sus_mort_posttest,
#         n_fit_sus_mort_postno,
#         n_fit_icap_cens,
#         n_fit_icap_mort,
#         n_fit_rec_neg_cens_posttest,
#         n_fit_rec_neg_cens_postno,
#         n_fit_rec_neg_mort,
#         n_fit_rec_pos_cens,
#         n_fit_rec_pos_mort,
#         n_fit_idead,
#         n_fit_endlive)
#         )

# obs_sample_sizes <- data.frame(data_cases, num_observations, row.names = NULL)

# write.csv(obs_sample_sizes, file = "results/obs_sample_sizes.csv", row.names=FALSE)

# d_fit_icap_cens$left_period_e - d_fit_icap_cens$left_age_e

# png("figures/icap_birth_relative_study.png")
# hist(d_fit_icap_mort$left_period_e - d_fit_icap_mort$left_age_e,breaks=100)
# dev.off()


# obs_sample_sizes_desc <- read.csv("../obs_sample_sizes_description.csv")
# print(xtable(obs_sample_sizes_desc),include.rownames=FALSE)


d_fit_rec_pos <- rbind(d_fit_rec_pos_mort,d_fit_rec_pos_cens)
n_fit_rec_pos <- nrow(d_fit_rec_pos)



################################################################
###
### reaggregating survival data into likelihood specific classes
###
################################################################

### never get sick, no draw for survival
d_sus <- rbind(d_fit_sus_cens_posttest,
                   d_fit_sus_mort_posttest,
                   d_fit_rec_neg_mort,
                   d_fit_rec_neg_cens_posttest)

### infected at capture, no draw
d_icap <- rbind(d_fit_icap_mort, d_fit_icap_cens)


### infected at recapture, draw from entry to recapture
d_recap  <- rbind(d_fit_rec_pos_cens, d_fit_rec_pos_mort)

### infected at mortality, draw from entry to mort
d_idead  <- d_fit_idead

### no test at mortality or censoring, so do 2 part approach
### first, draw whether they get infected
### second, draw the interval during which that ocurrs
### these are all right censored so only 1 line
d_fit_sus_draw <- rbind(
    d_fit_sus_cens_postno,
    d_fit_endlive
)

#these are mortalities, so if draw, then we need 2 rows to correct 
d_sus_mort_postno <- d_fit_sus_mort_postno

### no test at mortality or censoring, so do 2 part approach
### first, draw whether they get infected
### second, draw the interval during which that ocurrs
### this type the interval that it could occur is from
### recapture until right censoring 
### (rather than entry to right censoring)

d_fit_rec_neg_cens_postno

##########################################################################
###
### Formatting these six data types to fit into models with the 
### 2 row formatting, i.e. a row for censoring and a row for mortality
###
##########################################################################

d_fit_sus <- d_sus
d_fit_sus$left <- d_sus$left_age_e
d_fit_sus$right <- d_sus$right_age_r
d_fit_sus$surv_censor <- 1
for(i in 1:nrow(d_sus)) {
    if(!is.na(d_sus$right_period_s[i])){
        d_fit_sus <- rbind(d_fit_sus,data.frame(cbind(d_sus[i,],left = NA,right=NA,surv_censor = 0)))
        d_fit_sus$left[nrow(d_fit_sus)] <- d_sus$right_age_r[i]
        d_fit_sus$right[nrow(d_fit_sus)] <- d_sus$right_age_s[i]
    }
}
d_fit_sus <- d_fit_sus[order(d_fit_sus$lowtag),]
d_fit_sus$teststatus <- 0
d_fit_sus$age2date  <-  d_fit_sus$left_period_e - d_fit_sus$left_age_e
n_fit_sus <- nrow(d_fit_sus)

d_fit_icap <- d_icap
d_fit_icap$left <- d_icap$left_age_e
d_fit_icap$right <- d_icap$right_age_r
d_fit_icap$surv_censor <- 1
for(i in 1:nrow(d_icap)) {
    if(!is.na(d_icap$right_period_s[i])){
        d_fit_icap <- rbind(d_fit_icap,data.frame(cbind(d_icap[i,],left = NA,right=NA,surv_censor = 0)))
        d_fit_icap$left[nrow(d_fit_icap)] <- d_icap$right_age_r[i]
        d_fit_icap$right[nrow(d_fit_icap)] <- d_icap$right_age_s[i]
    }
}
d_fit_icap <- d_fit_icap[order(d_fit_icap$lowtag),]
d_fit_icap$teststatus <- 1
d_fit_icap$age2date <- d_fit_icap$left_period_e - d_fit_icap$left_age_e

d_fit_recap <- d_recap
d_fit_recap$left <- d_recap$left_age_e
d_fit_recap$right <- d_recap$right_age_r
d_fit_recap$surv_censor <- 1
for(i in 1:nrow(d_recap)) {
    if(!is.na(d_recap$right_period_s[i])){
        d_fit_recap <- rbind(d_fit_recap,data.frame(cbind(d_recap[i,],left = NA,right=NA,surv_censor = 0)))
        d_fit_recap$left[nrow(d_fit_recap)] <- d_recap$right_age_r[i]
        d_fit_recap$right[nrow(d_fit_recap)] <- d_recap$right_age_s[i]
    }
}
d_fit_recap <- d_fit_recap[order(d_fit_recap$lowtag),]
d_fit_recap$teststatus <- 1

d_fit_idead <- d_idead
d_fit_idead$left <- d_idead$left_age_e
d_fit_idead$right <- d_idead$right_age_r
d_fit_idead$surv_censor <- 1
for(i in 1:nrow(d_idead)) {
    if(!is.na(d_idead$right_period_s[i])){
        d_fit_idead <- rbind(d_fit_idead,data.frame(cbind(d_idead[i,],left = NA,right=NA,surv_censor = 0)))
        d_fit_idead$left[nrow(d_fit_idead)] <- d_idead$right_age_r[i]
        d_fit_idead$right[nrow(d_fit_idead)] <- d_idead$right_age_s[i]
    }
}
d_fit_idead <- d_fit_idead[order(d_fit_idead$lowtag),]
d_fit_idead$teststatus <- 1

d_fit_sus_draw$left <- d_fit_sus_draw$left_age_e
d_fit_sus_draw$right <- d_fit_sus_draw$right_age_r
d_fit_sus_draw$surv_censor <- 1

d_fit_sus_mort_postno <- d_sus_mort_postno
d_fit_sus_mort_postno$left <- d_sus_mort_postno$left_age_e
d_fit_sus_mort_postno$right <- d_sus_mort_postno$right_age_r
d_fit_sus_mort_postno$surv_censor <- 1
for(i in 1:nrow(d_sus_mort_postno)) {
        d_fit_sus_mort_postno <- rbind(d_fit_sus_mort_postno,data.frame(cbind(d_sus_mort_postno[i,],left = NA,right=NA,surv_censor = 0)))
        d_fit_sus_mort_postno$left[nrow(d_fit_sus_mort_postno)] <- d_sus_mort_postno$right_age_r[i]
        d_fit_sus_mort_postno$right[nrow(d_fit_sus_mort_postno)] <- d_sus_mort_postno$right_age_s[i]
}
d_fit_sus_mort_postno <- d_fit_sus_mort_postno[order(d_fit_sus_mort_postno$lowtag),]

d_fit_rec_neg_cens_postno$left <- d_fit_rec_neg_cens_postno$left_age_e
d_fit_rec_neg_cens_postno$right <- d_fit_rec_neg_cens_postno$right_age_r
d_fit_rec_neg_cens_postno$surv_censor <- 1

################################################################
### Arranging data for FOI
################################################################

d_fit_hunt <- rbind(d_fit_hunt_neg,
                    d_fit_hunt_pos)
n_fit_hunt <- nrow(d_fit_hunt)

d_fit_sus_foi <- d_fit_sus[d_fit_sus$surv_censor == 1,]
for(i in 1:nrow(d_fit_sus_foi)) {
    if(!is.na(d_fit_sus_foi$right_age_smonth[i])){
        d_fit_sus_foi$right[i] <- d_fit_sus_foi$right_age_smonth[i]
    }
}

# d_fit_sus_foi$left[d_fit_sus_foi$lowtag %in% 
#                     c(d_fit_rec_neg_mort$lowtag,
#                     d_fit_rec_neg_cens_posttest$lowtag)]  <- 
#         d_fit_sus_foi$ageweek_recap[d_fit_sus_foi$lowtag %in% 
#                     c(d_fit_rec_neg_mort$lowtag,
#                     d_fit_rec_neg_cens_posttest$lowtag)]

d_fit_sus_foi$left <- d_fit_sus_foi$left_age_month
d_fit_sus_foi$left[d_fit_sus_foi$lowtag %in% 
                    c(d_fit_rec_neg_mort$lowtag,
                    d_fit_rec_neg_cens_posttest$lowtag)]  <- 
        d_fit_sus_foi$agemonth_recap[d_fit_sus_foi$lowtag %in% 
                    c(d_fit_rec_neg_mort$lowtag,
                    d_fit_rec_neg_cens_posttest$lowtag)]

# d_fit_sus_foi$age2date <- d_fit_sus_foi$left_period_e - d_fit_sus_foi$left_age_e
# n_fit_sus_foi <- nrow(d_fit_sus_foi)

d_fit_sus_foi$age2date <- d_fit_sus_foi$emonth - d_fit_sus_foi$left_age_month
n_fit_sus_foi <- nrow(d_fit_sus_foi)

d_fit_icap_foi <- d_fit_icap[d_fit_icap$surv_censor == 1,]
d_fit_icap_foi$left <- d_fit_icap_foi$left_age_month
d_fit_icap_foi$right <- d_fit_icap_foi$right_age_rmonth
# d_fit_icap_foi$age2date <- d_fit_icap_foi$left_period_e - d_fit_icap_foi$left_age_e
d_fit_icap_foi$age2date <- d_fit_icap_foi$emonth - d_fit_icap_foi$left_age_month
n_fit_icap_foi <- nrow(d_fit_icap_foi)
d_fit_icap_foi$teststatus

d_fit_rec_neg_cens_postno_foi <- d_fit_rec_neg_cens_postno
d_fit_rec_neg_cens_postno_foi$teststatus <- 0
# d_fit_rec_neg_cens_postno_foi$right <- d_fit_rec_neg_cens_postno_foi$ageweek_recap
d_fit_rec_neg_cens_postno_foi$right <- d_fit_rec_neg_cens_postno_foi$agemonth_recap

# d_fit_rec_neg_cens_postno_foi$age2date <- d_fit_rec_neg_cens_postno_foi$left_period_e -
#                                           d_fit_rec_neg_cens_postno_foi$left_age_e
d_fit_rec_neg_cens_postno_foi$age2date <- d_fit_rec_neg_cens_postno_foi$emonth -
                                          d_fit_rec_neg_cens_postno_foi$left_age_month
n_fit_rec_neg_cens_postno_foi <- nrow(d_fit_rec_neg_cens_postno_foi)

d_fit_recap_foi <- d_fit_recap[d_fit_recap$surv_censor == 1,]
d_fit_recap_foi$teststatus <- 1
d_fit_recap_foi$left <- d_fit_recap_foi$left_age_month
d_fit_recap_foi$right <- d_fit_recap_foi$agemonth_recap
# d_fit_recap_foi$right <- d_fit_recap_foi$ageweek_recap

d_fit_recap_foi$age2date <- d_fit_recap_foi$emonth - d_fit_recap_foi$left_age_month
# d_fit_recap_foi$age2date <- d_fit_recap_foi$left_period_e - d_fit_recap_foi$left_age_e

n_fit_recap_foi <- nrow(d_fit_recap_foi)

d_fit_idead_foi <- d_fit_idead[d_fit_idead$surv_censor == 1,]
d_fit_idead_foi$teststatus <- 1
# d_fit_idead_foi$right <- d_fit_idead_foi$right_age_s
d_fit_idead_foi$left <- d_fit_idead_foi$left_age_month
d_fit_idead_foi$right <- d_fit_idead_foi$right_age_smonth
# d_fit_idead_foi$age2date <- d_fit_idead_foi$left_period_e - d_fit_idead_foi$left_age_e
d_fit_idead_foi$age2date <- d_fit_idead_foi$emonth - d_fit_idead_foi$left_age_month

n_fit_idead_foi <- nrow(d_fit_idead_foi)

#############################################
#############################################
#############################################
###
###  Parameters for fitting survival models
###
#############################################
#############################################
#############################################

d_fit_idead <- d_fit_idead[d_fit_idead$surv_censor == 0,]
d_fit_idead$age2date <- d_fit_idead$left_period_e - d_fit_idead$left_age_e
n_fit_idead <- nrow(d_fit_idead)

d_fit_rec_pos_cens$age2date <- d_fit_rec_pos_cens$left_period_e - d_fit_rec_pos_cens$left_age_e
n_fit_rec_pos_cens <- nrow(d_fit_rec_pos_cens)

d_fit_rec_pos_mort$age2date <- d_fit_rec_pos_mort$left_period_e - d_fit_rec_pos_mort$left_age_e
n_fit_rec_pos_mort <- nrow(d_fit_rec_pos_mort)

d_fit_sus_draw$age2date <- d_fit_sus_draw$left_period_e - d_fit_sus_draw$left_age_e
n_fit_sus_draw <- nrow(d_fit_sus_draw)

d_fit_sus_mort_postno$age2date <- d_fit_sus_mort_postno$left_period_e - d_fit_sus_mort_postno$left_age_e
n_fit_sus_mort_postno <- nrow(d_fit_sus_mort_postno)


d_fit_rec_neg_cens_postno$age2date <- d_fit_rec_neg_cens_postno$left_period_e - d_fit_rec_neg_cens_postno$left_age_e
n_fit_rec_neg_cens_postno <- nrow(d_fit_rec_neg_cens_postno)



# #fudging fix periodweek recap is 1 week off for this individual
# d_fit_rec_pos_cens$periodweek_recap <- d_fit_rec_pos_cens$periodweek_recap -1

### the distance between period week and recap period week is different than the distance (in intervals) of ageweek recap to age at entry
d_fit_rec_neg_cens_postno$periodweek_recap[which(d_fit_rec_neg_cens_postno$ageweek_recap - 
        d_fit_rec_neg_cens_postno$left_age_e !=
      d_fit_rec_neg_cens_postno$periodweek_recap - 
      d_fit_rec_neg_cens_postno$left_period_e)] <- d_fit_rec_neg_cens_postno$periodweek_recap[which(d_fit_rec_neg_cens_postno$ageweek_recap - 
        d_fit_rec_neg_cens_postno$left_age_e !=
      d_fit_rec_neg_cens_postno$periodweek_recap - 
      d_fit_rec_neg_cens_postno$left_period_e)] + 1

# d_fit_rec_neg_cens_postno$ageweek_recap - d_fit_rec_neg_cens_postno$left_age_e 
# d_fit_rec_neg_cens_postno$periodweek_recap - d_fit_rec_neg_cens_postno$left_period_e


### the distance between period week and recap period week is different than the distance (in intervals) of ageweek recap to age at entry
d_fit_rec_pos_mort$periodweek_recap[which(d_fit_rec_pos_mort$ageweek_recap - 
        d_fit_rec_pos_mort$left_age_e !=
      d_fit_rec_pos_mort$periodweek_recap - 
      d_fit_rec_pos_mort$left_period_e)] <- d_fit_rec_pos_mort$periodweek_recap[which(d_fit_rec_pos_mort$ageweek_recap - 
        d_fit_rec_pos_mort$left_age_e !=
      d_fit_rec_pos_mort$periodweek_recap - 
      d_fit_rec_pos_mort$left_period_e)] + 1

d_fit_rec_pos_mort$ageweek_recap - d_fit_rec_pos_mort$left_age_e 
d_fit_rec_pos_mort$periodweek_recap - d_fit_rec_pos_mort$left_period_e



