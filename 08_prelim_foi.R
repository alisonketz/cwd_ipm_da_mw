######################################################################################
###
### Preliminary constants for running in the model
###
######################################################################################

### Number of age effects for survival
nT_age_surv <- max(d_surv$right_age_s, na.rm = TRUE) - 1
nT_age_surv_aah_f <- intvl_step_yr_weekly * n_agef + 2
nT_age_short_f <- intvl_step_yr_weekly * (n_agef - 1) + 2
nT_age_surv_aah_m <- intvl_step_yr_weekly * n_agem + 1 
nT_age_short_m <- intvl_step_yr_weekly * (n_agem - 1) + 1
### adding 1 and 2 to account for the fact that deer that live 
### to the 9 or 7 years or older actually live more than
### 52 weeks in a year, due to leap years and timing of weeks beginning into end

nT_age_foi <- max(d_surv$right_age_smonth, na.rm = TRUE) - 1
nT_age_surv_aah_f_foi <- intvl_step_yr_monthly * n_agef + 2
nT_age_short_f_foi <- intvl_step_yr_monthly * (n_agef - 1) + 2
nT_age_surv_aah_m_foi <- intvl_step_yr_monthly * n_agem + 1 
nT_age_short_m_foi <- intvl_step_yr_monthly * (n_agem - 1) + 1


####################################################################################
###
### calculating age_lookup for AAH population model, 
### which is the ageclass given the number of weeks age
###
####################################################################################

# age_lookup_f <- c(rep(1:4, each = intvl_step_yr_weekly),
#                        rep(5, 2 * intvl_step_yr_weekly),
#                        rep(6, 3 * intvl_step_yr_weekly))
# age_lookup_f <- c(age_lookup_f,
#                   rep(7, nT_age_surv - length(age_lookup_f)))

# age_lookup_m <- c(rep(1:4, each = intvl_step_yr_weekly),
#                        rep(5, 2 * intvl_step_yr_weekly),
#                        rep(6, 3 * intvl_step_yr_weekly))
# age_lookup_m <- c(age_lookup_m,
#                   rep(6, nT_age_surv - length(age_lookup_m)))

age_lookup_f <- c(rep(1:4, each = intvl_step_yr_monthly),
                       rep(5, 2 * intvl_step_yr_monthly),
                       rep(6, 3 * intvl_step_yr_monthly))
age_lookup_f <- c(age_lookup_f,
                  rep(7, nT_age_foi - length(intvl_step_yr_monthly)))

age_lookup_m <- c(rep(1:4, each = intvl_step_yr_monthly),
                       rep(5, 2 * intvl_step_yr_monthly),
                       rep(6, 3 * intvl_step_yr_monthly))
age_lookup_m <- c(age_lookup_m,
                  rep(6, nT_age_foi - length(age_lookup_m)))


######################################################
### age to date conversion for FOI age/period effects
######################################################
death_end <- "2022-05-14"
first_birth <- "1992-05-15"
cwd_df$birthweek <- floor(interval(first_birth,
                            cwd_df$birth_date) / weeks(1)) + 1
cwd_df$birthmonth <- floor(interval(first_birth,
                            cwd_df$birth_date) / months(1)) + 1
cwd_df$weekkill <- floor(interval(study_origin,
                            cwd_df$kill_date) / weeks(1))
cwd_df$yearkill <- cwd_df$kill_year - year(study_origin) + 1

# interval("1994-05-15","1995-01-01") %/% weeks(1)
######################################################################

end_dates <- paste(2002:2022,"05","14",sep="-")
interval_cuts <- floor(interval(study_origin,end_dates)/weeks(1))
(intvl_step_yearly <- c(diff(interval_cuts),52))

period_lookup_foi <- rep(1, ceiling(interval("1992-05-15","2002-05-15") / weeks(1)))
for(t in 2:n_year){
  period_lookup_foi <- c(period_lookup_foi,
                         rep(t,intvl_step_yearly[t]))
}
#making sure the period effect lookup vector is sufficiently long
period_lookup_foi <- c(period_lookup_foi,rep(n_year,max(cwd_df$birthweek+cwd_df$ageweeks) - length(period_lookup_foi)))
(n_period_lookup <- length(period_lookup_foi))

nT_period_prestudy <- floor(interval("1992-05-15","2001-05-15") / weeks(1))
period_lookup_foi_hunt <- period_lookup_foi
period_lookup_foi <- period_lookup_foi[-(1:nT_period_prestudy)]
period_lookup_foi <- period_lookup_foi[1:nT_period_overall]
# max(d_surv$left_period_e+(d_surv$right_age_s - d_surv$left_age_e),na.rm=TRUE)
nT_period_overall_hunt <- length(period_lookup_foi_hunt)


# num_foi_cal <- table(period_lookup_foi_study)

#############################################################################################
### Creating adjacency matrix and hyper parameter
### values for the dcar_normal implementation
### for period effects for FOI 
#############################################################################################

#create num vector
num_period <- c(1, rep(2, n_year - 2), 1)

#create adjacency vector along both years
temp <- as.matrix(bandSparse(n = n_year, k = c(1), symmetric = T))
temp2 <- matrix(0, n_year, n_year)
for (i in 1:nrow(temp2)) {
  temp2[i, ] <- temp[i, ] * c(1:n_year)
}
adj_period <- t(temp2)[which(t(temp2) != 0)]
n_adj_period <- length(adj_period)
weights_period <- rep(1, length(adj_period))

#########################################################
###
### Age and period indexing for FOI for collared deer
###
#########################################################

# age_week_indx <- c(rep(1,intvl_step_yr_weekly),#fawns
#                    rep(2,intvl_step_yr_weekly),#1
#                    rep(3,intvl_step_yr_weekly),#2
#                    rep(4,intvl_step_yr_weekly),#3
#                    rep(5, intvl_step_yr_weekly * 2),#4-5
#                    rep(6, intvl_step_yr_weekly * 3),#6-8
#                    rep(7,nT_age_surv - 
#                          length(c(rep(1, intvl_step_yr_weekly),#fawns
#                                   rep(2, intvl_step_yr_weekly),#1
#                                   rep(3, intvl_step_yr_weekly),#2
#                                   rep(4, intvl_step_yr_weekly),#3
#                                   rep(5, intvl_step_yr_weekly * 2),#4-5
#                                   rep(6, intvl_step_yr_weekly * 3)))))#6_

# period_week_indx <- c(rep(1,51),#2017
#                    rep(2,52),#2018
#                    rep(3,52),#2019
#                    rep(4,52),#2020
#                    rep(5,52),#2021
#                    rep(6,nT_period_collar - length(c(rep(1,51),#2017
#                                               rep(2,52),#2018
#                                               rep(3,52),#2019
#                                               rep(4,52),#2020
#                                               rep(5,52))))#2022
#                    )

# period_week_indx_col <- period_week_indx + n_year_precollar - 1
