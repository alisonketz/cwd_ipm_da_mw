############################################################################################
############################################################################################
############################################################################################
###
### Model Statement
###
############################################################################################
############################################################################################
############################################################################################

modelcode <- nimbleCode({

  ##############################
  ### Priors
  ##############################
  
  beta_male ~ dnorm(0, .1)

  ##############################
  ### Force of infection model
  ##############################
  tau_age_foi_male  ~ dgamma(1, 1)
  tau_age_foi_female  ~ dgamma(1, 1)
  tau1_age_foi_male <- .0000001 * tau_age_foi_male
  tau1_age_foi_female <- .0000001 * tau_age_foi_female
  m_age_foi[1] ~ dnorm(0, tau1_age_foi_male)
  f_age_foi[1] ~ dnorm(0, tau1_age_foi_female)
  m_age_foi[2] ~ dnorm(0, tau1_age_foi_male)
  f_age_foi[2] ~ dnorm(0, tau1_age_foi_female)
  for (i in 3:n_ageclassm) {
    m_age_foi[i]~dnorm(2 * m_age_foi[i-1] - m_age_foi[i-2], tau_age_foi_male)
  }
  for (i in 3:n_ageclassf) {
    f_age_foi[i]~dnorm(2 * f_age_foi[i-1] - f_age_foi[i-2], tau_age_foi_female)
  }
  m_age_foi_mu <- mean(m_age_foi[1:n_ageclassm])
  f_age_foi_mu <- mean(f_age_foi[1:n_ageclassf])

  # Period effects 
  # tau_period_foi_male  ~ dgamma(1, 1)
  # tau_period_foi_female  ~ dgamma(1, 1)

  ###ICAR specification
  # f_period_foi[1:n_year] ~ dcar_normal(adj = adj_period[1:n_adj_period],
  #                                 weights = weights_period[1:n_adj_period],
  #                                 num = num_period[1:n_year],
  #                                 tau = tau_period_foi_female,
  #                                 zero_mean = 1)
  
  # m_period_foi[1:n_year] ~ dcar_normal(adj = adj_period[1:n_adj_period],
  #                                 weights = weights_period[1:n_adj_period],
  #                                 num = num_period[1:n_year],
  #                                 tau = tau_period_foi_male,
  #                                 zero_mean = 1)
  ### RW1 Specification
  # tau1_period_foi_f <- .0000001 * tau_period_foi_female
  # tau1_period_foi_m <- .0000001 * tau_period_foi_male
  # f_period_foi_temp[1] ~ dnorm(0, tau1_period_foi_f)
  # m_period_foi_temp[1] ~ dnorm(0, tau1_period_foi_m)
  # for (t in 2:n_year) {
  #   f_period_foi_temp[t] ~ dnorm(f_period_foi_temp[t - 1], tau_period_foi_female)
  #   m_period_foi_temp[t] ~ dnorm(m_period_foi_temp[t - 1], tau_period_foi_male)
  # }
  # f_period_foi_mu <- mean(f_period_foi_temp[1:n_year])
  # m_period_foi_mu <- mean(m_period_foi_temp[1:n_year])
  # for (t in 1:n_year) {
  #   f_period_foi[t] <- f_period_foi_temp[t] - f_period_foi_mu
  #   m_period_foi[t] <- m_period_foi_temp[t] - m_period_foi_mu
  # }
  
  ### removing first years prior to surveillance data
  # tau1_period_foi_f <- .0000001 * tau_period_foi_female
  # tau1_period_foi_m <- .0000001 * tau_period_foi_male
  # f_period_foi_temp[9] ~ dnorm(0, tau1_period_foi_f)
  # m_period_foi_temp[9] ~ dnorm(0, tau1_period_foi_m)
  # for (t in 10:n_year) {
  #   f_period_foi_temp[t] ~ dnorm(f_period_foi_temp[t - 1], tau_period_foi_female)
  #   m_period_foi_temp[t] ~ dnorm(m_period_foi_temp[t - 1], tau_period_foi_male)
  # }
  # f_period_foi_mu <- mean(f_period_foi_temp[9:n_year])
  # m_period_foi_mu <- mean(m_period_foi_temp[9:n_year])
  # for (t in 9:n_year) {
  #   f_period_foi[t] <- f_period_foi_temp[t] - f_period_foi_mu
  #   m_period_foi[t] <- m_period_foi_temp[t] - m_period_foi_mu
  # }


  ### random effect for East/West spatial model
  space[1] <- 0
  space[2]  <- space_temp * space_mix

  space_temp ~ dnorm(0, 1)
  space_mix ~ dunif(-1, 1)

  ############################################################
  ############################################################
  ### Age/Period Survival Model
  ############################################################
  ############################################################

  ####################################
  ### Susceptibles survival intercept
  ####################################
  
  # beta0_survival_sus ~ dnorm(0, .1)
  # beta0_survival_sus ~ T(dnorm(-6, .1),,0)
  beta0_sus_temp ~ dnorm(0, .1)
  sus_mix ~ dunif(-1, 1)
  beta0_survival_sus <- beta0_sus_temp * sus_mix

  ##################################
  ### Infected survival intercept
  ##################################

  # beta0_survival_inf ~ dnorm(0, .1)
  # beta0_survival_inf ~ T(dnorm(-6, .1),,0)
  beta0_inf_temp ~ dnorm(0, .1)
  inf_mix ~ dunif(-1, 1)
  beta0_survival_inf <- beta0_inf_temp * inf_mix

  ##################################
  ### Priors for Age and Period effects
  ##################################
  #Age effects
  for (k in 1:nknots_age) {
    # b_age_survival[k] ~ dnorm(0, tau_age_survival)
    ln_b_age_survival[k] ~ dnorm(0, tau_age_survival)
    b_age_survival[k] <- exp(ln_b_age_survival[k])
  }
  tau_age_survival ~ dgamma(1, 1)

  for (t in 1:nT_age_surv) {
    age_effect_survival_temp[t] <- inprod(b_age_survival[1:nknots_age],
                                     Z_age[t, 1:nknots_age])
    age_effect_survival[t] <-  age_effect_survival_temp[t] -
                               mu_age_effect_survival_temp
  }
  mu_age_effect_survival_temp <- mean(age_effect_survival_temp[1:nT_age_surv])

  #Period effects from collar data
  for (k in 1:nknots_period) {
    b_period_survival[k] ~ dnorm(0, tau_period_survival)
  }
  tau_period_survival ~ dgamma(1, 1)
  for (t in 1:nT_period_collar) {
    period_effect_surv[t] <- inprod(b_period_survival[1:nknots_period],
                                    Z_period[t, 1:nknots_period])
  }

  #Period effects from aah data
  tau_period_precollar ~ dgamma(1,1)
  for (k in 1:n_year_precollar) {
    period_annual_survival[k] ~ dnorm(0, tau_period_precollar)
  }

  period_effect_survival[1:nT_period_overall] <- set_period_effects_constant(
        n_year_precollar = n_year_precollar,
        nT_period_precollar = nT_period_precollar,
        nT_period_collar = nT_period_collar,
        nT_period_overall = nT_period_overall,
        yr_start = yr_start[1:n_year],
        yr_end = yr_end[1:n_year],
        period_effect_surv = period_effect_surv[1:nT_period_collar],
        period_annual_survival = period_annual_survival[1:n_year_precollar]
  )

  #######################################################################
  #######################################################################
  ## Likelihoods of FOI
  #######################################################################
  #######################################################################

  #######################################################################
  ###
  ###   User defined distribution for likelihood for
  ###   all harvested deer
  ###
  ###   d_fit_hunt
  ###
  #######################################################################

  y_hunt_foi ~ dFOIhunt(
          y = y_hunt_foi_teststatus[1:n_fit_hunt],
          n_cases = n_cases_foi_hunt[1:n_fit_hunt],
          n_samples = n_fit_hunt,
          a = a_fit_hunt[1:n_fit_hunt],
          sex = sex_fit_hunt[1:n_fit_hunt],
          age2date = age2date_fit_hunt[1:n_fit_hunt],
          f_age = f_age_foi[1:n_ageclassf],
          m_age = m_age_foi[1:n_ageclassm],
          age_lookup_f = age_lookup_f[1:n_age_lookup_f],
          age_lookup_m = age_lookup_m[1:n_age_lookup_m],
          period_lookup_foi = period_lookup_foi_hunt[1:nT_period_overall_hunt],
          f_period = f_period_foi[1:n_year],
          m_period = m_period_foi[1:n_year],
          space = space[1:n_study_area],
          sect = sect_hunt[1:n_fit_hunt]
  )

  ######################################################################
  ##
  ##   User defined distribution for likelihood for FOI of
  ##   collared deer that stay susceptible
  ##
  ##   d_fit_sus_foi
  ##
  ######################################################################

  for (i in 1:n_fit_sus_foi){
    y_sus_foi_teststatus[i] ~ dFOIcollar(left = left_sus_foi[i],
            right = right_sus_foi[i],
            sex = sex_sus_foi[i],
            age2date = age2date_sus_foi[i],
            f_age = f_age_foi[1:n_ageclassf],
            m_age = m_age_foi[1:n_ageclassm],
            age_lookup_f = age_lookup_f[1:n_age_lookup_f],
            age_lookup_m = age_lookup_m[1:n_age_lookup_m],
            period_lookup = period_lookup_foi[1:nT_period_overall],
            f_period = f_period_foi[1:n_year],
            m_period = m_period_foi[1:n_year],
            space = space[sect_sus_foi[i]])
  }

  #######################################################################
  ###
  ###   User defined distribution for likelihood for FOI of
  ###   collared deer that are infected at capture
  ###
  ###   d_fit_icap
  ###
  #######################################################################

  for (i in 1:n_fit_icap_foi){
    y_icap_foi_teststatus[i] ~ dFOIicap(left = left_icap_foi[i],
            sex = sex_icap_foi[i],
            age2date = age2date_icap_foi[i],
            f_age = f_age_foi[1:n_ageclassf],
            m_age = m_age_foi[1:n_ageclassm],
            age_lookup_f = age_lookup_f[1:n_age_lookup_f],
            age_lookup_m = age_lookup_m[1:n_age_lookup_m],
            period_lookup = period_lookup_foi[1:nT_period_overall],
            f_period = f_period_foi[1:n_year],
            m_period = m_period_foi[1:n_year],
            space = space[sect_icap_foi[i]])
  }



  #######################################################################
  ###
  ###   User defined distribution for likelihood for FOI of
  ###   collared deer that test neg at entry, then recap and test neg,
  ###    then no post mortem test (exclude that part)
  ###
  ###   d_fit_rec_neg_cens_postno_foi
  ###
  #######################################################################

  for (i in 1:n_fit_rec_neg_cens_postno_foi){
    y_rec_neg_cens_postno_foi_teststatus[i] ~ dFOIcollar(
      left = left_rec_neg_cens_postno_foi[i],
      right = right_rec_neg_cens_postno_foi[i],
      sex = sex_rec_neg_cens_postno_foi[i],
      age2date = age2date_rec_neg_cens_postno_foi[i],
      f_age = f_age_foi[1:n_ageclassf],
      m_age = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup_foi[1:nT_period_overall],
      f_period = f_period_foi[1:n_year],
      m_period = m_period_foi[1:n_year],
      space = space[sect_rec_neg_cens_postno_foi[i]])
  }


  #######################################################################
  ###
  ###   User defined distribution for likelihood for FOI of
  ###   collared deer that test neg at entry, then recap and test neg,
  ###    then no post mortem test (exclude that part)
  ###
  ###   d_fit_recap_foi
  ###
  #######################################################################

  for (i in 1:n_fit_recap_foi){
    y_recap_foi_teststatus[i] ~ dFOIcollar(
      left = left_recap_foi[i],
      right = right_recap_foi[i],
      sex = sex_recap_foi[i],
      age2date = age2date_recap_foi[i],
      f_age = f_age_foi[1:n_ageclassf],
      m_age = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup_foi[1:nT_period_overall],
      f_period = f_period_foi[1:n_year],
      m_period = m_period_foi[1:n_year],
      space = space[sect_recap_foi[i]])
  }
  #######################################################################
  ###
  ###   User defined distribution for likelihood for FOI of
  ###   collared deer that test neg at entry, then recap and test neg,
  ###    then no post mortem test (exclude that part)
  ###
  ###   d_fit_idead_foi
  ###
  #######################################################################

  for (i in 1:n_fit_idead_foi){
    y_idead_foi_teststatus[i] ~ dFOIcollar(
      left = left_idead_foi[i],
      right = right_idead_foi[i],
      sex = sex_idead_foi[i],
      age2date = age2date_idead_foi[i],
      f_age = f_age_foi[1:n_ageclassf],
      m_age = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup_foi[1:nT_period_overall],
      f_period = f_period_foi[1:n_year],
      m_period = m_period_foi[1:n_year],
      space = space[sect_idead_foi[i]])
  }

  #######################################################################
  #######################################################################
  ## Likelihoods of Survival
  #######################################################################
  #######################################################################


  #######################################################################
  ###
  ###   User defined distribution for likelihood for Survival of
  ###   collared deer that test neg at entry and postmortem
  ###
  ###   d_fit_sus
  ###
  #######################################################################

  for (i in 1:n_fit_sus){
    y_sus_surv[i] ~ dSurvival(
      left = left_sus[i],
      right = right_sus[i],
      sex = sex_sus[i],
      age2date = age2date_sus[i],
      age_effect = age_effect_survival[1:nT_age_surv],
      period_effect = period_effect_survival[1:nT_period_overall],
      nT_age = nT_age_surv,
      beta0 = beta0_survival_sus,
      beta_male = beta_male)
  }

  #######################################################################
  ###
  ###   User defined distribution for likelihood for Survival of
  ###   collared deer that test pos at entry
  ###
  ###   d_fit_icap
  ###
  #######################################################################

  for (i in 1:n_fit_icap){
    y_icap_surv[i] ~ dSurvival(
      left = left_icap[i],
      right = right_icap[i],
      sex = sex_icap[i],
      age2date = age2date_icap[i],
      age_effect = age_effect_survival[1:nT_age_surv],
      period_effect = period_effect_survival[1:nT_period_overall],
      nT_age = nT_age_surv,
      beta0 = beta0_survival_inf,
      beta_male = beta_male)
  }

##################################################################
###
###  User Defined Distribution Age-Period Survival Model
###  for infected at mortality, must draw when infected
###
###  d_fit_idead
###
##################################################################

#   for (i in 1:n_fit_idead) {
#     y_idead[i] ~ dSurvival_idead(
#         e_age = idead_left_age_e[i],
#         r_age = idead_right_age_r[i],
#         s_age = idead_right_age_s[i],
#         e_period = idead_left_period_e[i],
#         s_period = idead_right_period_s[i],
#         sex = idead_sex[i],
#         age2date = idead_age2date[i],
#         age_effect_survival = age_effect_survival[1:nT_age_surv],
#         period_effect_survival = period_effect_survival[1:nT_period_overall],
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi[1:n_ageclassf],
#         m_age_foi = m_age_foi[1:n_ageclassm],
#         f_period_foi = f_period_foi[1:n_year],
#         m_period_foi = m_period_foi[1:n_year],
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi[1:nT_period_overall],
#         space = space[sect_idead[i]],
#         age_lookup_f = age_lookup_f[1:nT_age_surv],
#         age_lookup_m = age_lookup_m[1:nT_age_surv]
#         )
#   }


# #   ##################################################################
# #   ###
# #   ###  User Defined Distribution Age-Period Survival Model
# #   ###  for infected at recapture and was then right censored
# #   ###  must draw when infected
# #   ###
# #   ###  d_fit_rec_pos_cens
# #   ###  (there's only 1)
# #   ###
# #   ##################################################################

#     y_rec_pos_cens ~ dSurvival_rec_pos_cens(
#         e_age = rec_pos_cens_left_age_e,
#         r_age = rec_pos_cens_right_age_r,
#         recap_age = rec_pos_cens_ageweek_recap,
#         e_period = rec_pos_cens_left_period_e,
#         recap_period = rec_pos_cens_periodweek_recap,
#         sex = rec_pos_cens_sex,
#         age2date = rec_pos_cens_age2date,
#         age_effect_survival = age_effect_survival[1:nT_age_surv],
#         period_effect_survival = period_effect_survival[1:nT_period_overall],
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi[1:n_ageclassf],
#         m_age_foi = m_age_foi[1:n_ageclassm],
#         f_period_foi = f_period_foi[1:n_year],
#         m_period_foi = m_period_foi[1:n_year],
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi[1:nT_period_overall],
#         space = space[sect_rec_pos_cens],
#         age_lookup_f = age_lookup_f[1:nT_age_surv],
#         age_lookup_m = age_lookup_m[1:nT_age_surv])


#   ##################################################################
#   ###
#   ###  User Defined Distribution Age-Period Survival Model
#   ###  for infected at recapture and was then right censored
#   ###  must draw when infected
#   ###
#   ###  d_fit_rec_pos_mort
#   ###
#   ##################################################################

#   for(i in 1:n_fit_rec_pos_mort){
#     y_rec_pos_mort[i] ~ dSurvival_rec_pos_mort(
#         e_age = rec_pos_mort_left_age_e[i],
#         r_age = rec_pos_mort_right_age_r[i],
#         s_age = rec_pos_mort_right_age_s[i],
#         recap_age = rec_pos_mort_ageweek_recap[i],
#         e_period = rec_pos_mort_left_period_e[i],
#         recap_period = rec_pos_mort_periodweek_recap[i],
#         sex = rec_pos_mort_sex[i],
#         age2date = rec_pos_mort_age2date[i],
#         age_effect_survival = age_effect_survival[1:nT_age_surv],
#         period_effect_survival = period_effect_survival[1:nT_period_overall],
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi[1:n_ageclassf],
#         m_age_foi = m_age_foi[1:n_ageclassm],
#         f_period_foi = f_period_foi[1:n_year],
#         m_period_foi = m_period_foi[1:n_year],
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi[1:nT_period_overall],
#         space = space[sect_rec_pos_mort[i]],
#         age_lookup_f = age_lookup_f[1:nT_age_surv],
#         age_lookup_m = age_lookup_m[1:nT_age_surv])
#   }

#   ###########################################################################
#   ###
#   ###  User Defined Distribution Age-Period Survival Model
#   ###  no test after initial capture, first drawing if gets infected
#   ###  then if it gets infected, accounting for the suscepting and infected
#   ###  portions in the joint likelihood
#   ###
#   ###  d_fit_sus_draw
#   ###
#   ############################################################################

#   for(i in 1:n_fit_sus_draw) {
#       y_sus_draw[i] ~ dSurvival_sus_draw(
#               e_age = sus_draw_left_age_e[i],
#               r_age = sus_draw_right_age_r[i],
#               e_period = sus_draw_left_period_e[i],
#               r_period = sus_draw_right_period_r[i],
#               sex = sus_draw_sex[i],
#               age2date = sus_draw_age2date[i],
#               age_effect_survival = age_effect_survival[1:nT_age_surv],
#               period_effect_survival = period_effect_survival[1:nT_period_overall],
#               nT_age_surv = nT_age_surv,
#               beta0_survival_inf = beta0_survival_inf,
#               beta0_survival_sus = beta0_survival_sus,
#               beta_male = beta_male,
#               f_age_foi = f_age_foi[1:n_ageclassf],
#               m_age_foi = m_age_foi[1:n_ageclassm],
#               f_period_foi = f_period_foi[1:n_year],
#               m_period_foi = m_period_foi[1:n_year],
#               nT_period_overall = nT_period_overall,
#               period_lookup_foi = period_lookup_foi[1:nT_period_overall],
#               space = space[sect_sus_draw[i]],
#               age_lookup_f = age_lookup_f[1:nT_age_surv],
#               age_lookup_m = age_lookup_m[1:nT_age_surv])
#   }

#   # ##################################################################
#   # ###
#   # ###  User Defined Distribution Age-Period Survival Model
#   # ###  no test after initial capture, first drawing if gets infected
#   # ###  then if it gets infected, accounting for the susceptible and infected
#   # ###  portions in the joint likelihood, these individuals died
#   # ###
#   # ###  d_fit_sus_mort_postno
#   # ###
#   # ##################################################################

#   for(i in 1:n_fit_sus_mort_postno) {
#       y_sus_mort_postno[i] ~ dSurvival_sus_mort_postno(
#         e_age = sus_mort_postno_left_age_e[i],
#         r_age = sus_mort_postno_right_age_r[i],
#         s_age = sus_mort_postno_right_age_s[i],
#         e_period = sus_mort_postno_left_period_e[i],
#         s_period = sus_mort_postno_right_period_s[i],
#         sex = sus_mort_postno_sex[i],
#         age2date = sus_mort_postno_age2date[i],
#         age_effect_survival = age_effect_survival[1:nT_age_surv],
#         period_effect_survival = period_effect_survival[1:nT_period_overall],
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi[1:n_ageclassf],
#         m_age_foi = m_age_foi[1:n_ageclassm],
#         f_period_foi = f_period_foi[1:n_year],
#         m_period_foi = m_period_foi[1:n_year],
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi[1:nT_period_overall],
#         space = space[sect_sus_mort_postno[i]],
#         age_lookup_f = age_lookup_f[1:nT_age_surv],
#         age_lookup_m = age_lookup_m[1:nT_age_surv])
#   }

#   #######################################################################
#   ###
#   ###   User defined distribution for likelihood for Survival of
#   ###   collared deer that test pos at entry
#   ###
#   ###  # must DRAW when infection occurs
#   ###   d_fit_rec_neg_cens_postno
#   ###
#   #######################################################################
  
#   for(i in 1:n_fit_rec_neg_cens_postno) {
#     y_rec_neg_cens_postno[i] ~ dSurvival_rec_neg_cens_postno(
#         e_age = rec_neg_cens_postno_left_age_e[i],
#         r_age = rec_neg_cens_postno_right_age_r[i],
#         recap_age = rec_neg_cens_postno_ageweek_recap[i],
#         e_period = rec_neg_cens_postno_left_period_e[i],
#         r_period = rec_neg_cens_postno_right_period_r[i],
#         recap_period = rec_neg_cens_postno_periodweek_recap[i],
#         sex = rec_neg_cens_postno_sex[i],
#         age2date = rec_neg_cens_postno_age2date[i],
#         age_effect_survival = age_effect_survival[1:nT_age_surv],
#         period_effect_survival = period_effect_survival[1:nT_period_overall],
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi[1:n_ageclassf],
#         m_age_foi = m_age_foi[1:n_ageclassm],
#         f_period_foi = f_period_foi[1:n_year],
#         m_period_foi = m_period_foi[1:n_year],
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi[1:nT_period_overall],
#         space = space[sect_rec_neg_cens_postno[i]],
#         age_lookup_f = age_lookup_f[1:nT_age_surv],
#         age_lookup_m = age_lookup_m[1:nT_age_surv])
#   }



  #######################################################
  #######################################################
  #######################################################
  ### Cause-specific mortality model
  #######################################################
  #######################################################
  #######################################################

  ### priors
  ### sex-specific hunt probability given mortalities
  beta0_cause ~ dnorm(0, .1)
  beta_cause_gun ~ dnorm(0, .1)
  beta_cause_ng ~ dnorm(0, .1)
  beta_cause_male ~ dnorm(0, .1)

  #######################################################
  ###
  ### Likelihood for cause of death
  ###
  #######################################################

  for (i in 1:records_cause) {
      p_cause[i]  <- ilogit(beta0_cause +
                            Z_cause_ng[interval_cause[i]] * beta_cause_ng +
                            Z_cause_gun[interval_cause[i]] * beta_cause_gun +
                            sex_cause[i] * beta_cause_male)
      mort_hh[i] ~ dbin(size = 1, prob = p_cause[i])
  }

  #######################################################
  ###
  ### Derived parameter for cause of death by hunter harvest
  ### given 
  ###
  #######################################################

  p_nogun_f <- ilogit(beta0_cause +
                      beta_cause_ng)
  p_gun_f <- ilogit(beta0_cause +
                    beta_cause_ng +
                    beta_cause_gun)
  p_nogun_m <- ilogit(beta0_cause +
                      beta_cause_ng +
                      beta_cause_male)
  p_gun_m <- ilogit(beta0_cause +
                    beta_cause_ng +
                    beta_cause_gun +
                    beta_cause_male)

  # #############################################################################
  # #############################################################################
  # #############################################################################
  # ## Age-at-harvest population model
  # #############################################################################
  # #############################################################################
  # #############################################################################

  # #######################################
  # ### Initial population
  # ########################################
  # # tau_pop ~ dgamma(1, 1)

  # #should this be different for pos/neg or m/f for study area?
  # for(i in 1:2) {
  #   tau_pop[i] ~ dgamma(10, 1)
  # }
  # ####################################
  # ### just setting the initial pop size
  # ####################################
  # # for(k in 1:n_study_area) {
  # #     for (a in 1:n_agef) {

  # #     #Initial population structure pop[sex,age,year] for susceptible deer
  # #     llpop_sus[k, 1, a, 1]  <- f_logpop_sus[k, a]
  # #     pop_sus[k, 1, a, 1] <- exp(llpop_sus[k, 1, a, 1])

  # #     #Initial population structure pop[study_area=k,sex=i,year=t,age=a]
  # #     llpop_inf[k, 1, a, 1] <- f_logpop_inf[k, a]
  # #     pop_inf[k, 1, a, 1] <- exp(llpop_inf[k, 1, a, 1])
  # #     }
  # #     for (a in 1:n_agem) {
  # #         ### Initial population structure pop
  # #         ### [study_area,sex,age,period(year)] for susceptible deer
  # #         llpop_sus[k, 2, a, 1]  <- m_logpop_sus[k, a]
  # #         pop_sus[k, 2, a, 1] <- exp(llpop_sus[k, 2, a, 1])

  # #         #Initial population structure pop for infected deer
  # #         llpop_inf[k, 2, a, 1] <- m_logpop_inf[k, a]
  # #         pop_inf[k, 2, a, 1] <- exp(llpop_inf[k, 2, a, 1])
  # #     }
  # # }

  # ##this one uses the hyperprior structure for initial pop size
  # for(k in 1:n_study_area) {
  #   for (a in 1:n_agef) {

  #     #Initial population structure pop[sex,age,year] for susceptible deer
  #     llpop_sus[k, 1, a, 1] ~ dnorm(f_logpop_sus[k, a], tau_pop[1])#tau_pop[1]
  #     pop_sus[k, 1, a, 1] <- exp(llpop_sus[k, 1, a, 1])

  #     #Initial population structure pop[study_area=k,sex=i,year=t,age=a]
  #     llpop_inf[k, 1, a, 1] ~ dnorm(f_logpop_inf[k, a], tau_pop[2])#tau_pop[1]
  #     pop_inf[k, 1, a, 1] <- exp(llpop_inf[k, 1, a, 1])
  #   }
  #   for (a in 1:n_agem) {
  #       ### Initial population structure pop
  #       ### [study_area,sex,age,period(year)] for susceptible deer
  #       llpop_sus[k, 2, a, 1] ~ dnorm(m_logpop_sus[k, a], tau_pop[1])#tau_pop[2]
  #       pop_sus[k, 2, a, 1] <- exp(llpop_sus[k, 2, a, 1])

  #       #Initial population structure pop for infected deer
  #       llpop_inf[k, 2, a, 1] ~ dnorm(m_logpop_inf[k, a], tau_pop[2])#tau_pop[2]
  #       pop_inf[k, 2, a, 1] <- exp(llpop_inf[k, 2, a, 1])
  #   }
  # }

  # ###########################
  # ###Reporting Rates
  # ###########################

  # report_overall ~ dbeta(report_hyp_all[1], report_hyp_all[2])
  # for(t in 1:12){#2002-2014
  #   report[t]  <- report_overall
  # }
  # for(t in 13:17){ #2015-2020
  #   report[t] ~ dbeta(report_hyp_y[t - 12, 1], report_hyp_y[t - 12, 2])
  # }
  # report[18:20]  <- report_overall #2021

  # ############################
  # #### Fecundity
  # ############################

  # mu_fec ~ dnorm(0, 1)
  # fec_prec_eps ~ dgamma(1, 1)

  # #Observations of fawns & does overall from the 3 counties
  # for(t in 1:n_year_fec_early) {
  #   fec_epsilon[t] ~ dnorm(0, fec_prec_eps)
  #   fec[t] <- exp(mu_fec + fec_epsilon[t])
  #   Nfawn[t] ~ dpois(fec[t] * Ndoe[t])
  # }

  # #for 2017:2021
  # for(t in (n_year_fec_early + 1):n_year){
  #   fec[t] ~ dgamma(obs_ct_fd_alpha[t], obs_ct_fd_beta[t])
  # }

  ##################################################
  ### Overall Survival Susceptibles
  ##################################################

  # sn_sus[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_aah(
  #     nT_age = nT_age_surv,
  #     nT_period_overall = nT_period_overall,
  #     nT_age_short_f = nT_age_short_f,
  #     nT_age_short_m = nT_age_short_m,
  #     nT_age_surv_aah_f = nT_age_surv_aah_f,
  #     nT_age_surv_aah_m = nT_age_surv_aah_m,
  #     beta0 = beta0_survival_sus,
  #     beta_male = beta_male,
  #     age_effect = age_effect_survival[1:nT_age_surv],
  #     period_effect = period_effect_survival[1:nT_period_overall],
  #     yr_start = yr_start[1:n_year],
  #     yr_end = yr_end[1:n_year],
  #     n_year = n_year,
  #     n_agef = n_agef,
  #     n_agem = n_agem)

  # ###################################################
  # #### Overall Survival CWD Infected
  # ###################################################

  # sn_inf[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_aah(
  #     nT_age = nT_age_surv,
  #     nT_period_overall = nT_period_overall,
  #     nT_age_short_f = nT_age_short_f,
  #     nT_age_short_m = nT_age_short_m,
  #     nT_age_surv_aah_f = nT_age_surv_aah_f,
  #     nT_age_surv_aah_m = nT_age_surv_aah_m,
  #     beta0 = beta0_survival_inf,
  #     beta_male = beta_male,
  #     age_effect = age_effect_survival[1:nT_age_surv],
  #     period_effect = period_effect_survival[1:nT_period_overall],
  #     yr_start = yr_start[1:n_year],
  #     yr_end = yr_end[1:n_year],
  #     n_year = n_year,
  #     n_agef = n_agef,
  #     n_agem = n_agem)

  # ###################################################
  # #### Hunting Survival Susceptibles
  # ###################################################

  # sh_sus[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_harvest(
  #     nT_age = nT_age_surv,
  #     nT_period_overall = nT_period_overall,
  #     nT_age_short_f = nT_age_short_f,
  #     nT_age_short_m = nT_age_short_m,
  #     nT_age_surv_aah_f = nT_age_surv_aah_f,
  #     nT_age_surv_aah_m = nT_age_surv_aah_m,
  #     beta0 = beta0_survival_sus,
  #     beta_male = beta_male,
  #     age_effect = age_effect_survival[1:nT_age_surv],
  #     period_effect = period_effect_survival[1:nT_period_overall],
  #     n_sex = n_sex,
  #     n_year = n_year,
  #     n_agef = n_agef,
  #     n_agem = n_agem,
  #     ng_start = ng_start[1:n_year],
  #     gun_start = gun_start[1:n_year],
  #     gun_end = gun_end[1:n_year],
  #     ng_end = ng_end[1:n_year],
  #     yr_start = yr_start[1:n_year],
  #     yr_end = yr_end[1:n_year],
  #     p_nogun_f = p_nogun_f,
  #     p_nogun_m = p_nogun_m,
  #     p_gun_f = p_gun_f,
  #     p_gun_m = p_gun_m
  #     )

  # ###################################################
  # #### Hunting Survival Infected
  # ###################################################

  # sh_inf[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_harvest(
  #     nT_age = nT_age_surv,
  #     nT_period_overall = nT_period_overall,
  #     nT_age_short_f = nT_age_short_f,
  #     nT_age_short_m = nT_age_short_m,
  #     nT_age_surv_aah_f = nT_age_surv_aah_f,
  #     nT_age_surv_aah_m = nT_age_surv_aah_m,
  #     beta0 = beta0_survival_sus,
  #     beta_male = beta_male,
  #     age_effect = age_effect_survival[1:nT_age_surv],
  #     period_effect = period_effect_survival[1:nT_period_overall],
  #     n_sex = n_sex,
  #     n_year = n_year,
  #     n_agef = n_agef,
  #     n_agem = n_agem,
  #     ng_start = ng_start[1:n_year],
  #     gun_start = gun_start[1:n_year],
  #     gun_end = gun_end[1:n_year],
  #     ng_end = ng_end[1:n_year],
  #     yr_start = yr_start[1:n_year],
  #     yr_end = yr_end[1:n_year],
  #     p_nogun_f = p_nogun_f,
  #     p_nogun_m = p_nogun_m,
  #     p_gun_f = p_gun_f,
  #     p_gun_m = p_gun_m
  #     )


  # ###########################################################
  # #### Annual Probability of Infection based on FOI hazards
  # ###########################################################

  # psi[1:n_study_area, 1:n_sex, 1:n_agef, 1:n_year] <-
  #     calc_infect_prob(age_lookup_f = age_lookup_f[1:nT_age_surv],
  #           age_lookup_m = age_lookup_m[1:nT_age_surv],
  #           n_agef = n_agef,
  #           n_agem = n_agem,
  #           yr_start = yr_start[1:n_year],
  #           yr_end = yr_end[1:n_year],
  #           f_age = f_age_foi[1:n_ageclassf],
  #           m_age = m_age_foi[1:n_ageclassm],
  #           f_period = f_period_foi[1:n_year],
  #           m_period = m_period_foi[1:n_year],
  #           nT_period_overall = nT_period_overall,
  #           period_lookup_foi= period_lookup_foi[1:nT_period_overall],
  #           n_year = n_year,
  #           n_sex = n_sex,
  #           n_study_area = n_study_area,
  #           space = space[n_study_area],
  #           nT_age_surv_aah_f = nT_age_surv_aah_f,
  #           nT_age_surv_aah_m = nT_age_surv_aah_m
  #           )

  # ##################################################################
  # ### Probability of Infection from birth pulse to end of harvest
  # ### based on FOI hazards
  # ##################################################################

  # psi_hat[1:n_study_area, 1:n_sex, 1:n_agef, 1:n_year] <-
  #     calc_infect_prob_hunt(age_lookup_f = age_lookup_f[1:nT_age_surv],
  #         age_lookup_m = age_lookup_m[1:nT_age_surv],
  #         n_agef = n_agef,
  #         n_agem = n_agem,
  #         yr_start = yr_start[1:n_year],
  #         yr_end = yr_end[1:n_year],
  #         ng_end = ng_end[1:n_year],
  #         f_age = f_age_foi[1:n_ageclassf],
  #         m_age = m_age_foi[1:n_ageclassm],
  #         f_period = f_period_foi[1:n_year],
  #         m_period = m_period_foi[1:n_year],
  #         nT_period_overall = nT_period_overall,
  #         period_lookup_foi = period_lookup_foi[1:nT_period_overall],
  #         n_year = n_year,
  #         n_sex = n_sex,
  #         n_study_area = n_study_area,
  #         space = space[n_study_area],
  #         nT_age_surv_aah_f = nT_age_surv_aah_f,
  #         nT_age_surv_aah_m = nT_age_surv_aah_m,
  #         fudge_factor = .5
  #         )

  # #################################################
  # ## Earn-a-buck correction factor
  # ## based on Van Deelen et al (2010)
  # #################################################

  # eab_antlerless_temp ~ dgamma(eab_antlerless_alpha,eab_antlerless_beta)
  # eab_antlered_temp ~ dgamma(eab_antlered_alpha,eab_antlered_beta)
  # for(t in 1:n_year) {
  #   eab_antlerless[t] <- eab_antlerless_temp^x_eab[t]
  #   eab_antlered[t]  <- eab_antlered_temp^x_eab[t]
  # }

  #####################################################################
  ##
  ## Population Process Model
  ## Population Projection
  ## pop_proj temporarily holds the projected age class
  ##
  #####################################################################
    #####################################
    ### version with projection matrix
    #####################################

    # for (k in 1:n_study_area) {

    #   ##################
    #   ### Susceptible
    #   ##################

    #   for (t in 2:n_year) {
        

    #     ##############################################################
    #     ##############################################################
    #     ### AAH population model version with a projection matrix
    #     ##############################################################
    #     ##############################################################

    #     #Female: project forward annually
    #     for (a in 1:(n_agef - 2)) {
    #       pop_sus_proj[k, 1, a, t] <- pop_sus[k, 1, a, t - 1] *
    #                                   sn_sus[1, a, t - 1] *
    #                                   (1 - psi[k, 1, a, t - 1])
    #     }

    #     #female accumulating age class
    #     pop_sus_proj[k, 1, n_agef - 1, t] <- pop_sus[k, 1, n_agef - 1, t - 1] *
    #                                         sn_sus[1, n_agef - 1, t - 1] *
    #                                         (1 - psi[k, 1, n_agef - 1, t - 1]) +
    #                                         pop_sus[k, 1,  n_agef, t - 1] *
    #                                         sn_sus[1, n_agef, t - 1] *
    #                                         (1 - psi[k, 1, n_agef, t - 1])

    #     #Female: set projection into population model matrix
    #     for (a in 2:n_agef) {
    #       pop_sus[k, 1, a, t] <- pop_sus_proj[k, 1, (a - 1), t]
    #     }

    #     #Female: fawn class = total #females * unisex fawns per female/2
    #     pop_sus[k, 1, 1, t] <- (sum(pop_sus_proj[k, 1, 1:(n_agef - 1), t]) +
    #                             sum(pop_inf_proj[k, 1, 1:(n_agef - 1), t])) *
    #                             fec[t] * .5
    #     ##########
    #     ### Males
    #     ##########

    #     #Male: project forward anually
    #     for (a in 1:(n_agem - 2)) {
    #       pop_sus_proj[k, 2, a, t] <- pop_sus[k, 2, a, t - 1] *
    #                                   sn_sus[2, a, t - 1] *
    #                                   (1 - psi[k, 2, a, t - 1])
    #     }

    #     #Male: accumulating age class
    #     pop_sus_proj[k, 2, n_agem - 1, t] <- pop_sus[k, 2, n_agem - 1, t - 1] *
    #                                         sn_sus[2, n_agem - 1, t - 1] *
    #                                         (1 - psi[k, 2, n_agem - 1, t - 1]) +
    #                                         pop_sus[k, 2,  n_agem, t - 1] *
    #                                         sn_sus[2, n_agem, t - 1] *
    #                                         (1 - psi[k, 2, n_agem, t - 1])

    #     #Male: set projection into population model matrix
    #     for (a in 2:n_agem) {
    #       pop_sus[k, 2, a, t] <- pop_sus_proj[k, 2, (a - 1), t]
    #     }

    #     # Male: fawn class = total #females * unisex fawns per female/2
    #     pop_sus[k, 2, 1, t] <- pop_sus[k, 1, 1, t]

    #     ###################################################
    #     ### Infected/Infectious
    #     ###################################################

    #     ###########
    #     ### Females
    #     ###########

    #     #Female: project forward anually
    #     for (a in 1:(n_agef - 2)) {
    #       pop_inf_proj[k, 1, a, t] <- pop_inf[k, 1, a, t - 1] *
    #                                   sn_inf[1, a, t - 1] +
    #                                   pop_sus[k, 1, a, t - 1] *
    #                                   sn_sus[1, a, t - 1] *
    #                                   psi[k, 1, a, t - 1]
    #     }
    #     ##Female: accumulating age = 9.5+ years
    #     pop_inf_proj[k, 1, n_agef - 1, t] <- pop_inf[k, 1, n_agef - 1, t - 1] *
    #                                          sn_inf[1, n_agef - 1, t - 1] +
    #                                          pop_inf[k, 1, n_agef, t - 1] *
    #                                          sn_inf[1, n_agef, t - 1] +
    #                                          pop_sus[k, 1, n_agef - 1, t - 1] *
    #                                          sn_sus[1, n_agef - 1, t - 1] *
    #                                          psi[k, 1, n_agef - 1, t - 1] +
    #                                          pop_sus[k, 1,  n_agef, t - 1] *
    #                                          sn_sus[1, n_agef, t - 1] *
    #                                          psi[k, 1, n_agef, t - 1]


    #     #Female: set projection into population model matrix
    #     for (a in 2:n_agef) {
    #       pop_inf[k, 1, a, t] <- pop_inf_proj[k, 1, (a - 1), t]
    #     }
    #     ##Female: fawn class
    #     ##there are no infected fawns at birth
    #     pop_inf[k, 1, 1, t] <- 0

    #     ##########
    #     ### Males
    #     ##########

    #     ### Male: project forward anually
    #     for (a in 1:(n_agem - 2)) {
    #         pop_inf_proj[k, 2, a, t] <- pop_inf[k, 2, a, t - 1] *
    #                                       sn_inf[2, a, t - 1] +
    #                                       pop_sus[k, 2, a, t - 1] *
    #                                       sn_sus[2, a, t - 1] *
    #                                       psi[k, 2, a, t - 1]
    #     }

    #     ### Male: accumulating age class
    #     pop_inf_proj[k, 2, n_agem - 1, t] <- pop_inf[k, 2, n_agem - 1, t - 1] *
    #                                          sn_inf[2, n_agem - 1, t - 1] +
    #                                          pop_inf[k, 2, n_agem, t - 1] *
    #                                          sn_inf[2, n_agem, t - 1] +
    #                                          pop_sus[k, 2, n_agem - 1, t - 1] *
    #                                          sn_sus[2, n_agem - 1, t - 1] *
    #                                          psi[k, 2, n_agem - 1, t - 1] +
    #                                          pop_sus[k, 2,  n_agem, t - 1] *
    #                                          sn_sus[2, n_agem, t - 1] *
    #                                          psi[k, 2, n_agem, t - 1] 

    #     ### Male: set projection into population model matrix
    #     for (a in 2:n_agem) {
    #       pop_inf[k, 2, a, t] <- pop_inf_proj[k, 2, (a - 1), t]
    #     }

    #     #Male: fawn class
    #     #there are no infected fawns at birth
    #     pop_inf[k, 2, 1, t] <- 0

    #   }#end t
    # }#end study_area

    # ######################################################################
    # ### Observation Model
    # ######################################################################

    # tau_obs ~ dgamma(1, 1)

    # # for (i in 1:n_sex) {
    # #   tau_obs[i] ~ dgamma(1, 1)
    # # }#end i

    # for(k in 1:n_study_area) {

    #   # for (i in 1:n_sex) {
    #   #   tau_obs[k, i] ~ dgamma(1, 1)
    #   # }#end i

    #   for (t in 1:n_year) {
    #     for (a in 1:n_agef) {
    #       harv_pop[k, 1, a, t] <- (pop_inf[k, 1, a, t] *
    #                               (1 - sh_inf[1, a, t]) +
    #                               pop_sus[k, 1, a, t] *
    #                               (1 - sh_sus[1, a, t]) *
    #                               (1 - psi_hat[k, 1, a, t]) +
    #                               pop_sus[k, 1, a, t] *
    #                               (1 - sh_inf[1, a, t]) *
    #                               psi_hat[k, 1, a, t]) *
    #                               report[t]
    #     }
    #     for (a in 1:n_agem) {
    #       harv_pop[k, 2, a, t] <- (pop_inf[k, 2, a, t] *
    #                               (1 - sh_inf[2, a, t]) +
    #                               pop_sus[k, 2, a, t] *
    #                               (1 - sh_sus[2, a, t]) *
    #                               (1 - psi_hat[k, 2, a, t]) +
    #                               pop_sus[k, 2, a, t] *
    #                               (1 - sh_inf[2, a, t]) *
    #                               psi_hat[k, 2, a, t]) *
    #                               report[t]
    #     }

    #     #Total Antlerless Harvest
    #     #adding in male fawns
    #     mu_obs[k, 1, t] <- (sum(harv_pop[k, 1, 1:n_agef, t]) +
    #                         harv_pop[k, 2, 1, t]) # eab_antlerless[t] *

    #     #Total Antlered Harvest
    #     #excludes male fawns 
    #     mu_obs[k, 2, t] <- sum(harv_pop[k, 2, 2:n_agem, t]) #eab_antlered[t] *

    #     ###################################
    #     #Likelihood for overall total
    #     ###################################

    #     for (j in 1:n_sex) {
    #       lobs[k, j, t] ~ dnorm(log(mu_obs[k, j, t]), tau_obs)# tau_obs[j]
    #     }#end i

    #     ###################################
    #     #Likelihood for overall total
    #     ###################################

    #     ### parameters for likelihood harvest data by antlerless group
    #     ### proportion of each age class
    #     ### Antlerless deer
    #     p_less[k, 1, t] <- harv_pop[k, 1, 1, t] / mu_obs[k, 1, t]#f
    #     p_less[k, 2, t] <- harv_pop[k, 1, 2, t] / mu_obs[k, 1, t]#1
    #     p_less[k, 3, t] <- harv_pop[k, 1, 3, t] / mu_obs[k, 1, t]#2
    #     p_less[k, 4, t] <- harv_pop[k, 1, 4, t] / mu_obs[k, 1, t]#3
    #     p_less[k, 5, t] <- sum(harv_pop[k, 1, 5:6, t]) / mu_obs[k, 1, t]#4-5
    #     p_less[k, 6, t] <- sum(harv_pop[k, 1, 7:9, t]) / mu_obs[k, 1, t]#6-8
    #     p_less[k, 7, t] <- harv_pop[k, 1, 10, t] / mu_obs[k, 1, t]#9+
    #     p_less[k, 8, t] <- 1 - sum(p_less[k, 1:7, t])#male f, antlerless

    #     ### harvest data bt antlered group
    #     p_ant[k, 1, t] <- harv_pop[k, 2, 2, t] / mu_obs[k, 2, t]#1
    #     p_ant[k, 2, t] <- harv_pop[k, 2, 3, t] / mu_obs[k, 2, t]#2
    #     p_ant[k, 3, t] <- harv_pop[k, 2, 4, t] / mu_obs[k, 2, t]#3
    #     p_ant[k, 4, t] <- sum(harv_pop[k, 2, 5:6, t]) / mu_obs[k, 2, t]#4-5
    #     p_ant[k, 5, t] <- 1 - sum(p_ant[k, 1:4, t]) #6+

    #   }# end t

    #   for(t in 1:n_year) {
    #       #antlerless, male fawns
    #       Cage_less[k, 1:(n_ageclassf + 1),t] ~ dmulti(prob = p_less[k,1:(n_ageclassf + 1),t],
    #                             size = sizeCage_f[k, t])
    #       #antlered
    #       Cage_ant[k, 1:(n_ageclassm - 1),t] ~ dmulti(prob = p_ant[k, 1:(n_ageclassm - 1),t],
    #                 size = sizeCage_m[k, t])
    #     }
    #   }# end t

})#end model statement