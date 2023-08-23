#######################################################################
###
###   Likelihoods for each kind of data listed below
###
#######################################################################

# d_fit_hunt
# d_fit_endlive
# d_fit_sus_cens_posttest
# d_fit_sus_cens_postno
# d_fit_idead
# d_fit_sus_mort_posttest
# d_fit_sus_mort_postno
# d_fit_rec_neg_mort
# d_fit_rec_neg_cens
# d_fit_rec_pos_mort
# d_fit_rec_pos_cens
# d_fit_icap_cens
# d_fit_icap_mort

#######################################################################
###
###   User defined distribution for likelihood for
###  harvested deer, revised for multiple deer
###
###   d_fit_hunt
###
#######################################################################

dFOIhunt <- nimble::nimbleFunction( # nolint
    run = function( # nolint
                   ### argument type declarations
                   x = integer(0),
                   n_cases = double(1),
                   n_samples = integer(0), # number of samples in dataset
                   a = double(1), # age (weeks) at harvest
                   sex = double(1),
                   age2date = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double(0)) {

        # start the loop through individuals
        sumllik <- 0
        for (i in 1:n_samples) {

            # intitialize scalars

            #################################################
            ### loop over ages and accumulate sums over 1:a-1
            ### have to loop separately for lam_inf
            #################################################

            if (sex[i] == 0) { # age loops for females
                lik_foi <- 0
                lam_foij <- 0
                for (j in 1:(a[i] - 1)) {
                    # sum up foi hazard from 1  to j
                    lam_foij <- exp(space[sect[i]] +
                        f_age_foi[age_lookup_f[j]] +
                        f_period_foi[period_lookup_foi[age2date[i] + j]])
                    # sum up like_temp (no sus hazard when j=1)
                    lik_foi <- lik_foi + lam_foij 
                }
                p <- 1 - exp(-lik_foi)
                lik_temp <- dbinom(x,1,p,log=TRUE)
            } else { # age loops for males
                lik_foi <- 0
                lam_foij <- 0
                for (j in 1:(a[i] - 1)) {
                    # sum up foi hazard from 1  to j
                    lam_foij <- exp(space[sect[i]] +
                        m_age_foi[age_lookup_f[j]] +
                        m_period_foi[period_lookup_foi[age2date[i] + j]])
                    # sum up like_temp (no sus hazard when j=1)
                    lik_foi <- lik_foi + lam_foij 
                }
                p <- 1 - exp(-lik_foi)
                lik_temp <- dbinom(x,1,p,log=TRUE)
            } # end if(sex)

            #######################################
            ### accumulate the joint likelihood
            #######################################
            # if(is.na(lik_temp)){stop("ack")}
            sumllik <- sumllik + lik_temp * n_cases[i]
        } # end the loop for individual i
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        }
    }
)

nimble::registerDistributions(list(
    dFOIhunt = list(
        BUGSdist = "dFOIhunt(n_cases,n_samples,a,sex,age2date,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_cases = double(1)",
            "n_samples = integer(0)",
            "a = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "f_period_foi = double(1)",
            "m_period_foi = double(1)",
            "period_lookup_foi = double(1)",
            "space = double(1)",
            "sect = double(1)",
            "returnType = double(0)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dFOIhunt",
    dFOIhunt,
    envir = .GlobalEnv
)

# starttime <- Sys.time()
# test <- dFOIhunt(
#     x = d_fit_hunt$teststatus,
#     n_cases = d_fit_hunt$n_cases,
#     n_samples = nrow(d_fit_hunt),
#     a = d_fit_hunt$ageweeks,
#     sex = d_fit_hunt$sex,
#     age2date = d_fit_hunt$birthweek - 1,
#     f_age_foi = f_age_foi,
#     m_age_foi = m_age_foi,
#     age_lookup_f = age_lookup_f,
#     age_lookup_m = age_lookup_m,
#     period_lookup_foi = period_lookup_foi_hunt,
#     f_period_foi = f_period_foi,
#     m_period_foi = m_period_foi,
#     space = c(0, -.55),
#     sect = d_fit_hunt$ew,
#     log = TRUE
# )
# (end <- Sys.time() - starttime)
# test
#     x = d_fit_hunt$teststatus
#     n_cases = d_fit_hunt$n_cases
#     n_samples = nrow(d_fit_hunt)
#     a = d_fit_hunt$ageweeks
#     sex = d_fit_hunt$sex
#     age2date = d_fit_hunt$birthweek - 1
#     f_age_foi = f_age_foi
#     m_age_foi = m_age_foi
#     age_lookup_f = age_lookup_f
#     age_lookup_m = age_lookup_m
#     period_lookup_foi = period_lookup_foi
#     f_period_foi = f_period_foi
#     m_period_foi = m_period_foi
#     space = c(0, -.55)
#     sect = d_fit_hunt$ew


########################################################################
###
###   User defined distribution for FOI from collared deer 
###   with antemortem tests (+) at capture
### 
########################################################################

dFOIicap <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x=double(0),
        left = double(0),
        sex = double(0),
        age2date = double(0),
        f_age = double(1),
        m_age = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup = double(1),
        f_period = double(1),
        m_period = double(1),
        space = double(0),
        log = double()
        ) {

    logL<- 0 #intialize log-likelihood
    gam <- nimNumeric(left - 1)

    for (k in 1:(left - 1)) {
        if (sex == 0) {
            gam[k] <- space +
                      f_age[age_lookup_f[k]] +
                      f_period[period_lookup[k + age2date - 1]]
        } else {
            gam[k] <- space +
                      m_age[age_lookup_m[k]] +
                      m_period[period_lookup[k + age2date - 1]]
        }
    }
    #total probability of getting infected
    p <- 1 - exp(-sum(exp(gam[1:(left-1)])))
    logL <- dbinom(x,1,p,log=TRUE)

    returnType(double(0))
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dFOIicap = list(
        BUGSdist = 'dFOIicap(left,sex,age2date,f_age,m_age,age_lookup_f,age_lookup_m,f_period,m_period,period_lookup,space)',
        types = c('p = double(0)',
                  'left=double(0)',
                  'sex=double(0)',
                  'age2date=double(0)',
                  'f_age=double(1)',
                  'm_age=double(1)',
                  'age_lookup_f=double(1)',
                  'age_lookup_m=double(1)',
                  'f_period=double(1)',
                  'm_period=double(1)',
                  'period_lookup=double(1)',
                  'space=double(0)'
                  ),
        discrete = TRUE
    )
))

assign('dFOIicap', dFOIicap, envir = .GlobalEnv)


#############################################################
###
###   User defined distribution for FOI
###
############################################################


dFOIcollar <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x=double(),
        left = double(0),
        right = double(0),
        sex = double(0),
        age2date = double(0),
        f_age = double(1),
        m_age = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period=double(1),
        m_period=double(1),
        space = double(0),
        log = double()
        ) {

    logL<-0 #intialize log-likelihood
    gam <-nimNumeric(right)
    for (k in left:(right-1)) {
        if(sex == 0){
            gam[k] <- space +
                      f_age[age_lookup_f[k]] +
                      f_period[period_lookup[k + age2date]]
        } else {
             gam[k] <- space +
                       m_age[age_lookup_m[k]] +
                       m_period[period_lookup[ k + age2date]]
        }
    }
    #total probability of getting infected
    p <- 1 - exp(-sum(exp(gam[left:(right - 1)])))
    logL <- dbinom(x, 1, p, log = TRUE)
    if(is.na(logL)){stop("ack")}    
    returnType(double(0))
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dFOIcollar = list(
        BUGSdist = 'dFOIcollar(left,right,age2date,sex,f_age,m_age,age_lookup_f,age_lookup_m,f_period,m_period,period_lookup,space)',
        types = c('p = double(0)',
                  'gam = double(1)',
                  'left=double(0)',
                  'right=double(0)',
                  'age2date=double(0)',
                  'sex=double(0)',
                  'f_age=double(1)',
                  'm_age=double(1)',
                  'age_lookup_f=double(1)',
                  'age_lookup_m=double(1)',
                  'f_period=double(1)',
                  'm_period=double(1)',
                  'period_lookup=double(1)',
                  'space=double(0)'
                  ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign('dFOIcollar', dFOIcollar, envir = .GlobalEnv)

i=70
        x =  d_fit_sus_foi$teststatus[i]
        left = d_fit_sus_foi$left[i]
        right = d_fit_sus_foi$right[i]
        sex =  d_fit_sus_foi$sex[i]
        age2date = d_fit_sus_foi$age2date[i]


# test=c()
# for(i in 1:n_fit_sus_foi){
# test[i] <- dFOIcollar(
#         x =  d_fit_sus_foi$teststatus[i],
#         left = d_fit_sus_foi$left[i],
#         right = d_fit_sus_foi$right[i],
#         sex =  d_fit_sus_foi$sex[i],
#         age2date = d_fit_sus_foi$age2date[i],
#         m_age = c(rnorm(1, -6, sd = .1),
#                   rnorm(1, -5.5, sd = .1),
#                   rnorm(1, -5, sd = .1),
#                   rnorm(1, -5.5, sd = .1),
#                   rnorm(1, -6, sd = .1),
#                   rnorm(1, -7.2, sd = .1)) - 1.5,
#         f_age = c(rnorm(1, -6, sd = .1),
#                   rnorm(1, -5.5, sd = .1),
#                   rnorm(1, -6, sd = .1),
#                   rnorm(1, -6.5, sd = .1),
#                   rnorm(1, -6.8, sd = .1),
#                   rnorm(1, -7.2, sd = .1),
#                   rnorm(1, -8, sd = .1)) - 1.5,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup=period_lookup_foi,
#         f_period=rep(0, n_year),
#         m_period=rep(0, n_year),
#         space = -.55,
#         log = TRUE
#         ) 
# }
# test
# i=1
# x = d_fit_recap_foi$teststatus[i] 
#       left = d_fit_recap_foi$left[i]
#       right = d_fit_recap_foi$right[i]
#       sex = d_fit_recap_foi$sex[i]
#       age2date = d_fit_recap_foi$age2date[i]

#   for (i in 1:n_fit_recap_foi){
#     dFOIcollar(x = d_fit_recap_foi$teststatus[i],
#       left = d_fit_recap_foi$left[i],
#       right = d_fit_recap_foi$right[i],
#       sex = d_fit_recap_foi$sex[i],
#       age2date = d_fit_recap_foi$age2date[i],
#       f_age = f_age_foi[1:n_ageclassf],
#       m_age = m_age_foi[1:n_ageclassm],
#       age_lookup_f = age_lookup_f[1:n_age_lookup_f],
#       age_lookup_m = age_lookup_m[1:n_age_lookup_m],
#       period_lookup = period_lookup_foi[1:nT_period_overall],
#       f_period = f_period_foi[1:n_year],
#       m_period = m_period_foi[1:n_year],
#       space = -.55)
#   }



##################################################################
###
###  User Defined Distribution
###  Age-Period Survival Model
###
##################################################################

dSurvival <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        left = double(0),
        right = double(0),
        sex = double(0),
        age2date = double(0),
        age_effect = double(1),
        period_effect = double(1),
        nT_age = double(0),
        beta0 = double(0),
        beta_male = double(0),
        log = double()
        ) {
    
    logL<-0 #intialize log-likelihood
   
    UCH <-nimNumeric(nT_age)

    for (k in left:(right - 1)) {
        UCH[k] <- exp(beta0 + 
                      age_effect[k] +
                      period_effect[k + age2date] +
                      beta_male * sex)
    }
    # total prob of surviving
    p <- exp(-sum(UCH[left:(right - 1)]))
    logL <- dbinom(x, 1, p, log = TRUE)

    returnType(double())
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSurvival = list(
        BUGSdist = 'dSurvival(left,right,sex,age2date,age_effect,period_effect,nT_age,beta0,beta_male)',
        types = c('left = double(0)',
                  'right = double(0)',
                  'sex = double(0)',
                  'age2date = double(0)',
                  'age_effect = double(1)',
                  'period_effect = double(1)',
                  'nT_age = double(0)',
                  'beta0 = double(0)',
                  'beta_male = double(0)'
                  ),
        discrete = TRUE
    )
))
 
## for a user-defined distribution
assign('dSurvival', dSurvival, envir = .GlobalEnv)

# test <- rep(NA,n_fit_sus)
#   for (i in 1:n_fit_sus){
#     test[i]  <- dSurvival(d_fit_sus$surv_censor[i],
#       left = d_fit_sus$left[i],
#       right = d_fit_sus$right[i],
#       sex = d_fit_sus$sex[i],
#       age2date = d_fit_sus$age2date[i],
#       age_effect = age_effect_survival_test[1:nT_age_surv],
#       period_effect = period_effect_survival_test[1:nT_period_overall],
#       nT_age = nT_age_surv,
#       beta0 = beta0_survival_sus,
#       beta_male = beta_male,
#       log=TRUE)
#   }
# test

##################################################################
###
###  User Defined Distribution Age-Period Survival Model
###  for infected at mortality, must draw when infected
###
###  d_fit_idead
###
##################################################################

dSurvival_idead <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e_age = double(0),
        r_age = double(0),
        s_age = double(0),
        e_period = double(0),
        s_period = double(0),
        sex = double(0),
        age2date = double(0),
        age_effect_survival = double(1),
        period_effect_survival = double(1),
        nT_age_surv = double(0),
        beta0_survival_inf = double(0),
        beta0_survival_sus = double(0),
        beta_male = double(0),
        f_age_foi = double(1),
        m_age_foi = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        nT_period_overall = double(0),
        period_lookup_foi = double(1),
        space = double(0),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        log = double()
        ) {

    ntemp <- s_age - e_age
    prob <- nimNumeric(ntemp)
    prob_out <- nimNumeric(ntemp)
    hazard <- nimNumeric(ntemp)
    if(sex == 0){

        hazard[1:ntemp] <- exp(f_age_foi[age_lookup_f[e_age:(s_age - 1)]] +
            f_period_foi[period_lookup_foi[(e_period):
                                           (s_period - 1)]] +
            space)

        prob[1] <- (1 - exp(-hazard[1]))
        for (j in 2:ntemp) {
            prob[j] <- (1 - exp(-hazard[j])) * exp(sum(-hazard[1:(j - 1)]))
        }
        prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
    } else {
        hazard[1:ntemp] <- exp(m_age_foi[age_lookup_m[e_age:(s_age - 1)]] +
            m_period_foi[period_lookup_foi[(e_period):
                                           (s_period - 1)]] +
            space)

        prob[1] <- (1 - exp(-hazard[1]))
        for (j in 2:ntemp) {
            prob[j] <- (1 - exp(-hazard[j])) * exp(sum(-hazard[1:(j - 1)]))
        }
        prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
    }
    age_add <- rcat(n = 1, prob[1:ntemp])

    ###############################################
    ### Survival likelihood  - Susceptible portion 
    ###############################################

    UCH_sus <-nimNumeric(nT_age_surv)
    for (k in e_age:(e_age + age_add - 1)) {
        UCH_sus[k] <- exp(beta0_survival_sus + 
                      age_effect_survival[k] +
                      period_effect_survival[k + age2date] +
                      beta_male * sex)
    }
    # total prob of surviving
    p <- exp(-sum(UCH_sus[e_age:(e_age + age_add - 1)]))
    lik_sus <- dbinom(1, 1, p, log = TRUE)

    ###############################################
    ### Survival likelihood  - Infected portion 
    ###############################################

    UCH_inf <-nimNumeric(nT_age_surv)
    for (k in (e_age + age_add):(s_age - 1)) {
        UCH_inf[k] <- exp(beta0_survival_inf + 
                      age_effect_survival[k] +
                      period_effect_survival[k + age2date] +
                      beta_male * sex)
    }
    # total prob of surviving
    p_inf_alive <- exp(-sum(UCH_inf[(e_age + age_add):(r_age - 1)]))
    lik_inf_alive <- dbinom(1, 1, p_inf_alive, log = TRUE)
    p_inf_dead <- 1 - exp(-sum(UCH_inf[(r_age):(s_age - 1)]))
    lik_inf_dead <- dbinom(0, 1, p_inf_dead, log = TRUE)

    #total log likelihood
    logL <- lik_sus + lik_inf_alive + lik_inf_dead

    returnType(double())
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSurvival_idead = list(
        BUGSdist = 'dSurvival_idead(e_age,
        r_age,
        s_age,
        e_period,
        s_period,
        sex,
        age2date,
        age_effect_survival,
        period_effect_survival,
        nT_age_surv,
        beta0_survival_inf,
        beta0_survival_sus,
        beta_male,
        f_age_foi,
        m_age_foi,
        f_period_foi,
        m_period_foi,
        nT_period_overall,
        period_lookup_foi,
        space,
        age_lookup_f,
        age_lookup_m
        )',
        types = c(
        "e_age = double(0)",
        "r_age = double(0)",
        "s_age = double(0)",
        "e_period = double(0)",
        "s_period = double(0)",
        "sex = double(0)",
        "age2date = double(0)",
        "age_effect_survival = double(1)",
        "period_effect_survival = double(1)",
        "nT_age_surv = double(0)",
        "beta0_survival_inf = double(0)",
        "beta0_survival_sus = double(0)",
        "beta_male = double(0)",
        "f_age_foi = double(1)",
        "m_age_foi = double(1)",
        "f_period_foi = double(1)",
        "m_period_foi = double(1)",
        "nT_period_overall = double(0)",
        "period_lookup_foi = double(1)",
        "space = double(0)",
        "age_lookup_f = double(1)",
        "age_lookup_m = double(1)"
                  ),
        discrete = TRUE
    )
))

## for a user-defined distribution
assign('dSurvival_idead', dSurvival_idead, envir = .GlobalEnv)


# i=1
# space_test <- c(0,-.5)
# starttime <- Sys.time()
# test <- dSurvival_idead(1,
#         e_age = d_fit_idead$left_age_e[i],
#         r_age = d_fit_idead$right_age_r[i],
#         s_age = d_fit_idead$right_age_s[i],
#         e_period = d_fit_idead$left_period_e[i],
#         s_period = d_fit_idead$right_period_s[i],
#         sex = d_fit_idead$sex[i],
#         age2date = d_fit_idead$age2date[i],
#         age_effect_survival = age_effect_survival_test,
#         period_effect_survival = period_effect_survival_test,
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi,
#         space = space_test[d_fit_idead$study_area[i]],
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         log = TRUE
# )
# (endtime <- Sys.time() - starttime)
# test


##################################################################
###
###  User Defined Distribution Age-Period Survival Model
###  for infected at recapture and was then right censored
###  must draw when infected
###
###  d_fit_rec_pos_cens
###
##################################################################

dSurvival_rec_pos_cens <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e_age = double(0),
        r_age = double(0),
        recap_age = double(0),
        e_period = double(0),
        recap_period = double(0),
        sex = double(0),
        age2date = double(0),
        age_effect_survival = double(1),
        period_effect_survival = double(1),
        nT_age_surv = double(0),
        beta0_survival_inf = double(0),
        beta0_survival_sus = double(0),
        beta_male = double(0),
        f_age_foi = double(1),
        m_age_foi = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        nT_period_overall = double(0),
        period_lookup_foi = double(1),
        space = double(0),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        log = double()
        ) {

    ntemp <- recap_age - e_age
    prob <- nimNumeric(ntemp)
    prob_out <- nimNumeric(ntemp)
    hazard <- nimNumeric(ntemp)
    # if(sex == 0){
        #  hazard[1:ntemp] <- -exp(f_age_foi[age_lookup_f[e_age:(recap_age - 1)]] +
        #     f_period_foi[period_lookup_foi[(e_period):
        #                                    (recap_period - 1)]] +
        #     space)
        
    #     prob[1] <- (1 - exp(hazard[1]))
    #     for (j in 2:ntemp) {
    #         prob[j] <- (1 - exp(hazard[j])) * exp(sum(hazard[1:(j - 1)]))
    #     }
    #     prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
    # } else {
        hazard[1:ntemp] <- exp(m_age_foi[age_lookup_m[e_age:(recap_age - 1)]] +
            m_period_foi[period_lookup_foi[(e_period):
                                           (recap_period - 1)]] +
            space)
        prob[1] <- (1 - exp(-hazard[1]))
        for (j in 2:ntemp) {
            prob[j] <- (1 - exp(-hazard[j])) * exp(sum(-hazard[1:(j - 1)]))
        }
        prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
    # }
    age_add <- rcat(n = 1, prob[1:ntemp])

    ###############################################
    ### Survival likelihood  - Susceptible portion 
    ###############################################

    UCH_sus <-nimNumeric(nT_age_surv)
    for (k in e_age:(e_age + age_add -1)) {
        UCH_sus[k] <- exp(beta0_survival_sus + 
                      age_effect_survival[k] +
                      period_effect_survival[k + age2date] +
                      beta_male * sex)
    }
    # total prob of surviving
    p <- exp(-sum(UCH_sus[e_age:(e_age + age_add - 1)]))
    lik_sus <- dbinom(1, 1, p, log = TRUE)

    ###############################################
    ### Survival likelihood  - Infected portion 
    ###############################################

    UCH_inf <-nimNumeric(nT_age_surv)
    for (k in (e_age + age_add):(r_age - 1)) {
        UCH_inf[k] <- exp(beta0_survival_inf + 
                      age_effect_survival[k] +
                      period_effect_survival[k + age2date] +
                      beta_male * sex)
    }
    # total prob of surviving
    p_inf_alive <- exp(-sum(UCH_inf[(e_age + age_add):(r_age - 1)]))
    lik_inf_alive <- dbinom(1, 1, p_inf_alive, log = TRUE)
   
    #total log likelihood
    logL <- lik_sus + lik_inf_alive

    returnType(double())
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSurvival_rec_pos_cens = list(
        BUGSdist = 'dSurvival_rec_pos_cens(e_age,
        r_age,
        recap_age,
        e_period,
        recap_period,
        sex,
        age2date,
        age_effect_survival,
        period_effect_survival,
        nT_age_surv,
        beta0_survival_inf,
        beta0_survival_sus,
        beta_male,
        f_age_foi,
        m_age_foi,
        f_period_foi,
        m_period_foi,
        nT_period_overall,
        period_lookup_foi,
        space,
        age_lookup_f,
        age_lookup_m
        )',
        types = c(
        "e_age = double(0)",
        "r_age = double(0)",
        "recap_age = double(0)",
        "e_period = double(0)",
        "recap_period = double(0)",
        "sex = double(0)",
        "age2date = double(0)",
        "age_effect_survival = double(1)",
        "period_effect_survival = double(1)",
        "nT_age_surv = double(0)",
        "beta0_survival_inf = double(0)",
        "beta0_survival_sus = double(0)",
        "beta_male = double(0)",
        "f_age_foi = double(1)",
        "m_age_foi = double(1)",
        "f_period_foi = double(1)",
        "m_period_foi = double(1)",
        "nT_period_overall = double(0)",
        "period_lookup_foi = double(1)",
        "space = double(0)",
        "age_lookup_f = double(1)",
        "age_lookup_m = double(1)"
                  ),
        discrete = TRUE
    )
))

## for a user-defined distribution
assign('dSurvival_rec_pos_cens', dSurvival_rec_pos_cens, envir = .GlobalEnv)

# i=1
# space_test <- c(0,-.5)
# starttime <- Sys.time()
# test <- dSurvival_rec_pos_cens(x = 1,
#         e_age = d_fit_rec_pos_cens$left_age_e[i],
#         r_age = d_fit_rec_pos_cens$right_age_r[i],
#         recap_age = d_fit_rec_pos_cens$ageweek_recap[i],
#         e_period = d_fit_rec_pos_cens$left_period_e[i],
#         recap_period = d_fit_rec_pos_cens$periodweek_recap[i],
#         sex = d_fit_rec_pos_cens$sex[i],
#         age2date = d_fit_rec_pos_cens$age2date[i],
#         age_effect_survival = age_effect_survival_test,
#         period_effect_survival = period_effect_survival_test,
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi,
#         space = space_test[d_fit_rec_pos_cens$study_area[i]],
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         log = TRUE
# )
# (endtime <- Sys.time() - starttime)
# test

##################################################################
###
###  User Defined Distribution Age-Period Survival Model
###  for infected at recapture and was then right censored
###  must draw when infected
###
###  d_fit_rec_pos_mort
###
##################################################################

dSurvival_rec_pos_mort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e_age = double(0),
        r_age = double(0),
        s_age = double(0),
        recap_age = double(0),
        e_period = double(0),
        recap_period = double(0),
        sex = double(0),
        age2date = double(0),
        age_effect_survival = double(1),
        period_effect_survival = double(1),
        nT_age_surv = double(0),
        beta0_survival_inf = double(0),
        beta0_survival_sus = double(0),
        beta_male = double(0),
        f_age_foi = double(1),
        m_age_foi = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        nT_period_overall = double(0),
        period_lookup_foi = double(1),
        space = double(0),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        log = double()
        ) {

    ntemp <- recap_age - e_age
    prob <- nimNumeric(ntemp)
    prob_out <- nimNumeric(ntemp)
    hazard <- nimNumeric(ntemp)
    if(sex == 0){
        hazard[1:ntemp] <- exp(f_age_foi[age_lookup_f[e_age:(recap_age - 1)]] +
            f_period_foi[period_lookup_foi[(e_period):
                                           (recap_period - 1)]] +
            space)
        prob[1] <- (1 - exp(hazard[1]))
        for (j in 2:ntemp) {
            prob[j] <- (1 - exp(hazard[j])) * exp(sum(hazard[1:(j - 1)]))
        }
        prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
    } else {

        hazard[1:ntemp] <- exp(m_age_foi[age_lookup_m[e_age:(recap_age - 1)]] +
            m_period_foi[period_lookup_foi[(e_period):
                                           (recap_period - 1)]] +
            space)
        prob[1] <- (1 - exp(-hazard[1]))
        for (j in 2:ntemp) {
            prob[j] <- (1 - exp(-hazard[j])) * exp(sum(-hazard[1:(j - 1)]))
        }
        prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
    }
    age_add <- rcat(n = 1, prob_out[1:ntemp])

    ###############################################
    ### Survival likelihood  - Susceptible portion 
    ###############################################

    UCH_sus <-nimNumeric(nT_age_surv)
    for (k in e_age:(e_age + age_add - 1)) {
        UCH_sus[k] <- exp(beta0_survival_sus + 
                      age_effect_survival[k] +
                      period_effect_survival[k + age2date] +
                      beta_male * sex)
    }
    # total prob of surviving
    p <- exp(-sum(UCH_sus[e_age:(e_age + age_add - 1)]))
    lik_sus <- dbinom(1, 1, p, log = TRUE)

    ###############################################
    ### Survival likelihood  - Infected portion 
    ###############################################

    UCH_inf <-nimNumeric(nT_age_surv)
    for (k in (e_age + age_add):(s_age - 1)) {
        UCH_inf[k] <- exp(beta0_survival_inf + 
                      age_effect_survival[k] +
                      period_effect_survival[k + age2date] +
                      beta_male * sex)
    }
    # total prob of surviving
    p_inf_alive <- exp(-sum(UCH_inf[(e_age + age_add):(r_age - 1)]))
    lik_inf_alive <- dbinom(1, 1, p_inf_alive, log = TRUE)
    p_inf_dead <- 1 - exp(-sum(UCH_inf[(r_age):(s_age - 1)]))
    lik_inf_dead <- dbinom(0, 1, p_inf_dead, log = TRUE)

    #total log likelihood
    logL <- lik_sus + lik_inf_alive + lik_inf_dead

    returnType(double())
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSurvival_rec_pos_mort = list(
        BUGSdist = 'dSurvival_rec_pos_mort(e_age,
        r_age,
        s_age,
        recap_age,
        e_period,
        recap_period,
        sex,
        age2date,
        age_effect_survival,
        period_effect_survival,
        nT_age_surv,
        beta0_survival_inf,
        beta0_survival_sus,
        beta_male,
        f_age_foi,
        m_age_foi,
        f_period_foi,
        m_period_foi,
        nT_period_overall,
        period_lookup_foi,
        space,
        age_lookup_f,
        age_lookup_m
        )',
        types = c(
        "e_age = double(0)",
        "r_age = double(0)",
        "s_age = double(0)",
        "recap_age = double(0)",
        "e_period = double(0)",
        "recap_period = double(0)",
        "sex = double(0)",
        "age2date = double(0)",
        "age_effect_survival = double(1)",
        "period_effect_survival = double(1)",
        "nT_age_surv = double(0)",
        "beta0_survival_inf = double(0)",
        "beta0_survival_sus = double(0)",
        "beta_male = double(0)",
        "f_age_foi = double(1)",
        "m_age_foi = double(1)",
        "f_period_foi = double(1)",
        "m_period_foi = double(1)",
        "nT_period_overall = double(0)",
        "period_lookup_foi = double(1)",
        "space = double(0)",
        "age_lookup_f = double(1)",
        "age_lookup_m = double(1)"
                  ),
        discrete = TRUE
    )
))

## for a user-defined distribution
assign('dSurvival_rec_pos_mort', dSurvival_rec_pos_mort, envir = .GlobalEnv)

# i=1
# space_test <- c(0,-.5)
# starttime <- Sys.time()
# test <- dSurvival_rec_pos_mort(x = 1,
#         e_age = d_fit_rec_pos_mort$left_age_e[i],
#         r_age = d_fit_rec_pos_mort$right_age_r[i],
#         s_age = d_fit_rec_pos_mort$right_age_s[i],
#         recap_age = d_fit_rec_pos_mort$ageweek_recap[i],
#         e_period = d_fit_rec_pos_mort$left_period_e[i],
#         recap_period = d_fit_rec_pos_mort$periodweek_recap[i],
#         sex = d_fit_rec_pos_mort$sex[i],
#         age2date = d_fit_rec_pos_mort$age2date[i],
#         age_effect_survival = age_effect_survival_test,
#         period_effect_survival = period_effect_survival_test,
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi,
#         space = space_test[d_fit_rec_pos_mort$study_area[i]],
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         log = TRUE
# )
# (endtime <- Sys.time() - starttime)
# test


# test <- c()
# starttime <- Sys.time()

# for(i in 1:n_fit_rec_pos_mort){
#     test[i] <- dSurvival_rec_pos_mort(x = 1,
#             e_age = d_fit_rec_pos_mort$left_age_e[i],
#             r_age = d_fit_rec_pos_mort$right_age_r[i],
#             s_age = d_fit_rec_pos_mort$right_age_s[i],
#             recap_age = d_fit_rec_pos_mort$ageweek_recap[i],
#             e_period = d_fit_rec_pos_mort$left_period_e[i],
#             recap_period = d_fit_rec_pos_mort$periodweek_recap[i],
#             sex = d_fit_rec_pos_mort$sex[i],
#             age2date = d_fit_rec_pos_mort$age2date[i],
#             age_effect_survival = age_effect_survival_test,
#             period_effect_survival = period_effect_survival_test,
#             nT_age_surv = nT_age_surv,
#             beta0_survival_inf = beta0_survival_inf,
#             beta0_survival_sus = beta0_survival_sus,
#             beta_male = beta_male,
#             f_age_foi = f_age_foi,
#             m_age_foi = m_age_foi,
#             f_period_foi = f_period_foi,
#             m_period_foi = m_period_foi,
#             nT_period_overall = nT_period_overall,
#             period_lookup_foi = period_lookup_foi,
#             space = space_test[d_fit_rec_pos_mort$study_area[i]],
#             age_lookup_f = age_lookup_f,
#             age_lookup_m = age_lookup_m,
#             log = TRUE
#     )
# }


# (endtime <- Sys.time() - starttime)
# test



###########################################################################
###
###  User Defined Distribution Age-Period Survival Model
###  no test after initial capture, first drawing if gets infected
###  then if it gets infected, accounting for the suscepting and infected
###  portions in the joint likelihood
###
###  sus_draw
###
############################################################################

dSurvival_sus_draw <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e_age = double(0),
        r_age = double(0),
        e_period = double(0),
        r_period = double(0),
        sex = double(0),
        age2date = double(0),
        age_effect_survival = double(1),
        period_effect_survival = double(1),
        nT_age_surv = double(0),
        beta0_survival_inf = double(0),
        beta0_survival_sus = double(0),
        beta_male = double(0),
        f_age_foi = double(1),
        m_age_foi = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        nT_period_overall = double(0),
        period_lookup_foi = double(1),
        space = double(0),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        log = double()
        ) {

    ntemp <- r_age - e_age
    hazard <- nimNumeric(ntemp)
    prob <- nimNumeric(ntemp)
    prob_out <- nimNumeric(ntemp)
    ################################
    ### First draw if gets infected
    ################################

    if(sex == 0){
        hazard[1:ntemp] <- exp(f_age_foi[age_lookup_f[e_age:(r_age - 1)]] +
            f_period_foi[period_lookup_foi[(e_period):
                         (r_period - 1)]] +
            space)

    } else {
        hazard[1:ntemp] <- exp(m_age_foi[age_lookup_m[e_age:(r_age - 1)]] +
            m_period_foi[period_lookup_foi[(e_period):
                         (r_period - 1)]] +
            space)
    }
    p_inf <- 1 - exp(-sum(hazard))
    gets_infected <- rbinom(1,1,p_inf)

    if(gets_infected == TRUE){

        prob[1] <- (1 - exp(-hazard[1]))
        for (j in 2:ntemp) {
            prob[j] <- (1 - exp(-hazard[j])) * exp(-sum(hazard[1:(j - 1)]))
        }
        prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
        age_add <- rcat(n = 1, prob[1:ntemp])

        ###############################################
        ### Survival likelihood  - Susceptible portion 
        ###############################################

        UCH_sus <-nimNumeric(nT_age_surv)
        for (k in e_age:(e_age + age_add - 1)) {
            UCH_sus[k] <- exp(beta0_survival_sus + 
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p <- exp(-sum(UCH_sus[e_age:(e_age + age_add - 1)]))
        lik_sus <- dbinom(1, 1, p, log = TRUE)

        ###############################################
        ### Survival likelihood  - Infected portion 
        ###############################################

        UCH_inf <-nimNumeric(nT_age_surv)
        for (k in (e_age + age_add):(r_age - 1)) {
            UCH_inf[k] <- exp(beta0_survival_inf + 
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p_inf_alive <- exp(-sum(UCH_inf[(e_age + age_add):(r_age - 1)]))
        lik_inf_alive <- dbinom(1, 1, p_inf_alive, log = TRUE)

        #total log likelihood
        logL <- lik_sus + lik_inf_alive

    } else {

        ###############################################
        ### Survival likelihood  - Susceptible portion 
        ###############################################

        UCH_sus <- nimNumeric(nT_age_surv)
        for (k in e_age:(r_age - 1)) {
            UCH_sus[k] <- exp(beta0_survival_sus +
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p <- exp(-sum(UCH_sus[e_age:(r_age - 1)]))
        logL <- dbinom(1, 1, p, log = TRUE)
    } #end else
    returnType(double())
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSurvival_sus_draw = list(
        BUGSdist = 'dSurvival_sus_draw(e_age,
        r_age,
        e_period,
        r_period,
        sex,
        age2date,
        age_effect_survival,
        period_effect_survival,
        nT_age_surv,
        beta0_survival_inf,
        beta0_survival_sus,
        beta_male,
        f_age_foi,
        m_age_foi,
        f_period_foi,
        m_period_foi,
        nT_period_overall,
        period_lookup_foi,
        space,
        age_lookup_f,
        age_lookup_m
        )',
        types = c(
        "e_age = double(0)",
        "r_age = double(0)",
        "e_period = double(0)",
        "r_period = double(0)",
        "sex = double(0)",
        "age2date = double(0)",
        "age_effect_survival = double(1)",
        "period_effect_survival = double(1)",
        "nT_age_surv = double(0)",
        "beta0_survival_inf = double(0)",
        "beta0_survival_sus = double(0)",
        "beta_male = double(0)",
        "f_age_foi = double(1)",
        "m_age_foi = double(1)",
        "f_period_foi = double(1)",
        "m_period_foi = double(1)",
        "nT_period_overall = double(0)",
        "period_lookup_foi = double(1)",
        "space = double(0)",
        "age_lookup_f = double(1)",
        "age_lookup_m = double(1)"
                  ),
        discrete = TRUE
    )
))

## for a user-defined distribution
assign('dSurvival_sus_draw', dSurvival_sus_draw, envir = .GlobalEnv)
# i=1
# space_test <- c(0,-.5)
# starttime <- Sys.time()
# test <- dSurvival_sus_draw(1,
#         e_age = d_fit_sus_draw$left_age_e[i],
#         r_age = d_fit_sus_draw$right_age_r[i],
#         e_period = d_fit_sus_draw$left_period_e[i],
#         r_period = d_fit_sus_draw$right_period_r[i],
#         sex = d_fit_sus_draw$sex[i],
#         age2date = d_fit_sus_draw$age2date[i],
#         age_effect_survival = age_effect_survival_test,
#         period_effect_survival = period_effect_survival_test,
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi,
#         space = space_test[d_fit_sus_draw$study_area[i]],
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         log = TRUE
# )
# (endtime <- Sys.time() - starttime)
# test



##################################################################
###
###  User Defined Distribution Age-Period Survival Model
###  no test after initial capture, first drawing if gets infected
###  then if it gets infected, accounting for the susceptible and infected
###  portions in the joint likelihood, these individuals died
###
###  d_fit_sus_mort_postno
###
##################################################################

dSurvival_sus_mort_postno <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e_age = double(0),
        r_age = double(0),
        s_age = double(0),
        e_period = double(0),
        s_period = double(0),
        sex = double(0),
        age2date = double(0),
        age_effect_survival = double(1),
        period_effect_survival = double(1),
        nT_age_surv = double(0),
        beta0_survival_inf = double(0),
        beta0_survival_sus = double(0),
        beta_male = double(0),
        f_age_foi = double(1),
        m_age_foi = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        nT_period_overall = double(0),
        period_lookup_foi = double(1),
        space = double(0),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        log = double()
        ) {

    ntemp <- s_age - e_age
    hazard <- nimNumeric(ntemp)
    prob <- nimNumeric(ntemp)
    prob_out <- nimNumeric(ntemp)

    ################################
    ### First draw if gets infected
    ################################

    if(sex == 0){
         hazard[1:ntemp] <- exp(f_age_foi[age_lookup_f[e_age:(s_age - 1)]] +
            f_period_foi[period_lookup_foi[(e_period):
                         (s_period - 1)]] +
            space)
    } else {
        hazard[1:ntemp] <- exp(m_age_foi[age_lookup_m[e_age:(s_age - 1)]] +
            m_period_foi[period_lookup_foi[(e_period):
                         (s_period - 1)]] +
            space)
    }
    p_inf <- 1 - exp(-sum(hazard))
    gets_infected <- rbinom(1,1,p_inf)

    if(gets_infected == TRUE){

        prob[1] <- (1 - exp(-hazard[1]))
        for (j in 2:ntemp) {
            prob[j] <- (1 - exp(-hazard[j])) * exp(-sum(hazard[1:(j - 1)]))
        }
        prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
        age_add <- rcat(n = 1, prob[1:ntemp])

        ###############################################
        ### Survival likelihood  - Susceptible portion 
        ###############################################

        UCH_sus <-nimNumeric(nT_age_surv)
        for (k in e_age:(e_age + age_add - 1)) {
            UCH_sus[k] <- exp(beta0_survival_sus + 
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p <- exp(-sum(UCH_sus[e_age:(e_age + age_add - 1)]))
        lik_sus <- dbinom(1, 1, p, log = TRUE)

        ###############################################
        ### Survival likelihood  - Infected portion 
        ###############################################

        UCH_inf <-nimNumeric(nT_age_surv)
        for (k in (e_age + age_add):(s_age - 1)) {
            UCH_inf[k] <- exp(beta0_survival_inf + 
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p_inf_alive <- exp(-sum(UCH_inf[(e_age + age_add):(r_age - 1)]))
        lik_inf_alive <- dbinom(1, 1, p_inf_alive, log = TRUE)
        p_inf_dead <- 1 - exp(-sum(UCH_inf[(r_age):(s_age - 1)]))
        lik_inf_dead <- dbinom(0, 1, p_inf_dead, log = TRUE)

        #total log likelihood
        logL <- lik_sus + lik_inf_alive + lik_inf_dead

    } else {

        ###############################################
        ### Survival likelihood  - Susceptible portion 
        ###############################################

        UCH_sus <-nimNumeric(nT_age_surv)
        for (k in e_age:(s_age - 1)) {
            UCH_sus[k] <- exp(beta0_survival_sus +
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p_sus_alive <- exp(-sum(UCH_sus[e_age:r_age]))
        lik_sus_alive <- dbinom(1, 1, p_sus_alive, log = TRUE)
        p_sus_dead <- 1 - exp(-sum(UCH_sus[(r_age):(s_age - 1)]))
        lik_sus_dead <- dbinom(0, 1, p_sus_dead, log = TRUE)

        logL <- lik_sus + lik_sus_alive + lik_sus_dead
    } #end else
    returnType(double())
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSurvival_sus_mort_postno = list(
        BUGSdist = 'dSurvival_sus_mort_postno(e_age,
        r_age,
        s_age,
        e_period,
        s_period,
        sex,
        age2date,
        age_effect_survival,
        period_effect_survival,
        nT_age_surv,
        beta0_survival_inf,
        beta0_survival_sus,
        beta_male,
        f_age_foi,
        m_age_foi,
        f_period_foi,
        m_period_foi,
        nT_period_overall,
        period_lookup_foi,
        space,
        age_lookup_f,
        age_lookup_m
        )',
        types = c(
        "e_age = double(0)",
        "r_age = double(0)",
        "e_period = double(0)",
        "r_period = double(0)",
        "sex = double(0)",
        "age2date = double(0)",
        "age_effect_survival = double(1)",
        "period_effect_survival = double(1)",
        "nT_age_surv = double(0)",
        "beta0_survival_inf = double(0)",
        "beta0_survival_sus = double(0)",
        "beta_male = double(0)",
        "f_age_foi = double(1)",
        "m_age_foi = double(1)",
        "f_period_foi = double(1)",
        "m_period_foi = double(1)",
        "nT_period_overall = double(0)",
        "period_lookup_foi = double(1)",
        "space = double(0)",
        "age_lookup_f = double(1)",
        "age_lookup_m = double(1)"
                  ),
        discrete = TRUE
    )
))

## for a user-defined distribution
assign('dSurvival_sus_mort_postno', dSurvival_sus_mort_postno, envir = .GlobalEnv)


# i=1
# space_test <- c(0,-.5)
# starttime <- Sys.time()
# test <- dSurvival_sus_mort_postno(1,
#         e_age = d_fit_sus_mort_postno$left_age_e[i],
#         r_age = d_fit_sus_mort_postno$right_age_r[i],
#         s_age = d_fit_sus_mort_postno$right_age_s[i],
#         e_period = d_fit_sus_mort_postno$left_period_e[i],
#         s_period = d_fit_sus_mort_postno$right_period_s[i],
#         sex = d_fit_sus_mort_postno$sex[i],
#         age2date = d_fit_sus_mort_postno$age2date[i],
#         age_effect_survival = age_effect_survival_test,
#         period_effect_survival = period_effect_survival_test,
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi,
#         space = space_test[d_fit_sus_mort_postno$study_area[i]],
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         log = TRUE
# )
# (endtime <- Sys.time() - starttime)
# test


##################################################################
###
###  User Defined Distribution Age-Period Survival Model
###  no test after initial capture, first drawing if gets infected
###  then if it gets infected, accounting for the suscepting and infected
###  portions in the joint likelihood
###
###  d_fit_rec_neg_cens_postno
###
##################################################################

dSurvival_rec_neg_cens_postno <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e_age = double(0),
        r_age = double(0),
        recap_age = double(0),
        e_period = double(0),
        r_period = double(0),
        recap_period = double(0),
        sex = double(0),
        age2date = double(0),
        age_effect_survival = double(1),
        period_effect_survival = double(1),
        nT_age_surv = double(0),
        beta0_survival_inf = double(0),
        beta0_survival_sus = double(0),
        beta_male = double(0),
        f_age_foi = double(1),
        m_age_foi = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        nT_period_overall = double(0),
        period_lookup_foi = double(1),
        space = double(0),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        log = double()
        ) {

    ntemp <- recap_age - e_age
    hazard <- nimNumeric(ntemp)
    prob <- nimNumeric(ntemp)
    prob_out <- nimNumeric(ntemp)

    ################################
    ### First draw if gets infected
    ################################

    if(sex == 0){
        hazard[1:ntemp] <- exp(f_age_foi[age_lookup_f[e_age:(recap_age - 1)]] +
            f_period_foi[period_lookup_foi[(e_period):
                         (recap_period - 1)]] +
            space)
    } else {
        hazard[1:ntemp] <- exp(m_age_foi[age_lookup_m[e_age:(recap_age - 1)]] +
            m_period_foi[period_lookup_foi[(e_period):
                         (recap_period - 1)]] +
            space)
    }
    p_inf <- 1 - exp(-sum(hazard))
    gets_infected <- rbinom(1,1,p_inf)

    if(gets_infected == TRUE){

        prob[1] <- (1 - exp(-hazard[1]))
        for (j in 2:ntemp) {
            prob[j] <- (1 - exp(-hazard[j])) * exp(-sum(hazard[1:(j - 1)]))
        }
        prob_out[1:ntemp] <- prob[1:ntemp] / sum(prob[1:ntemp])
        age_add <- rcat(n = 1, prob[1:ntemp])

        ###############################################
        ### Survival likelihood  - Susceptible portion 
        ###############################################

        UCH_sus <-nimNumeric(nT_age_surv)
        for (k in e_age:(e_age + age_add - 1)) {
            UCH_sus[k] <- exp(beta0_survival_sus + 
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p <- exp(-sum(UCH_sus[e_age:(e_age + age_add - 1)]))
        lik_sus <- dbinom(1, 1, p, log = TRUE)

        ###############################################
        ### Survival likelihood  - Infected portion 
        ###############################################

        UCH_inf <-nimNumeric(nT_age_surv)
        for (k in (e_age + age_add):(r_age - 1)) {
            UCH_inf[k] <- exp(beta0_survival_inf + 
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p_inf_alive <- exp(-sum(UCH_inf[(e_age + age_add):(r_age - 1)]))
        lik_inf_alive <- dbinom(1, 1, p_inf_alive, log = TRUE)

        #total log likelihood
        logL <- lik_sus + lik_inf_alive

    } else {

        ###############################################
        ### Survival likelihood  - Susceptible portion 
        ###############################################

        UCH_sus <-nimNumeric(nT_age_surv)
        for (k in e_age:(r_age - 1)) {
            UCH_sus[k] <- exp(beta0_survival_sus +
                        age_effect_survival[k] +
                        period_effect_survival[k + age2date] +
                        beta_male * sex)
        }
        # total prob of surviving
        p <- exp(-sum(UCH_sus[e_age:(r_age - 1)]))
        logL <- dbinom(1, 1, p, log = TRUE)
    } #end else
    returnType(double())
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSurvival_rec_neg_cens_postno = list(
        BUGSdist = 'dSurvival_rec_neg_cens_postno(e_age,
        r_age,
        recap_age,
        e_period,
        r_period,
        recap_period,
        sex,
        age2date,
        age_effect_survival,
        period_effect_survival,
        nT_age_surv,
        beta0_survival_inf,
        beta0_survival_sus,
        beta_male,
        f_age_foi,
        m_age_foi,
        f_period_foi,
        m_period_foi,
        nT_period_overall,
        period_lookup_foi,
        space,
        age_lookup_f,
        age_lookup_m
        )',
        types = c(
        "e_age = double(0)",
        "r_age = double(0)",
        "recap_age = double(0)",
        "e_period = double(0)",
        "r_period = double(0)",
        "recap_period = double(0)",
        "sex = double(0)",
        "age2date = double(0)",
        "age_effect_survival = double(1)",
        "period_effect_survival = double(1)",
        "nT_age_surv = double(0)",
        "beta0_survival_inf = double(0)",
        "beta0_survival_sus = double(0)",
        "beta_male = double(0)",
        "f_age_foi = double(1)",
        "m_age_foi = double(1)",
        "f_period_foi = double(1)",
        "m_period_foi = double(1)",
        "nT_period_overall = double(0)",
        "period_lookup_foi = double(1)",
        "space = double(0)",
        "age_lookup_f = double(1)",
        "age_lookup_m = double(1)"
                  ),
        discrete = TRUE
    )
))

## for a user-defined distribution
assign('dSurvival_rec_neg_cens_postno', dSurvival_rec_neg_cens_postno, envir = .GlobalEnv)

# i=2

#         e_age = d_fit_rec_neg_cens_postno$left_age_e[i]
#         r_age = d_fit_rec_neg_cens_postno$right_age_r[i]
#         recap_age = d_fit_rec_neg_cens_postno$ageweek_recap[i]
#         e_period = d_fit_rec_neg_cens_postno$left_period_e[i]
#         r_period = d_fit_rec_neg_cens_postno$right_period_r[i]
#         recap_period = d_fit_rec_neg_cens_postno$periodweek_recap[i]
#         sex = d_fit_rec_neg_cens_postno$sex[i]
#         age2date = d_fit_rec_neg_cens_postno$age2date[i]
#         age_effect_survival = age_effect_survival_test
#         period_effect_survival = period_effect_survival_test
#         nT_age_surv = nT_age_surv
#         beta0_survival_inf = beta0_survival_inf
#         beta0_survival_sus = beta0_survival_sus
#         beta_male = beta_male
#         f_age_foi = f_age_foi
#         m_age_foi = m_age_foi
#         f_period_foi = f_period_foi
#         m_period_foi = m_period_foi
#         nT_period_overall = nT_period_overall
#         period_lookup_foi = period_lookup_foi
#         space = space_test[d_fit_rec_neg_cens_postno$study_area[i]]
#         age_lookup_f = age_lookup_f
#         age_lookup_m = age_lookup_m


# i=1
# space_test <- c(0,-.5)
# starttime <- Sys.time()
# test <- dSurvival_rec_neg_cens_postno(1,
        # e_age = d_fit_rec_neg_cens_postno$left_age_e[i],
        # r_age = d_fit_rec_neg_cens_postno$right_age_r[i],
        # recap_age = d_fit_rec_neg_cens_postno$ageweek_recap[i],
        # e_period = d_fit_rec_neg_cens_postno$left_period_e[i],
        # r_period = d_fit_rec_neg_cens_postno$right_period_r[i],
        # recap_period = d_fit_rec_neg_cens_postno$periodweek_recap[i],
        # sex = d_fit_rec_neg_cens_postno$sex[i],
        # age2date = d_fit_rec_neg_cens_postno$age2date[i],
        # age_effect_survival = age_effect_survival_test,
        # period_effect_survival = period_effect_survival_test,
        # nT_age_surv = nT_age_surv,
        # beta0_survival_inf = beta0_survival_inf,
        # beta0_survival_sus = beta0_survival_sus,
        # beta_male = beta_male,
        # f_age_foi = f_age_foi,
        # m_age_foi = m_age_foi,
        # f_period_foi = f_period_foi,
        # m_period_foi = m_period_foi,
        # nT_period_overall = nT_period_overall,
        # period_lookup_foi = period_lookup_foi,
        # space = space_test[d_fit_rec_neg_cens_postno$study_area[i]],
        # age_lookup_f = age_lookup_f,
        # age_lookup_m = age_lookup_m,
#         log = TRUE
# )
# (endtime <- Sys.time() - starttime)
# test


# test <- c()
# starttime <- Sys.time()
# for(i in 1:n_fit_rec_neg_cens_postno){
#     test[i] <- dSurvival_rec_neg_cens_postno(1,
#         e_age = d_fit_rec_neg_cens_postno$left_age_e[i],
#         r_age = d_fit_rec_neg_cens_postno$right_age_r[i],
#         recap_age = d_fit_rec_neg_cens_postno$ageweek_recap[i],
#         e_period = d_fit_rec_neg_cens_postno$left_period_e[i],
#         r_period = d_fit_rec_neg_cens_postno$right_period_r[i],
#         recap_period = d_fit_rec_neg_cens_postno$periodweek_recap[i],
#         sex = d_fit_rec_neg_cens_postno$sex[i],
#         age2date = d_fit_rec_neg_cens_postno$age2date[i],
#         age_effect_survival = age_effect_survival_test,
#         period_effect_survival = period_effect_survival_test,
#         nT_age_surv = nT_age_surv,
#         beta0_survival_inf = beta0_survival_inf,
#         beta0_survival_sus = beta0_survival_sus,
#         beta_male = beta_male,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         nT_period_overall = nT_period_overall,
#         period_lookup_foi = period_lookup_foi,
#         space = space_test[d_fit_rec_neg_cens_postno$study_area[i]],
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         log = TRUE
# )

# }
# (endtime <- Sys.time() - starttime)
# test