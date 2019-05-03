library(data.table)
library(ggplot2)
library(rjags)
library(Rcpp)
library(gridExtra)

#source c++ VE to iVE conversion method
sourceCpp("./cumVEtoInsFast.cpp")

#### functions.r ####
#calculate empirical p-values and credible intervals from 2 matrixes (data_a and data_b) with the posterior samples that are to be compared
#cols argument can be used to specify the index of the column(s) to be compared
empiricalDifference <- function(
  data_a,
  data_b,
  cols = NULL
){
  if(is.null(cols)){
    cols <- c(1:ncol(data_a))
  }
  data_a <- data_a[, cols]
  data_b <- data_b[, cols]
  data_difference <- matrix(
    0,
    ncol=ncol(data_a),
    nrow=nrow(data_a)*nrow(data_b)
  )
  for(r in 1:nrow(data_a)){
    data_difference[c(1:nrow(data_a))+nrow(data_a)*(r-1), ] <- - sweep(data_b, 2, data_a[r,], "-")
  }
  data_difference_low <- sapply(
    c(1:ncol(data_difference)),
    function(x){
      quantile(
        data_difference[, x],
        probs=c(0.025)
      )
    }
  )
  data_difference_median <- sapply(
    c(1:ncol(data_difference)),
    function(x){
      quantile(
        data_difference[, x],
        probs=c(0.5)
      )
    }
  )
  data_difference_high <- sapply(
    c(1:ncol(data_difference)),
    function(x){
      quantile(
        data_difference[, x],
        probs=c(0.975)
      )
    }
  )
  data_pval <- matrix(
    0,
    ncol=ncol(data_a),
    nrow=nrow(data_b)
  )
  #count number of times data_a is higher than data_b
  for(r in 1:nrow(data_a)){
    data_pval[r, ] <- colSums(
      sweep(data_b, 2, data_a[r,], "<")
    )
  }
  data_pval <- 1 - (colSums(data_pval) / (nrow(data_a)*nrow(data_b)))
  return(
    list(
      "pvalue" = data_pval,
      "diff_low" = data_difference_low,
      "diff_high" = data_difference_high,
      "diff_median" = data_difference_median
    )
  )
}


#### prepare_extracted_data.r ####

vaccine_efficacy_data <- fread("rotavirus_vaccine_efficacy_extracted.csv")
vaccine_efficacy_data[,"study_id"] <- 0
vaccine_efficacy_data[,"follow_up_p"] <- 0
vaccine_efficacy_data[,"cum_cases_placebo"] <- 0
vaccine_efficacy_data[,"cum_cases_vaccine"] <- 0
vaccine_efficacy_data[,"persontime_placebo"] <- 0
vaccine_efficacy_data[,"persontime_vaccine"] <- 0
vaccine_efficacy_data[,"N_placebo_start"] <- 0
vaccine_efficacy_data[,"N_vaccine_start"] <- 0

#calculate follow-up time per-period in each unique study
#loop through studies within each stratum
for(m in unique(vaccine_efficacy_data[, mortality])){
  vaccine_efficacy_stratum <- vaccine_efficacy_data[mortality == m]
  vaccine_efficacy_stratum[,"study_id"] <- cumsum(
    vaccine_efficacy_stratum[,period] == "Period 1"
  )
  for(s in unique(vaccine_efficacy_stratum[, study_id])){
    vaccine_efficacy_study <- vaccine_efficacy_stratum[study_id == s]
    #loop through study-periods
    for(p in 1:nrow(vaccine_efficacy_study)){
      if(p == 1){
        vaccine_efficacy_study[p,"follow_up_p"] <- vaccine_efficacy_study[p,follow_up_end]
        vaccine_efficacy_study[p,"cum_cases_placebo"] <- vaccine_efficacy_study[p,cases_placebo]
        vaccine_efficacy_study[p,"cum_cases_vaccine"] <- vaccine_efficacy_study[p,cases_vaccine]
      } else {
        vaccine_efficacy_study[p,"follow_up_p"] <- vaccine_efficacy_study[p,follow_up_end]
        vaccine_efficacy_study[p,"cum_cases_placebo"] <- vaccine_efficacy_study[p-1,cum_cases_placebo] + vaccine_efficacy_study[p,cases_placebo]
        vaccine_efficacy_study[p,"cum_cases_vaccine"] <- vaccine_efficacy_study[p-1,cum_cases_vaccine] + vaccine_efficacy_study[p,cases_vaccine]
      }
    }
    vaccine_efficacy_study[,"N_placebo_start"] <- vaccine_efficacy_study[1,N_placebo]
    vaccine_efficacy_study[,"N_vaccine_start"] <- vaccine_efficacy_study[1,N_vaccine]
    
    #calculate follow-up time, use N at start in denominator, as studies do not consistently censor at follow-up
    vaccine_efficacy_study[,"persontime_placebo"] <- vaccine_efficacy_study[,N_placebo_start] * vaccine_efficacy_study[,follow_up_p]
    vaccine_efficacy_study[,"persontime_vaccine"] <- vaccine_efficacy_study[,N_vaccine_start] * vaccine_efficacy_study[,follow_up_p]
    
    vaccine_efficacy_stratum <- rbindlist(
      list(
        vaccine_efficacy_stratum[study_id != s],
        vaccine_efficacy_study
      )
    )
  }
  
  vaccine_efficacy_data <- rbindlist(
    list(
      vaccine_efficacy_data[mortality != m],
      vaccine_efficacy_stratum
    )
  )
}

#add 0.5 if no cases have been observed in vaccinated arm (cannot divide by 0), validated method
vaccine_efficacy_data[cum_cases_vaccine == 0, "N_placebo_start"] <- vaccine_efficacy_data[cum_cases_vaccine == 0, N_placebo_start] + 0.5
vaccine_efficacy_data[cum_cases_vaccine == 0, "N_vaccine_start"] <- vaccine_efficacy_data[cum_cases_vaccine == 0, N_vaccine_start] + 0.5
vaccine_efficacy_data[cum_cases_vaccine == 0, "cases_placebo"] <- vaccine_efficacy_data[cum_cases_vaccine == 0, cases_placebo] + 0.5
vaccine_efficacy_data[cum_cases_vaccine == 0, "cases_vaccine"] <- vaccine_efficacy_data[cum_cases_vaccine == 0, cases_vaccine] + 0.5
vaccine_efficacy_data[cum_cases_vaccine == 0, "cum_cases_placebo"] <- vaccine_efficacy_data[cum_cases_vaccine == 0, cum_cases_placebo] + 0.5
vaccine_efficacy_data[cum_cases_vaccine == 0, "cum_cases_vaccine"] <- vaccine_efficacy_data[cum_cases_vaccine == 0, cum_cases_vaccine] + 0.5

#calculate the relative risk (relative rate not possible with this data)
vaccine_efficacy_data[,"rr"] <- (
  vaccine_efficacy_data[,cum_cases_vaccine]/vaccine_efficacy_data[,N_vaccine_start]
)/(
  vaccine_efficacy_data[,cum_cases_placebo]/vaccine_efficacy_data[,N_placebo_start]
)
vaccine_efficacy_data[,"rr_high"] <- exp(
  log(vaccine_efficacy_data[,rr])
  +1.96*sqrt(
    (
      (vaccine_efficacy_data[,N_vaccine_start] - vaccine_efficacy_data[,cum_cases_vaccine])/vaccine_efficacy_data[,cum_cases_vaccine]
    )/vaccine_efficacy_data[,N_vaccine_start] +(
      (vaccine_efficacy_data[,N_placebo_start] - vaccine_efficacy_data[,cum_cases_placebo])/vaccine_efficacy_data[,cum_cases_placebo]
    )/vaccine_efficacy_data[,N_placebo_start]
  )
)

vaccine_efficacy_data[,"rr_low"] <- exp(
  log(vaccine_efficacy_data[,rr])
  -1.96*sqrt(
    (
      (vaccine_efficacy_data[,N_vaccine_start] - vaccine_efficacy_data[,cum_cases_vaccine])/vaccine_efficacy_data[,cum_cases_vaccine]
    )/vaccine_efficacy_data[,N_vaccine_start] +(
      (vaccine_efficacy_data[,N_placebo_start] - vaccine_efficacy_data[,cum_cases_placebo])/vaccine_efficacy_data[,cum_cases_placebo]
    )/vaccine_efficacy_data[,N_placebo_start]
  )
)

vaccine_efficacy_data[, "vaccine_efficacy"] <- "cumulative"
vaccine_efficacy_data[, "mortality_factor"] <- factor(
  vaccine_efficacy_data[, mortality],
  # levels=c("Low", "Medium", "High")
  levels=c("VE_data", "iVE_data")
)
#show extracted data by mortality stratum
#load models for pooled analysis in jcodes object (in jags language)
source("rjags_models_pooled.r")

#### r_models_pooled.r ####
#load models for pooled analaysis in rcodes object (in R language)
rcodes <- list(
  "waning_none" = function(t, mean_rr, alpha, beta){
    return(
      rep(
        exp(mean_rr),
        length(t)
      )
    )
  },
  "waning_linear" = function(t, mean_rr, alpha, beta){
    return(
      exp(mean_rr) + exp(alpha)*t
    )
  },
  "waning_power" = function(t, mean_rr, alpha, beta){
    return(
      exp(mean_rr + exp(alpha)*log(t))
    )
  },
  "waning_power_01" = function(t, mean_rr, alpha, beta){
    return(
      (
        exp(alpha)*t^exp(beta)
      )/(
        exp(mean_rr)+exp(alpha)*t^exp(beta)
      )
    )
  },
  "waning_power2" = function(t, mean_rr, alpha, beta){
    return(
      exp(mean_rr + alpha*log(t))
    )
  },
  "waning_power3" = function(t, mean_rr, alpha, beta){
    return(
      exp(mean_rr + alpha*t)
    )
  },
  "waning_sigmoid" = function(t, mean_rr, alpha, beta){
    return(
      exp( mean_rr ) * (1/( 1 + exp(alpha)*exp( -exp(beta)*t) ) )
    )
  },
  "waning_sigmoid_01" = function(t, mean_rr, alpha, beta){
    return(
      exp(mean_rr)/(
        exp(mean_rr) + exp(alpha)*exp( -exp(beta)*t )
      )
    )
  },
  "waning_gamma_01" = function(t, mean_rr, alpha, beta){
    return(
      pgamma( t, exp(alpha), exp(beta) )
    )
  }
)

#select model (waning_power)
model <- c(
  "waning_none", "waning_linear", "waning_power",
  "waning_power_01", "waning_power2", "waning_power3",
  "waning_sigmoid", "waning_sigmoid_01", "waning_gamma_01"
)[3]
#set mcmc options
mcmc.chains <- 2
mcmc.length <- 10000
thin <- 4 ## change thins to 10 if you have issues getting the P values


##### run_rjags_model.r #####
#run the models with rjags. This will create 3 data-tables: posterior_low,
# posterior_medium, and posterior_high (one for each mortality stratum)
model_code <- jcodes[[model]]
#for(m in c("Low", "Medium", "High")){
for(m in c("VE_data", "iVE_data")){
  jdat <- list(
    cases_placebo=round(vaccine_efficacy_data[mortality == m, cum_cases_placebo]),
    cases_vaccine=round(vaccine_efficacy_data[mortality == m, cum_cases_vaccine]),
    person_weeks_placebo=vaccine_efficacy_data[mortality == m, persontime_placebo],
    person_weeks_vaccine=vaccine_efficacy_data[mortality == m, persontime_vaccine],
    timept=vaccine_efficacy_data[mortality == m, follow_up_end],
    ID=vaccine_efficacy_data[mortality == m, study_id],
    uID=unique(vaccine_efficacy_data[mortality == m, study_id])
  )
  
  jmodel <- jags.model(
    textConnection(model_code),
    data=jdat,
    n.chains=mcmc.chains,
    n.adapt=1000
  )
  jposterior <- coda.samples(
    jmodel,
    c("mean_rr", "alpha", "beta"),
    n.iter=mcmc.length,
    thin=thin
  )
  #get DIC value
  jdic <- dic.samples(
    jmodel,
    n.iter=mcmc.length,
    thin=thin,
    type="pD"
  )
  
  assign(
    paste0("posterior_", tolower(m)),
    jposterior
  )
  assign(
    paste0("dic_", tolower(m)),
    jdic
  )
}
#months in which VE and iVE should be estimated
times <- seq(0.5,36,0.5)

#process modelled results
# generate VE and iVE for each saved draw from the joint posterior distribution)
# see Appendix A in the main paper for more information about the conversion process
# this will create 7 data-tables:
# 3 tables with cumulative VE for each stratum, i.e.: ve_cumulative_low
# 3 tables with instantaneous VE for each stratum, i.e.: ve_instantaneous_low
# 1 table with all VE (iVE and VE) combined: vaccine_efficacy_central

###### process_modelled_results.r ######

for(m in c("VE_data", "iVE_data")){
  #combine chains
  jposterior <- do.call(
    rbind,
    get(
      paste0("posterior_", tolower(m))
    )
  )
  veCumulative <- rcodes[[model]]
  
  ve_cumulative <- matrix(
    0,
    ncol = length(times),
    nrow = nrow(jposterior)
  )
  ve_instantaneous <- ve_cumulative
  for(r in 1:nrow(jposterior)){
    ve_cumulative[r, ] <- 1 - veCumulative(
      t = times,
      mean_rr = jposterior[r, "mean_rr"],
      alpha = if("alpha" %in% colnames(jposterior)){
        jposterior[r, "alpha"]
      } else {
        0
      },
      beta = if("beta" %in% colnames(jposterior)){
        jposterior[r, "beta"]
      } else {
        0
      }
    )
    #convert to instantaneous measure, assume no seasonality in average baseline rate
    ve_instantaneous[r, ] <- cumVEtoInsFast(
      ve_cumulative[r, ],
      0
    )
  }
  
  vaccine_efficacy_central <- rbindlist(
    list(
      data.table(
        time = times,
        vaccine_efficacy = rep(
          "cumulative",
          length(times)
        ),
        median = sapply(
          times+1,
          function(x){
            median(
              ve_cumulative[, x]
            ) 
          }
        ),
        low = sapply(
          times+1,
          function(x){
            quantile(
              ve_cumulative[, x],
              0.025
            ) 
          }
        ),
        high = sapply(
          times+1,
          function(x){
            quantile(
              ve_cumulative[, x],
              0.975
            ) 
          }
        )
      ),
      data.table(
        time = times,
        vaccine_efficacy = rep(
          "instantaneous",
          length(times)
        ),
        median = sapply(
          times+1,
          function(x){
            median(
              ve_instantaneous[, x]
            ) 
          }
        ),
        low = sapply(
          times+1,
          function(x){
            quantile(
              ve_instantaneous[, x],
              0.025
            ) 
          }
        ),
        high = sapply(
          times+1,
          function(x){
            quantile(
              ve_instantaneous[, x],
              0.975
            ) 
          }
        )
      )
    )
  )
  vaccine_efficacy_central[, "mortality"] <- m
  assign(
    paste0(
      "vaccine_efficacy_central_",
      tolower(m)
    ),
    vaccine_efficacy_central
  )
  assign(
    paste0(
      "ve_cumulative_",
      tolower(m)
    ),
    ve_cumulative
  )
  assign(
    paste0(
      "ve_instantaneous_",
      tolower(m)
    ),
    ve_instantaneous
  )
}
vaccine_efficacy_central <- rbindlist(
  list(
    vaccine_efficacy_central_ive_data,
    # vaccine_efficacy_central_medium,
    vaccine_efficacy_central_ve_data
  )
)

vaccine_efficacy_central[, "mortality_factor"] <- factor(
  vaccine_efficacy_central[, mortality],
  levels=c("VE_data", "iVE_data")
)

#calculate empirical p-values and credible intervals for set timepoints
# compares columns 'cols' in data_a with those in data_b
#default: timepoint 1 and 24 (2 weeks and 12 months follow-up)
#change data_a, data_b, and cols if needed
empiricalDifference(
  data_a = ve_instantaneous_ve_data[c(1:5000),],
  data_b = ve_instantaneous_ive_data[c(1:5000),], ### there's no medium group in Kaja's data
  cols = c(1,24,60)
)



# These are the plots - starting with Waning (Power)

vaccine_efficacy_data_VE <- vaccine_efficacy_data[ which(vaccine_efficacy_data$mortality=='VE_data'), ]
vaccine_efficacy_data_iVE <- vaccine_efficacy_data[ which(vaccine_efficacy_data$mortality=='iVE_data'), ]
vaccine_efficacy_central_VE <- vaccine_efficacy_central[ which(vaccine_efficacy_central$mortality=='VE_data' & vaccine_efficacy_central$vaccine_efficacy== "cumulative"), ]
vaccine_efficacy_central_iVE <- vaccine_efficacy_central[ which(vaccine_efficacy_central$mortality=='iVE_data' & vaccine_efficacy_central$vaccine_efficacy== "cumulative"), ]
vaccine_efficacy_central_iVE_calc <- vaccine_efficacy_central[ which(vaccine_efficacy_central$mortality=='VE_data' & vaccine_efficacy_central$vaccine_efficacy== "instantaneous"), ]

VEdat <- ggplot()+
  geom_errorbar(data=vaccine_efficacy_data_VE, aes(x=follow_up_end, ymin=1-rr_high, ymax=1-rr_low))+
  geom_point(data=vaccine_efficacy_data_VE, aes(x=follow_up_end, y=1-rr, size=(cases_vaccine + cases_placebo)))+
  coord_cartesian(ylim=c(0, 1), xlim=c(0.2, 36))+
  ggtitle("Cumulative Vaccine Efficacy") +
  labs(x="", y="Vaccine efficacy (prop)\n")+
  scale_size_continuous(range = c(0.5,2.5))+
  guides(size="none")+
  geom_hline(yintercept=0)


if(!is.null(vaccine_efficacy_central_VE)){
  VEdat <- VEdat +geom_ribbon(data=vaccine_efficacy_central_VE,aes(x=time, ymin=low, ymax=high), alpha=0.2)+
    geom_line(data=vaccine_efficacy_central_VE,aes(x=time,y=median))+
    theme_bw()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          plot.title = element_text(colour = "grey39"))
}


##### Instantaneous Vaccine Efficacy #####

iVEdat <- ggplot()+
  geom_errorbar(data=vaccine_efficacy_data_iVE, aes(x=follow_up_end, ymin=1-rr_high, ymax=1-rr_low))+
  geom_point(data=vaccine_efficacy_data_iVE, aes(x=follow_up_end, y=1-rr, size=(cases_vaccine + cases_placebo)))+
  coord_cartesian(ylim=c(0, 1), xlim=c(0.2, 36))+
  ggtitle("Instantaneous Vaccine Efficacy") +
  #facet_grid(vaccine_efficacy~mortality_factor)+
  labs(x="Months since completed vaccine schedule", y="")+
  scale_size_continuous(range = c(0.5,2.5))+
  guides(size="none")+
  geom_hline(yintercept=0)


if(!is.null(vaccine_efficacy_central_iVE)){
  iVEdat <- iVEdat +geom_ribbon(data=vaccine_efficacy_central_iVE,aes(x=time, ymin=low, ymax=high), alpha=0.2)+
    geom_line(data=vaccine_efficacy_central_iVE,aes(x=time,y=median))+
    theme_bw()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          plot.title = element_text(colour = "grey39"))
}

##### Instantaneous Vaccine Efficacy Calculated from cumulative VE #####

iVEdat_cal <- ggplot()+
  geom_ribbon(data=vaccine_efficacy_central_iVE_calc,aes(x=time, ymin=low, ymax=high), alpha=0.2)+
  geom_line(data=vaccine_efficacy_central_iVE_calc,aes(x=time,y=median))+
  ggtitle("Instantaneous Vaccine Efficacy Derived\nfrom Cumulative Vaccine Efficacy") +
  theme_bw()+
  labs(x="", y="")+
  geom_hline(yintercept=0)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        plot.title = element_text(colour = "grey39"))

grid.arrange(VEdat,iVEdat,iVEdat_cal, ncol=3)

```
