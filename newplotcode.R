### subset data to get separate plots ###
vaccine_efficacy_data_VE <- vaccine_efficacy_data[ which(vaccine_efficacy_data$mortality=='VE_data'), ]
vaccine_efficacy_data_iVE <- vaccine_efficacy_data[ which(vaccine_efficacy_data$mortality=='iVE_data'), ]
vaccine_efficacy_central_VE <- vaccine_efficacy_central[ which(vaccine_efficacy_central$mortality=='VE_data' & vaccine_efficacy_central$vaccine_efficacy== "cumulative"), ]
vaccine_efficacy_central_iVE <- vaccine_efficacy_central[ which(vaccine_efficacy_central$mortality=='iVE_data' & vaccine_efficacy_central$vaccine_efficacy== "cumulative"), ]
vaccine_efficacy_central_iVE_calc <- vaccine_efficacy_central[ which(vaccine_efficacy_central$mortality=='VE_data' & vaccine_efficacy_central$vaccine_efficacy== "instantaneous"), ]


#### Cumulative Vaccine Efficacy ####
VEdat <- ggplot()+
  geom_errorbar(data=vaccine_efficacy_data_VE, aes(x=follow_up_end, ymin=1-rr_high, ymax=1-rr_low))+
  geom_point(data=vaccine_efficacy_data_VE, aes(x=follow_up_end, y=1-rr, size=(cases_vaccine + cases_placebo)))+
  coord_cartesian(ylim=c(0, 1), xlim=c(0.2, 36))+
  #facet_grid(vaccine_efficacy~mortality_factor)+
  labs(x="Months since completed vaccine schedule", y="Vaccine efficacy (prop)\n")+
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
print(VEdat+ggtitle("Cumulative Vaccine Efficacy"))

##### Instantaneous Vaccine Efficacy #####

iVEdat <- ggplot()+
  geom_errorbar(data=vaccine_efficacy_data_iVE, aes(x=follow_up_end, ymin=1-rr_high, ymax=1-rr_low))+
  geom_point(data=vaccine_efficacy_data_iVE, aes(x=follow_up_end, y=1-rr, size=(cases_vaccine + cases_placebo)))+
  coord_cartesian(ylim=c(0, 1), xlim=c(0.2, 36))+
  #facet_grid(vaccine_efficacy~mortality_factor)+
  labs(x="Months since completed vaccine schedule", y="Vaccine efficacy (prop)\n")+
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
print(iVEdat+ggtitle("Instantaneous Vaccine Efficacy"))

##### Instantaneous Vaccine Efficacy Calculated from cumulative VE #####

iVEdat_cal <- ggplot()+
geom_ribbon(data=vaccine_efficacy_central_iVE_calc,aes(x=time, ymin=low, ymax=high), alpha=0.2)+
  geom_line(data=vaccine_efficacy_central_iVE_calc,aes(x=time,y=median))+
  theme_bw()+
  labs(x="Months since completed vaccine schedule", y="Vaccine efficacy (prop)\n")+
  geom_hline(yintercept=0)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        plot.title = element_text(colour = "grey39"))
print(iVEdat_cal+ggtitle("Instantaneous Vaccine Efficacy Derived\nfrom Cumulative Vaccine Efficacy"))


