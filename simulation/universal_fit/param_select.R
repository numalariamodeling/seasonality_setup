library(tidyverse)
library(stringr)
library(glue)
library(lubridate)
library(patchwork)

iodir <- "./"

load_case_data <- function(projectdir) {
  # Function to load routine case data
  case_fname <- file.path(projectdir, 'burkina_cases/routine_seb_agg_confpres.csv')
  case_df <- data.table::fread(case_fname) %>% 
    rename(date = Date, repincd = case) %>%
    mutate(date = as.Date(date), year = year(date), 
           month = month(date)) %>%
    arrange(DS_Name, date)
  return(case_df)
}

# Global variables and prior distributions
batch <- "batch_1"
start_year = 1960
rnd <- 2
DS <- "BK"
params_name <- c(
  "Temperature_Shift",
  "CONST_Multiplier",
  "TEMPR_Multiplier",
  "WATEV_Multiplier"
)
prior_mean <- c(0, 0, 0, 0)
prior_sd <- c(1.5, 1.5, 1.5, 1.5)
lbound <- c(-5, -3, -3, -3)
ubound <- -lbound
space <- c("n", "l", "l", "l") # log normal for the multipliers
npar <- 500
nsel <- 25


# Generate particles for new round
# Round 0 - No particle selection, generate particles using priors
# Other rounds - Conduct particle selection, then add samples by perturbing
# selected particles
if (rnd == 0) {
  set.seed(4326)
  mat <- matrix(NA, nrow = npar, ncol = length(params_name))
  for (i in 1:length(params_name)) {
    tmp <- rnorm(npar, prior_mean[i], prior_sd[i])
    tmp <- ifelse(tmp < lbound[i], lbound[i] + (lbound[i] - tmp), tmp)
    tmp <- ifelse(tmp > ubound[i], ubound[i] - (tmp - ubound[i]), tmp)
    if (space[i] == "l") tmp <- 10^tmp
    mat[,i] <- tmp
  }
  
  colnames(mat) <- params_name
  df <- as.data.frame(mat)
  df$samp_id <- 1:nrow(df)
  perf <- data.frame(round = c(""), parameter=c(""), mean= c(""), marginal_sd=c(""))
  
} else {
  set.seed(4326+rnd)
  simoutdir <- file.path(iodir, glue("simulation_output/{DS}/{batch}/rnd{rnd-1}"))
  samp_df <- file.path(iodir, glue("simulation_output/{DS}/{batch}/rnd{rnd-1}.csv")) |>
    data.table::fread()
  # Normalize over 5 routine case data by annual sum as reference data
  rcases <- load_case_data("/projects/b1139/malaria-bf-hbhi/IO") |>
    filter(DS_Name == "Sapone") %>%
    filter(age == 'ov5')
  rcases1 <- rcases |>
    filter(DS_Name == "Sapone") |>
    group_by(year) |>
    mutate(norm_repincd = repincd/sum(repincd)) |>
    group_by(month) |>
    summarise(norm_repincd = mean(norm_repincd))
  
  # Calculate monthly treated cases in simulation
  simcases <- data.table::fread(file.path(simoutdir, "PfPR_ClinicalIncidence_monthly.csv"))
  simcases <- simcases |>
    filter(agebin == 5) |>
    group_by(Sample_ID, year=Year, month) |>
    summarise(pop = mean(Pop),
              case = mean(Cases)) |>
    mutate(date = paste(year, month) |> ym())
  
  simtreated <- data.table::fread(file.path(simoutdir, "events.csv")) |>
    filter((Age %/% 365 >= 5)) |>
    mutate(year = (Time - 1) %/% 365,
           month = (Time - year * 365 - 1) %/% 30 + 1) |>
    filter(month <= 12)
  
  simtreated <- simtreated |>
    group_by(Sample_ID, year, month, Run_Number) |>
    summarise(simincd = n()) |>
    summarise(simincd = mean(simincd), .groups = "drop") |>
    mutate(year = year + start_year)
  
  simtreated1 <- simtreated %>%
    left_join(simcases, by = c("Sample_ID", "year", "month")) %>%
    mutate(simincd = simincd / pop)
  
  # Normalize monthly simulated treated cases
  simcases1 <- simtreated1 |>
    group_by(Sample_ID, year) |>
    mutate(norm_simincd = simincd / sum(simincd)) |>
    group_by(Sample_ID, month) |>
    summarise(norm_simincd = mean(norm_simincd,na.rm=T))
  
  ic <- data.table::fread(file.path(simoutdir, "InsetChart.csv"))
  
  refpcr <- data.table::fread("/projects/b1139/FE_tmh6260/project_template/simulation_inputs/sapone/reference_data/pcr_prevalence_allAge.csv")
  refpcr <- refpcr %>% mutate(variable="PCR Parasite Prevalence") %>% rename(time=sim_day,value=pcr_prevalence)
  
  refpcr <- refpcr %>% rowwise() %>% mutate(month=unlist(strsplit(date,'/'))[1],
                                            year=unlist(strsplit(date,'/'))[3])
  
  score4 <- ic %>% 
    mutate(date = as.Date(time, origin="1960-01-01")) %>% 
    mutate(month=format.Date(date,format="%m-%y")) %>%
    group_by(Sample_ID, Run_Number,month) %>%
    summarize(EIR = sum(`Daily EIR`)) %>%
    group_by(Sample_ID) %>% summarize(max_EIR = max(EIR)) %>%
    group_by(Sample_ID) %>% summarize(score_EIR=5*(max_EIR>=100))
  
  score3 <- ic %>% mutate(year = (time - 1) %/% 365,
                          month = (time - year * 365 - 1) %/% 30 + 1) %>%
    mutate(year=year+start_year) %>%
    filter(month <= 12) %>%
    filter(paste(month,year,sep='-') %in% paste(refpcr$month,refpcr$year,sep='-')) %>%
    group_by(Sample_ID,year,month) %>%
    summarize(prevalence=mean(`PCR Parasite Prevalence`)) %>%
    mutate(month=paste(month),year=paste(year)) %>%
    left_join(refpcr)  %>%
    rowwise() %>%
    mutate(prev_score=sqrt(abs(prevalence-value)^2)) %>%
    group_by(Sample_ID) %>% summarize(prev_score=mean(prev_score))
  
  
  # Score simulated vs reference
  # Score 1 is a "shape score": seasonality shape matching (MSE)
  # Score 2 is an "intensity score": annual incidence should be close to
  # 2.5 per person (MAE)
  # Combined score is 200 shape score to 1 intensity score (to make sure
  # shape score is about the same scale as intensity score)
  score1 <- simcases1 |> left_join(rcases1, by = "month") |>
    summarise(shape_score = mean((norm_simincd - norm_repincd)^2) |> sqrt()) |>
    arrange(shape_score)
  score2 <- simcases |>
    group_by(Sample_ID, year) |>
    summarise(case = mean(case)) |>
    summarise(case = mean(case)) |>
    mutate(intensity_score = abs(case - 2.5))
  
  scores <- left_join(score1, score2, by="Sample_ID") |>
    left_join(score3) |>
    left_join(score4) |>
    mutate(final_score = shape_score * 200 + intensity_score + prev_score + score_EIR) |>
    arrange(final_score) |>
    mutate(rank = 1:n())
  
  score_sel <- scores$Sample_ID[scores$rank %in% 1:nsel]
  
  lookup <- scores %>% 
    filter(Sample_ID %in% score_sel) %>% 
    left_join(samp_df %>% rename(Sample_ID = samp_id)) %>% 
    mutate(param_set = paste(rank,"- Sample #:",Sample_ID,"\nTemp Shift:",round(Temperature_Shift,2),
                             '\nC:     TR:    WV:\n',round(CONST_Multiplier,2),round(WATEV_Multiplier,2),round(TEMPR_Multiplier,2), sep=' '))
  
  # Plot top 25 sims against reference
  
  age_pop <- data.table::fread(file.path(simoutdir,"age_population.csv"))
  age_pop <- age_pop %>% filter(AgeYears %in% c(5,50)) %>%
    mutate(age_group=AgeYears) %>% 
    mutate(age_group=ifelse(age_group==5,"Under 5",age_group)) %>%
    mutate(age_group=ifelse(age_group==50,"Over 5",age_group)) %>%
    group_by(age_group,Time,Sample_ID) %>% 
    summarize(pop = mean(NumIndividuals)) 
  
  iv <- data.table::fread(file.path(simoutdir, "event_counts.csv")) 
  
  iv <- iv %>% 
    filter(Sample_ID %in% score_sel[1:nsel]) %>%
    mutate(age_group=ifelse(Age_Year<5,"Under 5","Over 5"))%>%
    group_by(Sample_ID,Run_Number,Time,Event_Name,age_group) %>% summarize(count=sum(Individual_ID)) %>%
    group_by(Sample_ID,Time,Event_Name,age_group) %>%
    summarize(count=mean(count)) %>%
    left_join(samp_df %>% rename(Sample_ID = samp_id)) %>%
    left_join(scores) %>%
    left_join(age_pop)
  
  ev2 <-  iv %>% 
    pivot_wider(names_from = Event_Name, values_from = count) %>%
    ggplot(aes(x=as.Date(Time,origin=paste(start_year,"-01-01",sep="")))) +
    facet_wrap(~age_group, scales="free", ncol=1) +
    geom_point(aes(y=Bednet_Got_New_One/pop, color="Got New Bednet"), alpha=0.3) +
    geom_path(aes(y=Bednet_Using/pop, group=rank, color="Using Bednet"), alpha=0.3) +
    geom_path(aes(y=Bednet_Discarded/pop, group=rank, color="Discarded Bednet"), alpha=0.3) +
    geom_path(aes(y=Received_Treatment/pop, group=rank, color="Received Treatment"), alpha=0.3) +
    geom_point(aes(y=Received_SMC/pop,group=rank,color="Received_SMC"),alpha=0.3)+
    theme_minimal(base_size=10) +
    theme(legend.position="right")+
    ggtitle("Simulated Interventions") +
    labs(color=NULL) +
    guides(color=guide_legend(ncol = 1, override.aes = list(alpha=1, linewidth=5), title.position = "top", title.hjust = 0.5)) + 
    ylab("Fraction of Age Group Population") +xlab("") +
    scale_x_date(date_breaks="1 year", date_labels = "%b\n%Y", date_minor_breaks = "3 month", expand=c(0,0)) +
    scale_color_brewer(palette="Set1")
  
  channels = c("Air Temperature","Rainfall","Adult Vectors","Daily Bites per Human", "Infectious Vectors", "Daily EIR","PCR Parasite Prevalence")
  
  
  p2<- ggplot() +
    geom_point(aes(x=month, y=norm_repincd, shape="Data"), data=rcases1, color="black") +
    geom_line(aes(x=month, y=norm_simincd, group=Sample_ID, color=rank),
              data=simcases1 |> left_join(lookup) |> filter(!is.na(rank)),
              alpha = 0.2) +
    theme_minimal(base_size=10) +
    theme(legend.position="right", legend.box="vertical", legend.direction = "vertical") +
    ggtitle("Normalized Clinical Incidence - Ages 5-50")  +
    scale_color_distiller(palette="YlGnBu") +
    labs(color="Parameter Rank", shape=NULL) +
    guides(color=guide_colorbar(title.position = "top", title.hjust = 0.5, direction = "vertical")) +
    scale_x_continuous(breaks=seq(1,12,1), labels=month.abb, minor_breaks = NULL) +
    xlab("")
  
  ic2 <- ic %>% 
    filter(Sample_ID %in% score_sel[1:nsel]) %>%
    left_join(scores %>% select(Sample_ID,rank)) %>%
    left_join(samp_df |> rename(Sample_ID=samp_id)) %>%
    select(-c(Run_Number,day,year))  %>%
    group_by(Sample_ID,time,rank, Temperature_Shift, WATEV_Multiplier, CONST_Multiplier, TEMPR_Multiplier) %>% 
    summarize_all(mean)  %>%
    gather("variable","value",-c(1:7)) %>%
    ungroup() %>%
    ggplot(aes(x=as.Date(time,origin=paste(start_year,"-01-01",sep='')),y=value)) +
    facet_wrap(~factor(variable,levels=channels), ncol=1, scales="free_y") +
    geom_path(aes(group=Sample_ID, color=rank),alpha=0.2) +
    theme_minimal(base_size=10) +
    geom_point(data=refpcr, aes(x=as.Date(date, format="%m/%d/%Y"), y=value,shape="Data"), color="black") +
    theme(legend.position="none", legend.direction = 'horizontal') +
    xlab("") + ylab("") +
    labs(color="Parameter Set Ranking", shape=NULL) +
    scale_x_date(date_breaks="1 year", date_labels="%b\n%Y", date_minor_breaks = "1 month", limits = c(as.Date("2015-01-01"),NA),expand=c(0,0)) + 
    guides(color=guide_colorbar(title.position = "top", title.hjust = 0.5, direction = "horizontal")) +
    ggtitle("Simulation Outputs - All Ages") +
    scale_color_distiller(palette = "YlGnBu")
  
  my_layout <- "CCCBBB
                CCCBBB
                CCCAAA
                CCCAAA
                CCCAAA
                CCCAAA
                CCCAAA"
  
  ev2 + p2 + ic2 + plot_layout(design = my_layout)
  
  ggsave(glue("simulation_output/{DS}/{batch}/rnd{rnd-1}/IC_stacked_{nsel}fits_rnd{rnd-1}.png"), width=18,height=9)
  ggsave(glue("simulation_output/{DS}/{batch}/rnd{rnd-1}/IC_stacked_{nsel}fits_rnd{rnd-1}.pdf"), width=18,height=9)
  
  ev1 <-  iv %>% filter(Sample_ID == score_sel[1]) %>%
    pivot_wider(names_from = Event_Name, values_from = count) %>%
    ggplot(aes(x=as.Date(Time,origin=paste(start_year,"-01-01",sep="")))) +
    facet_wrap(~age_group, scales="free", ncol=1) +
    geom_point(aes(y=Bednet_Got_New_One/pop, color="Got New Bednet"), alpha=0.7) +
    geom_path(aes(y=Bednet_Using/pop, group=rank, color="Using Bednet"), alpha=0.7) +
    geom_path(aes(y=Bednet_Discarded/pop, group=rank, color="Discarded Bednet"), alpha=0.7) +
    geom_path(aes(y=Received_Treatment/pop, group=rank, color="Received Treatment"), alpha=0.7) +
    geom_point(aes(y=Received_SMC/pop,group=rank,color="Received_SMC"),alpha=0.7)+
    theme_minimal(base_size=10) +
    theme(legend.position="right")+
    ggtitle("Simulated Interventions") +
    labs(color=NULL) +
    guides(color=guide_legend(ncol = 1, override.aes = list(alpha=1, linewidth=5), title.position = "top", title.hjust = 0.5)) + 
    ylab("Fraction of Age Group Population") +xlab("") +
    scale_x_date(date_breaks="1 year", date_labels = "%b\n%Y", date_minor_breaks = "3 month", expand=c(0,0)) +
    scale_color_brewer(palette="Set1")
  
  
  
  p1<- ggplot() +
    geom_point(aes(x=month, y=norm_repincd, shape="Data"), data=rcases1, color="black") +
    geom_line(aes(x=month, y=norm_simincd, group=Sample_ID, color=rank),
              data=simcases1 |> left_join(lookup) |> filter(!is.na(rank)) |> filter(rank==1),
              alpha = 1) +
    theme_minimal(base_size=10) +
    theme(legend.position="right", legend.box="vertical", legend.direction = "vertical") +
    ggtitle("Normalized Clinical Incidence - Ages 5-50")  +
    scale_color_distiller(palette="YlGnBu") +
    labs(color="Parameter Rank", shape=NULL) +
    guides(color=guide_colorbar(title.position = "top", title.hjust = 0.5)) +
    scale_x_continuous(breaks=seq(1,12,1), labels=month.abb, minor_breaks = NULL) +
    xlab("")
  
  ic1 <- ic %>% 
    filter(Sample_ID %in% score_sel[1]) %>%
    left_join(scores %>% select(Sample_ID,rank)) %>%
    left_join(samp_df |> rename(Sample_ID=samp_id)) %>%
    select(-c(Run_Number,day,year))  %>%
    group_by(Sample_ID,time,rank, Temperature_Shift, WATEV_Multiplier, CONST_Multiplier, TEMPR_Multiplier) %>% 
    summarize_all(mean)  %>%
    gather("variable","value",-c(1:7)) %>%
    ungroup() %>%
    ggplot(aes(x=as.Date(time,origin=paste(start_year,"-01-01",sep='')),y=value)) +
    facet_wrap(~factor(variable,levels=channels), ncol=1, scales="free_y") +
    geom_path(aes(group=Sample_ID, color=rank),alpha=1) +
    theme_minimal(base_size=10) +
    geom_point(data=refpcr, aes(x=as.Date(date, format="%m/%d/%Y"), y=value,shape="Data"), color="black") +
    theme(legend.position="none", legend.box = "horizontal") +
    xlab("") + ylab("") +
    labs(color="Parameter Set Ranking", shape=NULL) +
    scale_x_date(date_breaks="1 year", date_labels="%b\n%Y", date_minor_breaks = "1 month", limits = c(as.Date("2015-01-01"),NA),expand=c(0,0)) + 
    guides(color=guide_colorbar(title.position = "top", title.hjust = 0.5)) +
    ggtitle("Simulation Outputs - All Ages") +
    scale_color_distiller(palette = "YlGnBu")
  
  
  
  my_layout <- "CCCBBB
                CCCBBB
                CCCAAA
                CCCAAA
                CCCAAA
                CCCAAA
                CCCAAA"
  
  ev1 + p1 + ic1 + plot_layout(design = my_layout)
  
  ggsave(glue("simulation_output/{DS}/{batch}/rnd{rnd-1}/IC_stacked_bestfit_rnd{rnd-1}.png"), width=18,height=9)
  ggsave(glue("simulation_output/{DS}/{batch}/rnd{rnd-1}/IC_stacked_bestfit_rnd{rnd-1}.pdf"), width=18,height=9)
  
  
  # Retrieve the params for the top 25 sample ID
  samp_sel <- scores |>
    filter(rank %in% 1:nsel) |>
    select(samp_id = Sample_ID, rank, final_score, shape_score, intensity_score) |>
    left_join(samp_df, by = "samp_id")
  logscale_par <- params_name[space == "l"]
  samp_sel[,logscale_par] <- log10(samp_sel[,logscale_par])
  
  # Calculate multivariate mean and covariance matrix
  cat("Mean\n")
  samp_sel[,params_name] %>%
    colMeans %>%
    print
  
  cat("Marginal SD\n") # This should shrink round over round
  samp_sel[,params_name] %>%
    as.matrix() %>%
    cov %>%
    diag %>%
    sqrt %>%
    print
  
  cat("Correlation matrix\n")
  samp_sel[,params_name] %>%
    as.matrix() %>%
    cor %>%
    print
  
  V <- samp_sel[,params_name] %>%
    as.matrix() %>%
    cov
  
  perf <- read.csv(glue("simulation_output/{DS}/{batch}/performance.csv"))
  
  a<-samp_sel[,params_name] %>% data.frame() %>% 
    summarize_all(mean) %>% gather("parameter","value") %>% 
    mutate(round=rnd-1) %>% mutate(variable="mean")
  
  samp_sel[,params_name] %>% as.data.frame() %>%
    cov() %>% diag() %>% sqrt() -> b
  
  b <- data.frame(value=b, parameter=params_name, row.names = NULL) %>% mutate(round=rnd-1, variable="marginal_sd")
  
  pp <- rbind.data.frame(a,b)
  
  pp %>% pivot_wider(values_from = value,names_from = variable) -> pp
  
  perf <- rbind(pp,perf %>% select(-c(X,X.1)))
  
  pp1 <- perf %>% ggplot(aes(x=round)) +
    geom_point(aes(y=ifelse(parameter=="Temperature_Shift",mean,mean))) +
    geom_line(aes(y=ifelse(parameter=="Temperature_Shift",mean,mean))) +
    geom_segment(aes(xend=round, 
                     y=ifelse(parameter=="Temperature_Shift",mean-marginal_sd,(mean-marginal_sd)),
                     yend=ifelse(parameter=="Temperature_Shift",mean+marginal_sd,(mean+marginal_sd)))) +
    facet_wrap(~parameter, ncol=1) +
    theme_minimal() +
    theme(strip.text=element_text(size=12))+
    ylab("Mean ± Marginal SD") + xlab("Fitting Round")
  
  pp2 <- perf %>% ggplot(aes(x=round)) +
    geom_point(aes(y=ifelse(parameter=="Temperature_Shift",mean,10^mean))) +
    geom_line(aes(y=ifelse(parameter=="Temperature_Shift",mean,10^mean))) +
    geom_segment(aes(xend=round, 
                     y=ifelse(parameter=="Temperature_Shift",mean-marginal_sd,10^(mean-marginal_sd)),
                     yend=ifelse(parameter=="Temperature_Shift",mean+marginal_sd,10^(mean+marginal_sd)))) +
    facet_wrap(~parameter, scales="free_y",ncol=1) +
    theme_minimal() +
    theme(strip.text=element_text(size=12)) +
    ylab("Mean ± Marginal SD") + xlab("Fitting Round")
  
  pp1 | pp2
  
  ggsave(glue("simulation_output/{DS}/{batch}/performance_plot.png"), width=8,height=4)
  ggsave(glue("simulation_output/{DS}/{batch}/performance_plot.pdf"), width=8,height=4)
  
  
  
  eir1 <- ic %>% filter(Sample_ID %in% score_sel[1:nsel]) %>%
    mutate(date=as.Date(time,origin="1960-01-01")) %>%
    mutate(month=month(date), Year=format.Date(date,format="%Y")) %>%
    filter(Year %in% seq(2010,2019)) %>%
    group_by(Sample_ID) %>% mutate(max_bites=max(`Daily Bites per Human`)) %>%
    #ungroup() %>% filter(max_bites<=200) %>%
    group_by(Sample_ID, Run_Number, month, Year) %>%
    summarize(monthly_EIR = sum(`Daily EIR`)) %>%
    group_by(Sample_ID,month,Year) %>%
    summarize(monthly_EIR = mean(monthly_EIR)) %>%
    #group_by(Sample_ID) %>% mutate(max_EIR=max(monthly_EIR)) %>%
    #ungroup() %>% filter(max_EIR<=100) %>%
    left_join(lookup) %>%
    ggplot(aes(x=as.Date(paste(month,"01",Year,sep="/"),format="%m/%d/%Y"),
               y=monthly_EIR, group=Sample_ID)) +
    geom_line(aes(color=rank), alpha=0.7) +
    scale_color_distiller(palette="Spectral") +
    theme_minimal() +
    theme(strip.text = element_blank(), legend.position="right") +
    facet_wrap(~trunc(as.integer(paste(Year))/5), ncol=1,scales="free") +
    scale_x_date(date_minor_breaks = '1 month', expand = c(0,0)) +
    xlab("") + ylab("Monthly EIR (mean of sum of Daily EIR)") +
    labs(color="Rank") +
    ggtitle(paste("Simulated monthly EIR vs. Time for top",nsel,"parameter sets",sep=" ")) 
  
  
  vr <- data.table::fread(file.path(simoutdir,"VectorStats.csv"))
  v1 <- vr %>% 
    filter(Sample_ID %in% score_sel[1:nsel]) %>%
    filter(Time >= max(vr$Time)-365) %>%
    group_by(Time, Species, Sample_ID) %>%
    summarize(AvailableHabitat=mean(AvailableHabitat),
              VectorPopulation=mean(VectorPopulation)) %>% 
    group_by(Time,Sample_ID) %>%
    mutate(TotalVectorPopulation=sum(VectorPopulation)) %>%
    left_join(scores) %>%
    ggplot(aes(x=as.Date(Time, origin="1960-01-01"), y=VectorPopulation/TotalVectorPopulation, color=Species)) +
    geom_path(aes(group=interaction(Species,Sample_ID), color=rank), alpha=0.5, size=1.5) +
    facet_wrap(~Species, nrow=1) +
    theme_minimal() +
    scale_x_date(date_breaks = "2 months", date_minor_breaks="1 month", date_labels="%b\n%Y", expand = c(0,0)) +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    ylab("Fraction of Vector Population") +
    xlab("") + 
    scale_color_distiller(palette="Spectral") +
    labs(color="Parameter Ranking") +
    theme(legend.position="none", strip.text=element_text(size=12))
  
  my_layout <- "BBBB
                  AAAA
                  AAAA"
  
  eir1 + v1 + plot_layout(design = my_layout)    
  
  ggsave(glue("simulation_output/{DS}/{batch}/rnd{rnd-1}/vector_EIR_plot.png"), width=8,height=10)
  ggsave(glue("simulation_output/{DS}/{batch}/rnd{rnd-1}/vector_EIR_plot.pdf"), width=8,height=10)
  
  
  
  
  #scale_x_continuous(breaks=seq(1,12,1),labels = month.abb, 
  # minor_breaks = NULL)
  
  # For each of the top 25 particles, retain them for next round
  # and add more samples by perturbing them with 2x the covariance 
  # matrix above
  nsamp_per_particle <- (npar - nsel)/nsel
  perturb_samp <- samp_sel[,params_name] %>% as.matrix()
  for (i in 1:nrow(samp_sel)) {
    m <- samp_sel[i,params_name] %>% as.matrix()
    samp <- mvtnorm::rmvnorm(nsamp_per_particle, m, 2 * V)
    colnames(samp) <- params_name
    
    perturb_samp <- rbind(perturb_samp, samp)
  }
  
  # If exceeding the bounds, adjust them
  for (j in 1:length(params_name)) {
    tmp <- perturb_samp[,j]
    tmp <- ifelse(tmp < lbound[j], lbound[j] - (lbound[j] - tmp), tmp)
    tmp <- ifelse(tmp > ubound[j], ubound[j] - (tmp - ubound[j]), tmp)
    if (space[j] == "l") tmp <- 10^tmp
    perturb_samp[,j] <- tmp
  }
  
  df <- as.data.frame(perturb_samp)
  df$samp_id <- 1:nrow(df)
  print(head(df))
}

# Output
outdir <- file.path(iodir, glue("simulation_output/{DS}/{batch}"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = T)
write.csv(x = perf, glue("simulation_output/{DS}/{batch}/performance.csv"))
fname <- glue("rnd{rnd}.csv")
data.table::fwrite(df, file.path(outdir, fname))
print(paste0(fname, ' saved under ', outdir))

