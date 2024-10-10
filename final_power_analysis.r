options("scipen"=9, "digits"=4)
library(dplyr)
library(MASS) # includes rnegbin for generating negative binomial distributions
library(ggplot2)
library(rlang)
library(corrplot)
library(Hmisc)
library(tidyverse)
library(viridis)
library(fabricatr)
library(DeclareDesign)
library(gmodels)
library(glmmTMB)
library(future)

plan(multicore) 

df <- read.csv("~/DMCA-exploratory-code-and-data/data/sim_df.csv")

df <- df[df$day_num < 23,]
df <- df[df$day_num > -23,]
df$day_num.factor <- factor(df$day_num)

handler <- function(data) {

    data <- data %>% arrange(twitter_user_id_unique, day_num.factor)

    data <- data %>%
        mutate(NTD = if_else(user_suspended == "True", 0, NTD))
        
    return(data)
    
}



## Integer acquired from London Integers
## https://londonintegers.com/int/40107429
set.seed(40107429)


config.df <- data.frame(
    pa.label = "DMCA design",
    n.min = 10000,
    n.max = 18000,
    n.days = 44,
    
    effect.magnitude = -0.03, #assumption: post-DMCA tweets decrease by 3%
    effect.slope   = -0.003 #assumption: post-DMCA tweets decrease by 0.3% per day

)

diagnose.experiment <- function(n.size, config.df, sims.count=200, bootstrap.sims.count=100){
    
    design <- declare_model(
        data = resample_data(df, N = c(twitter_user_id = n.size), unique_labels = TRUE)
    ) +
    
    declare_potential_outcomes(
        NTD_X_0 = rnegbin(n.size*config.df$n.days, mu=y.mu, theta=y.theta),
        NTD_X_1 = rnegbin(n.size*config.df$n.days, mu=y.mu + y.mu*config.df$effect.magnitude + y.mu * config.df$effect.slope*day_num , theta=y.theta)

    ) +

    declare_measurement(
            NTD = reveal_outcomes(NTD ~ X)
        
    )  +

    declare_step(

        handler = handler
        
    ) +
    
    declare_inquiry(
            NTD_PATE = config.df$effect.magnitude,
            NTD_SLOPE = config.df$effect.slope
    ) +
    declare_estimator(
        NTD ~ X + day_num + X:day_num + (1|twitter_user_id_unique),
        .method = glmmTMB,
        family = "nbinom2",
        ziformula = ~ user_suspended,
        .summary = broom.mixed::tidy,
        term = c("X", "X:day_num"),
        inquiry = c("NTD_PATE", "NTD_SLOPE")
    )


    tryCatch({
    # Run the diagnose_design function
    diagnosis <- diagnose_design(design, sims = sims.count, 
                                 bootstrap_sims = bootstrap.sims.count)
    diagnosis
    
    }, error = function(e) {
    # Save the problematic dataset and error message
    error.df <- data

    write.csv(error.df, "/home/lawright/data/error_data.csv", row.names = FALSE)
    
    
    return(NULL)  # Return NULL or handle as necessary
    })
    
    
}

min.diagnosis.power <- function(diagnosis){
    min(diagnosis$diagnosands_df['power'])
}


# Iterate linearly for a certain level of statistical power
# within the constraints of a configuration file
# at a certain sample size increment. Useful for
# illustrating ideas, or for comparing estimators with
# very different statistical power, where the binary search
# will optimize for the worst estimator but not show useful
# indormation about more efficient estimators
#
#` @param config.df The configuration file in question
#` @diagnosis.method The method that conducts a single DeclareDesign diagnosis and returns the diagnosis
#` @iteration.interval when iterating, use this interval between sample sizes

iterate.for.power <- function(config.df, diagnosis.method = diagnose.experiment, 
                             iteration.interval){  
    max.sample.size = config.df$n.max
    min.sample.size = config.df$n.min
    current.sample.size = min.sample.size
    
    iteration.count = ceiling((max.sample.size - min.sample.size) / iteration.interval)

    ## Initialize first iteration
    print(paste("min:", min.sample.size, "max:", max.sample.size, "current:", current.sample.size))
    flush.console()

    ptm = proc.time()
    ddf <- diagnosis.method(current.sample.size, config.df)
    ddf$diagnosands$n <- current.sample.size
    diagnoses.df = ddf$diagnosands
    current.power <- min.diagnosis.power(ddf)
    time.elapsed <- proc.time() -  ptm
    print(paste("     seconds:", as.integer(time.elapsed['elapsed'])))
    
    for(i in seq(1, iteration.count)){
        current.sample.size = current.sample.size + iteration.interval
        print(paste("min:", min.sample.size, "max:", max.sample.size, "current:", current.sample.size))
        flush.console()
    
        ptm = proc.time()
        ## conduct simulations

        tryCatch({
            # Run the diagnose_design function
            ddf <- diagnosis.method(current.sample.size, config.df)
            
          }, error = function(e) {
            # Save the problematic dataset and error message
            filename <- paste0("error_data_iteration_", iteration, ".RDS")
            saveRDS(design, file = filename)  # Save the design object
            
            cat("Error in iteration", iteration, ":", conditionMessage(e), "\n")
            return(NULL)  # Return NULL or handle as necessary
          })        

        
        ddf <- diagnosis.method(current.sample.size, config.df)
        ddf$diagnosands$n <- current.sample.size       
        ## append simulation results to dataframe
        diagnoses.df <- rbind(diagnoses.df, ddf$diagnosands)
        time.elapsed <- proc.time() -  ptm
        print(paste("     seconds:", as.integer(time.elapsed['elapsed'])))
    }
#    diagnoses.df[diagnoses.df$inquiry == "NTD_PATE" & diagnoses.df$term == "X" | diagnoses.df$inquiry == "NTD_SLOPE" & diagnoses.df$term == "X:day_num", ]
    diagnoses.df
}



interval = 1000
power.iterate.df <- iterate.for.power(config.df, iteration.interval = interval)

plot <- ggplot(power.iterate.df, aes(n, power, color=term)) +
    ## CHART SUBSTANCE
    geom_line() +
    geom_point() +
    ## LABELS AND COSMETICS
    geom_hline(yintercept=0.95, size=0.25) +
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1), labels=scales::percent) +
    scale_x_continuous(breaks = seq(config.df$n.min,config.df$n.max+10,interval)) +
    scale_color_viridis(discrete=TRUE) +
    xlab("sample size") +
    ggtitle("Statistical Power Associated with Estimators")

ggsave("power_plot.png")

write.csv(power.iterate.df, "~/DMCA-exploratory-code-and-data/data/loop_power_results_final.csv", row.names = FALSE)
