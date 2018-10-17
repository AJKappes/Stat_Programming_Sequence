library(ggplot2)
library(reshape2)
library(Synth)
library(dplyr, warn.conflicts = FALSE)

df <- read.csv('synth_smoking.csv')

# b #

dataprep_b_input <- dataprep(
  foo = df,
  predictors = c('lnincome', 'beer', 'age15to24', 'retprice'),
  dependent = 'cigsale',
  unit.variable = 'state',
  time.variable = 'year',
  special.predictors = list(list('cigsale', 1975, 'mean'),
                         list('cigsale', 1980, 'mean'),
                         list('cigsale', 1988, 'mean')),
  treatment.identifier = 3,
  controls.identifier = c(1, 2, seq(4, 39, 1)),
  time.predictors.prior = c(seq(1985, 1988, 1)),
  time.optimize.ssr = c(seq(1970, 1988, 1)),
  time.plot = c(sort(unique(df$year)))
)

synth_b_output <- synth(dataprep_b_input)

gaps_b <- dataprep_b_input$Y1plot-(dataprep_b_input$Y0plot%*%synth_b_output$solution.w)
synth_b_te <- mean(gaps_b[20:31])

# path.plot(synth_b_output, dataprep_b_input,
#           Ylab = c("Cigarette Sales"),
#           Xlab = c("Year"),
#           tr.intake = 1989,
#           Legend = c('treated', 'synthetic'),
#           Legend.position = c('top'),
#           Main = 'Synthetic control, treatment state 3 year > 1988, matched 1970, 1980, 1988'
#          )

# c #

dataprep_c_input <- dataprep(
  foo = df,
  predictors = c('lnincome', 'beer', 'age15to24', 'retprice'),
  dependent = 'cigsale',
  unit.variable = 'state',
  time.variable = 'year',
  special.predictors = list(list('cigsale', 1975:1988, 'mean')),
                            #list('cigsale', 1980, 'mean'),
                            #list('cigsale', 1988, 'mean')),
  treatment.identifier = 3,
  controls.identifier = c(1, 2, seq(4, 39, 1)),
  time.predictors.prior = c(seq(1970, 1988, 1)),
  time.optimize.ssr = c(seq(1970, 1988, 1)),
  time.plot = c(sort(unique(df$year)))
)

synth_c_output <- synth(dataprep_c_input)

gaps_c <- dataprep_c_input$Y1plot-(dataprep_c_input$Y0plot%*%synth_c_output$solution.w)
synth_c_te <- mean(gaps_c[20:31])

# path.plot(synth_c_output, dataprep_c_input,
#           Ylab = c("Cigarette Sales"),
#           Xlab = c("Year"),
#           tr.intake = 1989,
#           Legend = c('treated', 'synthetic'),
#           Legend.position = c('top'),
#           Main = 'Synthetic control, treatement state 3 year > 1988, matched 1970-88'
#           )

loss <- cbind(synth_b_output$loss.w, synth_c_output$loss.w)
colnames(loss) <- c('1975, 1980, 1998', '1970-88')
rownames(loss) <- 'w_weights'

# d #
placebo_list <- sort(unique(df$state))

placebo_fun <- function(s){
  
    placebo_input <- dataprep(
    foo = df,
    predictors = c('lnincome', 'beer', 'age15to24', 'retprice'),
    dependent = 'cigsale',
    unit.variable = 'state',
    time.variable = 'year',
    special.predictors = list(list('cigsale', 1975, 'mean'),
                              list('cigsale', 1980, 'mean'),
                              list('cigsale', 1988, 'mean')),
    treatment.identifier = s,
    controls.identifier = c(placebo_list[placebo_list != s]),
    time.predictors.prior = c(seq(1985, 1988, 1)),
    time.optimize.ssr = c(seq(1970, 1988, 1)),
    time.plot = c(sort(unique(df$year)))
  )
    
    placebo_output <- synth(placebo_input)
    placebo_gaps <- placebo_input$Y1plot-(placebo_input$Y0plot%*%placebo_output$solution.w)
    
    return(placebo_gaps)
  
}

placebo_df <- placebo_fun(1)

i <- 2
while (i < length(placebo_list) + 1){
  
  placebo_df <- cbind(placebo_df, placebo_fun(i))
  i = i + 1

}

placebo_df <- data.frame(cbind(placebo_df, sort(unique(df$year))))

placebo_names <- c('s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10',
                   's11', 's12', 's13', 's14', 's15', 's16', 's17', 's18', 's19',
                   's20', 's21', 's22', 's23', 's24', 's25', 's26', 's27', 's28',
                   's29', 's30', 's31', 's32', 's33', 's34', 's35', 's36', 's37', 
                   's38', 's39', 'year')

colnames(placebo_df) <- placebo_names
rownames(placebo_df) <- NULL
state <- c(colnames(placebo_df)[1:39])

placebo_plot_df <- melt(placebo_df, id = c('year'))

placebo_filter <- placebo_plot_df %>% group_by(variable) %>%
  filter(variable == 's3') %>% ungroup()

placebo_map <- ggplot() +
  geom_line(aes(x=year, y=value, group = variable), data = placebo_plot_df,
            colour = alpha('grey', 0.6)) +
  geom_line(aes(x=year, y=value, colour = variable), data = placebo_filter) +
  scale_color_manual(values = 'black') + 
  labs(title = "Fitted placebo gaps between treatment and synthetic control",
       x = 'Years', y = 'Placebo gaps', color = 'State') +
  theme(plot.title = element_text(size = 12, hjust = .5))

placebo_yr2000 <- placebo_df[31, -40]
pr_gaps <- round(length(which(placebo_yr2000 < as.numeric(placebo_yr2000[3])))
                 / length(placebo_yr2000), 3)

# (e) #

placebo_fun_mspe <- function(s){
  
  placebo_input <- dataprep(
    foo = df,
    predictors = c('lnincome', 'beer', 'age15to24', 'retprice'),
    dependent = 'cigsale',
    unit.variable = 'state',
    time.variable = 'year',
    special.predictors = list(list('cigsale', 1975, 'mean'),
                              list('cigsale', 1980, 'mean'),
                              list('cigsale', 1988, 'mean')),
    treatment.identifier = s,
    controls.identifier = c(placebo_list[placebo_list != s]),
    time.predictors.prior = c(seq(1985, 1988, 1)),
    time.optimize.ssr = c(seq(1970, 1988, 1)),
    time.plot = c(sort(unique(df$year)))
  )
  
  placebo_output <- synth(placebo_input)
  
  return(placebo_output$loss.v)
  
}

placebo_mspe_df <- data.frame(placebo_fun_mspe(1))

i <- 2
while (i < length(placebo_list) + 1){
  
  placebo_mspe_df <- cbind(placebo_mspe_df, placebo_fun_mspe(i))
  i = i + 1
  
}

colnames(placebo_mspe_df) <- placebo_names[-40]
cond_mspe <- placebo_mspe_df[which(placebo_mspe_df <= 20)]
cond_select <- c(colnames(cond_mspe))
cond_df <- cbind(placebo_df[cond_select], sort(unique(df$year)))
colnames(cond_df)[29] <- 'year'

cond_plot_df <- melt(cond_df, id = c('year'))

cond_filter <- cond_plot_df %>% group_by(variable) %>%
  filter(variable == 's3') %>% ungroup()

cond_map <- ggplot() +
  geom_line(aes(x=year, y=value, group = variable), data = cond_plot_df,
            colour = alpha('grey', 0.6)) +
  geom_line(aes(x=year, y=value, colour = variable), data = cond_filter) +
  scale_color_manual(values = 'black') + 
  labs(title = "Fitted placebo gaps between treatment and synthetic control, state MSPE <= 20",
       x = 'Years', y = 'Placebo gaps', color = 'State') +
  theme(plot.title = element_text(size = 12, hjust = .5))

cond_yr2000 <- cond_df[31, -40]
pr_cond <- round(length(which(cond_yr2000 < as.numeric(cond_yr2000[3])))
                 / length(cond_yr2000), 3)

