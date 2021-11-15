#Load packages 
library(tidyverse)
library(data.table)
library(qcc)


# definations for lambda sequence
lambda_seq_01 <- c(0.1)
lambda_seq_03  <- c(0.1, 0.5, 0.9)
lambda_seq_all <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
lam_data <- readRDS("lam_data.RDS") 



get_nsigma <- function(l) {
  
  lam_data[which(l >= lam_data$lam & l <  lam_data$lam +0.09 ),2]
  
}

# Sample Mean SD tibble for simulation result generation
sample_mean_sd_RI <-readRDS("sample_mean_sd.RDS")  
sample_TE <-readRDS("sample_TE.RDS")  



# Programming related functions -----------------------------------------------------

group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}


unnest_dt_2 <- function(dt, col, id){
  stopifnot(is.data.table(dt))
  by <- substitute(id)
  col <- substitute(unlist(col, recursive = FALSE))
  dt[, eval(col), by = eval(by)]
}




# simulation related functions ------------------------------------------------------

#Generates normal distributed for given days * n with given mean and SD 
gen_sim_norm_data <- function(days, n, res_mean, res_sd) {
  
  
  tibble(day = rep(1:days, each = n),
         id_day = rep(1:n, each = 1, len = n * days) , 
         result_num = rnorm(n * days, res_mean,  res_sd )  )
  
  
}

#Generates test name based simulation data with adding TEa values with virtual days
generate_sim_data <- function(days, n, source_df, TE_df) {

  generated_data <- source_df %>%
    nest(ref_data = c(mean_res, sd_res))  |> 
    mutate(ri = map(ref_data, ~(gen_sim_norm_data(days,n, .x$mean_res, .x$sd_res)) )) |> 
    select(test_name, ri) |> 
    unnest(cols = c(ri)) |> 
    left_join(source_df, by = "test_name") |> 
    left_join(select(TE_df, test_name,  tea), by = "test_name") |> 
    select(test_name, tea, low_trunc, high_trunc , mean_res, sd_res, everything() )
  
}



add_error <- function(df, tea, tea_n, error_point, error_direction) {
  
  
  te_seq <-  seq(0, tea, by = (tea  - 0) / (tea_n-1))
  
  # Create data.table according to SE sequence 
  temp <- data.table(SE = te_seq)
  
  # Add results into each SE column
  temp[, result_data := list(data.table(df))]
  
  temp_2 <- unnest_dt_2(temp, result_data, list(SE))

  temp_2[,result_final := fifelse(id_day < error_point,result_num,result_num  +  result_num * SE/100 * (error_direction))] |> 
    tibble() |> 
    mutate(tea = tea, error_point = error_point, error_direction = error_direction) |> 
    select(test_name, tea, SE, error_point, error_direction, everything()) |> 
    tibble()
  
  
  
}



# EWMA related functions --------------------------------------------------

calculate_ewma_multi_lambda_multi_day <- function(ewma_data, lambda_list = c(0.1, 0.9)) {
  
  temp <- tibble(lambda =  lambda_list) |>  
    mutate(ewma_result = map(lambda_list, ~calculate_ewma_multi_day(ewma_data, .x))) |> 
    unnest(cols = c(ewma_result))
  
  temp 
}


calculate_ewma_multi_day<- function(ewma_data, lambda  = 0.1) {
  
  ewma_data |>  
    group_by(day, mean_res, sd_res) |> 
    nest() |>  
    mutate(ewma_result = map(data, ~calculate_ewma(.x, lam = lambda, mean_res, sd_res, day) ))  |>  
    select(-data) |> 
    unnest(cols = ewma_result)
  
  
}


calculate_ewma <- function(df,  lam, mean_test, sd_test, cur_SE = NULL, cur_day = NULL, test_name = NULL)    {
  
  
  if (FALSE)  print(paste0(test_name, " Lambda: ", lam,  " Day :", cur_day)) # debuging
  
  ewma_result <- ewma(df$result_final, plot = FALSE, center = mean_test, std.dev = sd_test, sizes = 1, lambda = lam, 
                      nsigmas = get_nsigma(lam) )
  
  x <-    ewma_result[["x"]]
  y <-    ewma_result[["y"]]
  LCL <-    ewma_result[["limits"]][, "LCL"]
  UCL <-    ewma_result[["limits"]][, "UCL"]
  
  old_id_day <- df$id_day
  
  
  error_high <- as.integer( y > UCL)
  error_low <- as.integer( y < LCL)
  
  ewma_data <- data.table(old_id_day, y, LCL, UCL,error_low, error_high)
  setnames(ewma_data, "old_id_day", "id_day")
  ewma_data
 
  
} 




# Calculates ANPed MNPed --------------------------------------------------
calculate_ewma_stats <- function(ewma_result_df) {
  
  
  temp_df <-  ewma_result_df |> 
    mutate(error = if_else(error_direction == 1, as.logical(error_high),as.logical(error_low)) ) |> 
    select(-error_low, -error_high) |> 
    filter(id_day >= error_point) |> 
    filter(error) 

  if (nrow(temp_df) == 0) {
    
    stats <- ewma_result %>%
      filter(id_day >= error_point) %>% 
      select(test_name, tea, SE, error_point, error_direction, lambda) %>% 
      mutate(days_with_error  = 0, ANPed  = NA_real_, MNPed  = NA_real_, min_min_err  = NA_real_, max_min_err = NA_real_) %>% 
      distinct()
    
    # test_name   tea    SE error_point error_direction lambda days_with_error ANPed MNPed min_min_err max_min_err
    
    
  } else {
    
    
    stat_daily <- temp_df %>% 
      group_by(test_name, tea, SE, error_point, error_direction, lambda, day) |> 
      filter( row_number() == min( row_number() )) %>% 
      mutate(min_error = id_day - error_point ) %>% 
      ungroup()
    
    stat_all <- stat_daily %>%
      group_by(test_name, tea, SE, error_point, error_direction, lambda) |> 
      summarise(
        days_with_error = n(), ANPed = mean(min_error), MNPed = median(min_error),
        min_min_err = min(min_error), max_min_err = max(min_error)
      ) %>%
      ungroup()
    
  }
  
  
  stat_all
}





# START HERE --------------------------------------------------------------------

main <- function()  {
  
  
  example_distrubiton_data <- sample_mean_sd_RI |>
    rename(low_trunc = RI_low, high_trunc = RI_high)
  
  
  # Generate simulation data
  sim_data <- generate_sim_data(days = 200, n = 250, example_distrubiton_data, sample_TE)  
  
  #Real patient data must be same with sim_data
  sim_data
  
  
  ast_data <- sim_data |> filter(test_name   == "AST")
  
  
  # ADD SE ERROR WITH percentage of TEa 
  # tea_n: Defines SE count
  
  error_added_df <- sim_data |> 
    group_by(test_name, tea ) |> 
    mutate(id = cur_group_id()) |> 
    group_by(id, tea) |> 
    nest() |> 
    ungroup() |> 
    mutate(err_added = map2(data, tea, ~add_error(.x, .y,tea_n = 5, error_point = 30, error_direction = 1))) |> 
    select(err_added) |> 
    unnest(cols = c(err_added))
  
  
  # TRANSFORMATON CAN BE APPLIED AT THAT POINT
  
  # TRUNCATE RESULTS
  truncated_df <- error_added_df |> 
    filter(result_final >= low_trunc & result_final <= high_trunc )
  
  # TRANSFORMATON CAN ALSO BE APPLIED AT THIS POINT
  
  
  
  # EWMA ANALYSIS for multiple test with multiple days with multiple lambda
  
  ewma_result <- truncated_df |> 
    group_by(test_name, tea, SE, error_point, error_direction) |> 
    nest()  |> 
    mutate(ewma_result = map(data,  ~calculate_ewma_multi_lambda_multi_day(.x, lambda_seq_01))) |> 
    select(-data) |> 
    unnest(cols = c(ewma_result)) |> 
    ungroup()
  
  
  # MNPed CALCULATION
  
  
  stats <- calculate_ewma_stats(ewma_result)
  
  
  stats |> 
    arrange(test_name,   tea,    SE, error_point, error_direction, lambda)
  
}
















