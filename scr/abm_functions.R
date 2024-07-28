# Functions for Running the ABM and cleaning data

#-------------------------Data Cleaning-----------------------------------------
extract_percentages <- function(eeodf_wide, year){
  total <- eeodf_wide[4,12]
  eeodf_perc <- eeodf_wide[-4,]
  eeodf_perc[4:6,] <- eeodf_perc
  for(i in 2:ncol(eeodf_perc[4:6,])) {       
    eeodf_perc[4:6, i] <- eeodf_perc[4:6, i] / eeodf_perc[4:6,]$total_job *100
  }
  eeodf_perc$value <- c(rep("count", 3), rep("perc", 3))
  eeodf_perc <- eeodf_perc[-12]
  eeodf_perc <- pivot_longer(
    eeodf_perc, 
    cols = c("male_white", "male_black", "male_nat_am", "male_asian",
             "male_lat_other_more","fem_white", "fem_black", "fem_nat_am", "fem_asian",
             "fem_lat_other_more"), 
    names_to = "identity", 
    values_to = "n")
  eeodf_perc <- pivot_wider(
    eeodf_perc, 
    id_cols = c("level", "identity"), 
    values_from = "n", 
    names_from = "value")
  eeodf_perc$eth <- NA
  eeodf_perc$gen <- NA
  eeodf_perc$year <- as.factor(year)
  eeodf_perc$identity <- as.character(eeodf_perc$identity)
  perc_total <- eeodf_perc %>% 
    group_by(identity) %>% 
    summarise(perc_total = sum(count)/total*100)
  eeodf_perc <- merge.data.frame(eeodf_perc, perc_total)
  return(eeodf_perc)
}

extract_gen_eth <- function(eeodf_long){
  for(i in 1:length(eeodf_long$identity)) {
    if (str_detect(eeodf_long$identity[i], "fem") == TRUE) {
      eeodf_long$gen[i] <- "woman"
    }
    if (str_detect(eeodf_long$identity[i], "male") == TRUE){
      eeodf_long$gen[i] <- "man"
    }
    if (str_detect(eeodf_long$identity[i], "white") == TRUE) {
      eeodf_long$eth[i] <- "whi"
    }
    if (str_detect(eeodf_long$identity[i], "black") == TRUE){
      eeodf_long$eth[i] <- "bla"
    }
    if (str_detect(eeodf_long$identity[i], "nat_am") == TRUE){
      eeodf_long$eth[i] <- "nat"
    }
    if (str_detect(eeodf_long$identity[i], "lat_other_more") == TRUE){
      eeodf_long$eth[i]<- "oth"
    }
    if (str_detect(eeodf_long$identity[i], "asian") == TRUE){
      eeodf_long$eth[i] <- "asi"
    }
  }
  return(eeodf_long)
}


clean_emp <- function(emp_in){
  emp_out <- as.data.frame(table(emp_in[2:3]))
  emp_out <- emp_out %>% 
    group_by(level) %>% 
    mutate(perc = Freq/sum(Freq)) %>% 
    ungroup() %>%
    group_by(identity) %>% 
    mutate(total_sum =sum(Freq)) %>% 
    ungroup() %>% 
    mutate(perc_total = total_sum/sum(Freq))
  emp_out$gender <- str_split_fixed(emp_out$identity, "_",2)[,1]
  emp_out$ethnicity <- str_split_fixed(emp_out$identity, "_",2)[,2]
  emp_out$group <- "sim"
  return(emp_out)
}
#---------------------Meeting Probabilities-------------------------------------
meeting_probabilities <- function(employee_df, dist_bonus, connection_matrix){
  N <- nrow(employee_df) #number of employees
  prob_meet_matrix <- matrix(
    rnorm(N, 1/N, 0.003), # probability for each emp of meeting another with variance 
    ncol = N, 
    nrow = N,
    dimnames = list(employee_df$ID, employee_df$ID))
  for (row in 1:nrow(prob_meet_matrix)){
    emp_i_ID <- rownames(prob_meet_matrix)[row]
    for (col in (row+1):ncol(prob_meet_matrix)-1){
      emp_j_ID <- colnames(prob_meet_matrix)[col]
      if (employee_df$level[
        employee_df$ID==emp_i_ID] == employee_df$level[employee_df$ID==emp_j_ID]){
        prob_meet_matrix[emp_i_ID,emp_j_ID] <- prob_meet_matrix[emp_i_ID, emp_j_ID]+(
          prob_meet_matrix[emp_i_ID, emp_j_ID]/2)
        prob_meet_matrix[emp_j_ID,emp_i_ID] <- prob_meet_matrix[emp_j_ID, emp_i_ID]+(
          prob_meet_matrix[emp_j_ID, emp_i_ID]/2)
      }
    }
  }
  diag(prob_meet_matrix) <- 0
  #add distance bonus
  graph <- graph_from_matrix_and_df(employee_df, connection_matrix)
  dists <- distances(graph)
  dists[dists == Inf] <- 0
  dists <- round(dists*(dist_bonus^dists), digits = 5)
  prob_meet_matrix <- prob_meet_matrix+dists
  return(prob_meet_matrix)
}

softmax_matrix <- function(prob_meet_matrix){
  for (col in 1:ncol(prob_meet_matrix)){
    prob_meet_matrix[,col] <- exp(prob_meet_matrix[,col])/ sum(exp(prob_meet_matrix[,col]))
  }
  return(prob_meet_matrix)
}

ERG_meeting <- function(prob_meet_matrix, employee_df, ERG_meets_list){
  ERGs <- c("none", "AGN", "BGN", "W@G", "HOLA", "GAIN")
  ERGs_iden <- list(unique(employee_df$identity), 
                    unique(employee_df$identity[employee_df$ethnicity == "asian"]), 
                    unique(employee_df$identity[employee_df$ethnicity == "black"]), 
                    unique(employee_df$identity[employee_df$gender == "fem"]), 
                    unique(employee_df$identity[employee_df$ethnicity=="lat_other_more"]),
                    unique(employee_df$identity[employee_df$ethnicity=="nat_am"]))
  ERG_meet <- sample(ERGs, 1, replace = T)
  ERG_meets_list <- append(ERG_meets_list, ERG_meet)
  partici_iden <- ERGs_iden[ERGs == ERG_meet][[1]]
  partici_IDs <- employee_df$ID[employee_df$identity %in% partici_iden]
  for (row in 1:nrow(prob_meet_matrix)){
    emp_i_ID <- rownames(prob_meet_matrix)[row]
    for (col in (row+1):ncol(prob_meet_matrix)-1){
      emp_j_ID <- colnames(prob_meet_matrix)[col]
      if (emp_i_ID %in% partici_IDs & emp_j_ID %in% partici_IDs){
        prob_meet_matrix[row,col] <- prob_meet_matrix[row, col]*2
        prob_meet_matrix[col,row] <- prob_meet_matrix[col, row]*2
      }
    }
  }
  return(list(prob_meet_matrix, ERG_meets_list))                 
}

#--------------------Connection Determination-----------------------------------
calc_H_i <- function(IH_i, N, N_i){
  w_i = N_i/N
  H_i = (IH_i+(w_i/(1-w_i)))*(1-w_i)
  return(H_i)
}

calc_H_i_multi <- function(emp_ID, employee_df, IH_i){
  N <- nrow(employee_df)
  H_i_m <- matrix(nrow = 3, ncol = 2, dimnames = list(colnames(employee_df[3:5]), c("same", "diff")))
  for (row in 1:nrow(H_i_m)){
    type <- row.names(H_i_m)[row]
    N_i <- nrow(
      employee_df[employee_df[type] == as.character(employee_df[employee_df$ID == emp_ID,type]),])
    H_i_m[type, "same"] <- calc_H_i(IH_i, N, N_i)
    H_i_m[type, "diff"] <- 1-H_i_m[type, "same"]
    if (is.na(H_i_m[type, "same"])){
      print(sprintf("No diversity in %s", type))
      break
    }
  }
  return(H_i_m)
}

connect_determ <- function(emp_i_ID, emp_j_ID, employee_df, IH_i){
  emp_row_i <- employee_df[employee_df$ID==emp_i_ID,3:5]
  emp_row_j <- employee_df[employee_df$ID==emp_j_ID,3:5]
  hi_i <- calc_H_i_multi(emp_i_ID, employee_df, IH_i)
  hi_j <- calc_H_i_multi(emp_j_ID, employee_df, IH_i)
  check <- t(as.matrix(emp_row_i == emp_row_j))
  check <- cbind(check, ifelse(check == F, T, F))
  prob_i <- sum(hi_i[check])/3
  prob_j <- sum(hi_j[check])/3
  
  if (is.na(sum(hi_i))){
    connect <- NA
    return(connect)
  }
  connect_i <- sample(c(0,1), 1, replace = T, prob = c(1-prob_i, prob_i))
  connect_j <- sample(c(0,1), 1, replace = T, prob = c(1-prob_j, prob_j))
  if (connect_i + connect_j == 2){
    connect <- 1
  }
  else{
    connect <- 0
  }
  return(connect)
}

exp_decay <- function(connection_matrix){
  for (row in 1:nrow(connection_matrix)){
    for (col in (row+1):ncol(connection_matrix)-1){
      connection_matrix[row, 
                        col] <- connection_matrix[col, 
                                                  row] <- (connection_matrix[row, col])*(1-0.05)^1
    }
  }
  return(connection_matrix)
}

#------------------------Hiring Functions---------------------------------------
iden_proportion <- function(employee_df){
  hire_props <- employee_df[2:3] %>% 
    count(level, identity) %>% 
    group_by(level) %>% 
    mutate(prob = n/sum(n)) %>% 
    ungroup() 
  hire_props <- pivot_wider(hire_props[c(1:2,4)], names_from = identity, values_from = prob)
  hire_props <- hire_props %>% remove_rownames %>% column_to_rownames(var="level")
  hire_props[is.na(hire_props)] <- 0.001
  idens <- c("male_white", "male_black", "male_asian", "male_nat_am", "male_lat_other_more", 
             "fem_white", "fem_black", "fem_asian", "fem_nat_am", "fem_lat_other_more")
  for (type in idens){
    if (!type %in% colnames(hire_props)){
      vals <- c(0.01, 0.01, 0.01)
      hire_props <- cbind(hire_props, vals)
      colnames(hire_props)[ncol(hire_props)] <- type
    }
  }
  return(hire_props)
}

iden_probability <- function(employee_df){
  hire_props <- iden_proportion(employee_df)
  hire_probs <- hire_props
  for (row in rownames(hire_props)){
    max_prop <- (which.max(hire_props[row,]))
    for (col in colnames(hire_props)){
      if (col != max_prop){
        prop <- hire_props[row,col]
        hire_probs[row,col] <- rtruncnorm(1, a = (prop-0.001), mean = (prop+0.001), sd = 0.05)
        if (hire_probs[row, col] < 0.01){
          hire_probs[row, col] <- 0.01
        }
      }
    }
    if (1-rowSums(hire_probs[row,-max_prop])>0){
      hire_probs[row, max_prop] <- 1-rowSums(hire_probs[row,-max_prop])
    }
    else {
      hire_probs[row, max_prop] <- 0
    }
  }
  return(hire_probs)
}

new_hire_fun <- function(day, leaving_df, employee_df){
  N <- nrow(employee_df)
  hire_probs <- iden_probability(employee_df)
  new_hire <- data.frame(
    ID = paste0("e", as.character(((N+1)+(nrow(leaving_df)*((day/30)-1))):(
      ((N+1)+(nrow(leaving_df)*((day/30)-1)))+nrow(leaving_df)-1))),
    level = NA, identity = NA, gender = NA, ethnicity = NA)
  for (row in 1:nrow(leaving_df)){
    new_hire$level[row] <- leaving_df[row,"level"]
    new_hire$identity[row] <- sample(colnames(hire_probs), 1, 
                                     prob = hire_probs[leaving_df[row,"level"],])
    new_hire$gender[row] <- str_split_fixed(new_hire$identity[row], "_",2)[,1]
    new_hire$ethnicity[row] <- str_split_fixed(new_hire$identity[row], "_",2)[,2]
  }
  return(new_hire)
}

update_connection_matrix <- function(connection_matrix, leaving_df, new_hire){
  connection_matrix <- connection_matrix[!rownames(connection_matrix)%in%leaving_df$ID,]
  connection_matrix <- connection_matrix[,!colnames(connection_matrix)%in%leaving_df$ID] 
  new_rows <- matrix(0, nrow = nrow(new_hire), ncol = ncol(connection_matrix),
                     dimnames = list(new_hire$ID, colnames(connection_matrix)))
  connection_matrix <- rbind(connection_matrix, new_rows)
  new_cols <- matrix(0, nrow = nrow(connection_matrix), ncol = nrow(new_hire), 
                     dimnames = list(rownames(connection_matrix), new_hire$ID))
  connection_matrix <- cbind(connection_matrix, new_cols)
  return(connection_matrix)
}


#---------------------Calculations etc. ----------------------------------------
homophily_calc <- function(Homophily_all, employee_df, connection_matrix, day){
  N <- nrow(employee_df)
  n_iden <- length(unique(employee_df$identity))
  n_eth <- length(unique(employee_df$ethnicity))
  Homophily <- data.frame(day = NA,
                          type=c(rep("identity",n_iden), rep("gender",2), rep("ethnicity",n_eth)),
                          n = 0, type_value = NA, si = 0, di = 0, wi = 0, Hi = 0, IHi = 0)
  Homophily$type_value[Homophily$type == "identity"] <- unique(employee_df$identity)
  Homophily$type_value[Homophily$type == "gender"] <- c("male", "fem")
  Homophily$type_value[Homophily$type == "ethnicity"] <- unique(employee_df$ethnicity) 
  
  for (type in unique(Homophily$type)){
    for (row in rownames(connection_matrix)){
      type_value <- as.character(employee_df[employee_df$ID == row, type])
      si <- 0
      di <- 0
      for (col in colnames(connection_matrix)){
        if (connection_matrix[row,col]>0){
          if (type_value == employee_df[employee_df$ID == col, type]){
            si <- si+1
          }
          else{
            di <- di+1
          }
        }
      }
      Homophily$si[Homophily$type == type & Homophily$type_value == type_value] <- 
        Homophily$si[Homophily$type == type & Homophily$type_value == type_value] + si
      Homophily$di[Homophily$type == type & Homophily$type_value == type_value] <- 
        Homophily$di[Homophily$type == type & Homophily$type_value == type_value] + di
      Homophily$n[Homophily$type == type & Homophily$type_value == type_value] <- 
        Homophily$n[Homophily$type == type & Homophily$type_value == type_value] + 1
    }
  }
  for (row in 1:nrow(Homophily)){
    Homophily$wi[row] <- Homophily$n[row]/N
    Homophily$si[row] <- Homophily$si[row]/Homophily$n[row]
    Homophily$di[row] <- Homophily$di[row]/Homophily$n[row]
    Homophily$Hi[row] <- Homophily$si[row]/(Homophily$si[row]+Homophily$di[row])
    Homophily$IHi[row] <- (Homophily$Hi[row]-Homophily$wi[row])/(1-Homophily$wi[row])
  }
  Homophily$day <- day
  Homophily_all <- rbind(Homophily_all, Homophily)
  return(Homophily_all)
}

graph_from_matrix_and_df <- function(v_attribute_df, connection_matrix){
  connection_matrix[connection_matrix>0] <- 1
  edges <- as_edgelist(graph_from_adjacency_matrix(connection_matrix, mode = "undirected"))
  graph <- graph_from_data_frame(d=edges, v = v_attribute_df, directed = F)
  return(graph)
}

#------------------------Complete ABM function----------------------------------

abm_DEI_function <- function(employee_df, 
                             days_run, 
                             IH_i, 
                             dist_bonus){
  N = nrow(employee_df)
  # Empty matrix, for storing connections
  connection_matrix <- matrix(0,
                              ncol = N, 
                              nrow = N,
                              dimnames = list(employee_df$ID, employee_df$ID))
  # Empty data frame for storing employees that leave 
  left_df <- data.frame(ID = as.character(), 
                        level = as.character(), 
                        identity= as.character(), 
                        gender= as.character(), 
                        ethnicity= as.character(), 
                        sum= as.numeric(), 
                        day = as.numeric())
  # Empty list for storing which ERGs meet 
  ERG_meets_list <- list()
  # Empty data frame for storing average homophily index values for all groups
  Homophily_all <- data.frame(day = as.numeric(),type = as.character(),n = as.numeric(), 
                              type_value = as.character(), si = as.numeric(), di = as.numeric(), 
                              wi = as.numeric(), Hi = as.numeric(), IHi = as.numeric())
  
  # Start the loop, let it run for days_run
  for (day in 1:days_run){
    # Each day: define meeting probabilities (end with softmax funciton)
    prob_meet_matrix <- meeting_probabilities(employee_df, dist_bonus, connection_matrix)
    # Every five days: an ERG meets, and the attending employees get additional meeting bonuses
    if ((day/5)%% 1 == 0){
      ERG_output <- ERG_meeting(prob_meet_matrix, employee_df, ERG_meets_list)
      prob_meet_matrix <- ERG_output[[1]]
      ERG_meets_list <- ERG_output[[2]]
    }
    prob_meet_matrix <- softmax_matrix(prob_meet_matrix)
    # Each day: every employee meets another, it is determined if they connect (with decay)
    for (emp_i in employee_df$ID){
      emp_j <- sample(employee_df$ID, 1, replace = T, prob = prob_meet_matrix[,emp_i])
      conn <- connect_determ(emp_i, emp_j, employee_df, IH_i)
      if (is.na(conn)){
        print(sprintf("ABM has stopped at day %s and emp %s", day, emp_i))
        status <- "stop"
        break
      }
      status <- "continue"
      connection_matrix[emp_j,
                        emp_i] <- connection_matrix[emp_i,
                                                    emp_j] <- connection_matrix[emp_j,emp_i]+conn
    }
    if (status == "stop"){
      write_csv(employee_df, sprintf("out/employee_list_day%s_IH%s.csv", day, IH_i))
      write.csv(connection_matrix, file = sprintf("out/conn_matrix_day%s_IH%s.csv", day, IH_i))
      break
    }
    connection_matrix <- exp_decay(connection_matrix)
    # Once a year: the employee df, with its connection matrix is saved to csv files. 
    if ((day/365)%%1 == 0){
      write_csv(employee_df, sprintf("out/employee_list_day%s_IH%s.csv", day, IH_i))
      write.csv(connection_matrix, file = sprintf("out/conn_matrix_day%s_IH%s.csv", day, IH_i))
    }
    # Once a month: determine who leaves (5% of the employees) and who is hired
    if ((day/30)%%1 == 0){
      # Determine exits
      print(sprintf("day %s", day))
      Homophily_all <- homophily_calc(Homophily_all, employee_df, connection_matrix, day)
      colsums <- data.frame("ID" = as.vector(names(colSums(connection_matrix))), 
                            "sum" = as.vector(colSums(connection_matrix)))
      employee_df <- merge(employee_df, colsums)
      leaving_df <- top_frac(employee_df, -0.05)
      leaving_df$day <- day
      left_df <- rbind(left_df, leaving_df)
      # Determine hires
      new_hire <- new_hire_fun(day, leaving_df, employee_df)
      employee_df <- employee_df[!employee_df$ID%in%leaving_df$ID,-6]
      employee_df <- rbind(employee_df, new_hire)
      connection_matrix <- update_connection_matrix(connection_matrix, leaving_df, new_hire)
    }
  }
  write.csv(ERG_meets_list, sprintf("out/output_ERG_IH%s.csv", IH_i))
  write.csv(left_df, sprintf("out/output_left_IH%s.csv", IH_i))
  write.csv(Homophily_all, sprintf("out/output_hom_IH%s.csv", IH_i))
  return(print(sprintf("The ABM has run for %s days", day))
)
}