################
# Markov model #
################

## Model set-up ----

t_names <- c("without_drug", "with_drug")
n_treatments <- length(t_names)

s_names  <- c("Asymptomatic_disease", "Progressive_disease", "Dead")
n_states <- length(s_names)

n_cycles <- 46
Initial_age <- 55

cAsymp <- 500
cDeath <- 1000
cDrug <- 1000
cProg <- 3000
uAsymp <- 0.95
uProg <- 0.75
oDr <- 0.06
cDr <- 0.06
tpDcm <- 0.15

trans_c_matrix <-
  matrix(c(0, 0, cDeath,
           0, 0, cDeath),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

state_c_matrix <-
  matrix(c(cAsymp, cProg, 0,
           cAsymp + cDrug, cProg, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

state_q_matrix <-
  matrix(c(uAsymp, uProg, 0,
           uAsymp, uProg, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))


# Transition probabilities ---- 

# Transition probabilities
p_matrix <- array(data = 0,
                  dim = c(n_states, n_states, n_treatments),
                  dimnames = list(from = s_names,
                                  to = s_names,
                                  t_names))

# Is there a cost on entering the state?
is_trans_cost <-
  matrix(c(0, 0, 0,
           0, 0, 1,
           0, 0, 0),
         byrow = TRUE,
         nrow = n_states,
         dimnames = list(s_names,
                         s_names))


# Store population output for each cycle 

pop <- array(data = NA,
             dim = c(n_states, n_cycles, n_treatments),
             dimnames = list(state = s_names,
                             cycle = NULL,
                             treatment = t_names))

pop["Asymptomatic_disease", cycle = 1, ] <- 1000
pop["Progressive_disease", cycle = 1, ] <- 0
pop["Dead", cycle = 1, ] <- 0

trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_treatments),
               dimnames = list(state = s_names,
                               cycle = NULL,
                               treatment = t_names))

trans[, cycle = 1, ] <- 0


# Sum costs and QALYs for each cycle at a time for each drug 

cycle_empty_array <-
  array(NA,
        dim = c(n_treatments, n_cycles),
        dimnames = list(treatment = t_names,
                        cycle = NULL))

cycle_state_costs <- cycle_trans_costs <- cycle_empty_array
cycle_costs <- cycle_QALYs <- cycle_empty_array
LE <- LYs <- cycle_empty_array
cycle_QALE <- cycle_empty_array

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)


# Time-dependent probability matrix ----

p_matrix_cycle <- function(p_matrix, age, cycle) {
  
  tpProg <- 0.01
  tpDcm <- 0.15
  tpDn_lookup <-
    c("(34,44]" = 0.0017,
      "(44,54]" = 0.0044,
      "(54,64]" = 0.0138,
      "(64,74]" = 0.0379,
      "(74,84]" = 0.0912,
      "(84,100]" = 0.1958)
  Effect <- 0.5
  
  
  age_grp <- cut(age, breaks = c(34,44,54,64,74,84,100))
  
  tpDn <- tpDn_lookup[age_grp]
  
  # Matrix containing transition probabilities for without_drug
  
  p_matrix["Asymptomatic_disease", "Progressive_disease", "without_drug"] <- tpProg*cycle
  
  p_matrix["Asymptomatic_disease", "Dead", "without_drug"] <- tpDn
  
  p_matrix["Asymptomatic_disease", "Asymptomatic_disease", "without_drug"] <- 1 - tpProg*cycle - tpDn
  
  p_matrix["Progressive_disease", "Dead", "without_drug"] <- tpDcm + tpDn
  
  p_matrix["Progressive_disease", "Progressive_disease", "without_drug"] <- 1 - tpDcm - tpDn
  
  p_matrix["Dead", "Dead", "without_drug"] <- 1
  
  # Matrix containing transition probabilities for with_drug
  
  p_matrix["Asymptomatic_disease", "Progressive_disease", "with_drug"] <- tpProg*(1 - Effect)*cycle
  
  p_matrix["Asymptomatic_disease", "Dead", "with_drug"] <- tpDn
  
  p_matrix["Asymptomatic_disease", "Asymptomatic_disease", "with_drug"] <-
    1 - tpProg*(1 - Effect)*cycle - tpDn
  
  p_matrix["Progressive_disease", "Dead", "with_drug"] <- tpDcm + tpDn
  
  p_matrix["Progressive_disease", "Progressive_disease", "with_drug"] <- 1 - tpDcm - tpDn
  
  p_matrix["Dead", "Dead", "with_drug"] <- 1
  
  return(p_matrix)
}



## Run model ----

for (i in 1:n_treatments) {
  
  age <- Initial_age
  
  for (j in 2:n_cycles) {
    
    p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)
    
    pop[, cycle = j, treatment = i] <-
      pop[, cycle = j - 1, treatment = i] %*% p_matrix[, , treatment = i]
    
    trans[, cycle = j, treatment = i] <-
      pop[, cycle = j - 1, treatment = i] %*% (is_trans_cost * p_matrix[, , treatment = i])
    
    age <- age + 1
  }
  
  cycle_state_costs[i, ] <-
    (state_c_matrix[treatment = i, ] %*% pop[, , treatment = i]) * 1/(1 + cDr)^(1:n_cycles - 1)
  
  # Discounting at previous cycle
  cycle_trans_costs[i, ] <-
    (trans_c_matrix[treatment = i, ] %*% trans[, , treatment = i]) * 1/(1 + cDr)^(1:n_cycles - 2)
  
  cycle_costs[i, ] <- cycle_state_costs[i, ] + cycle_trans_costs[i, ]
  
  LE[i, ] <- c(1,1,0) %*% pop[, , treatment = i]
  
  LYs[i, ] <- LE[i, ] * 1/(1 + oDr)^(1:n_cycles - 1)
  
  cycle_QALE[i, ] <-
    state_q_matrix[treatment = i, ] %*% pop[, , treatment = i]
  
  cycle_QALYs[i, ] <- cycle_QALE[i, ] * 1/(1 + oDr)^(1:n_cycles - 1)
  
  total_costs[i] <- sum(cycle_costs[treatment = i, -1])
  total_QALYs[i] <- sum(cycle_QALYs[treatment = i, -1])
}




## Plot results ----

# Incremental costs and QALYs of with_drug vs to without_drug
c_incr <- total_costs["with_drug"] - total_costs["without_drug"]
q_incr <- total_QALYs["with_drug"] - total_QALYs["without_drug"]


# Incremental cost effectiveness ratio 
ICER <- c_incr/q_incr

plot(x = q_incr, y = c_incr,
     xlim = c(0, 1100),
     ylim = c(0, 10e6),
     pch = 16, cex = 1.5,
     xlab = "QALY difference",
     ylab = "Cost difference (£)",
     frame.plot = FALSE)
abline(a = 0, b = 30000) # Willingness-to-pay threshold



#############################################

# Probability Sensitivity Analysis (PSA)


ce_markov <- function(start_pop,
                      p_matrix,
                      state_c_matrix,
                      trans_c_matrix,
                      state_q_matrix,
                      is_trans_cost,
                      n_cycles = 46,
                      init_age = 55,
                      s_names = NULL,
                      t_names = NULL) {
  
  n_states <- length(start_pop)
  n_treat <- dim(p_matrix)[3]
  
  pop <- array(data = NA,
               dim = c(n_states, n_cycles, n_treat),
               dimnames = list(state = s_names,
                               cycle = NULL,
                               treatment = t_names))
  trans <- array(data = NA,
                 dim = c(n_states, n_cycles, n_treat),
                 dimnames = list(state = s_names,
                                 cycle = NULL,
                                 treatment = t_names))
  
  for (i in 1:n_states) {
    pop[i, cycle = 1, ] <- start_pop[i]
  }
  cycle_costs <- array(NA,
                       dim = c(n_treat, n_cycles),
                       dimnames = list(treatment = t_names,
                                       cycle = NULL))
  cycle_QALYs <- array(NA,
                       dim = c(n_treat, n_cycles),
                       dimnames = list(treatment = t_names,
                                       cycle = NULL))
  
  total_costs <- setNames(rep(NA, n_treat), t_names)
  total_QALYs <- setNames(rep(NA, n_treat), t_names)
  
  for (i in 1:n_treat) {
    
    age <- init_age
    
    for (j in 2:n_cycles) {
      
      p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)
      
      #browser()
      # Matrix multiplication
      pop[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% p_matrix[, , treatment = i]
      
      trans[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% (is_trans_cost * p_matrix[, , treatment = i])
      
      age <- age + 1
    }
    
    cycle_state_costs[i, ] <-
      (state_c_matrix[treatment = i, ] %*% pop[, , treatment = i]) * 1/(1 + cDr)^(1:n_cycles - 1)
    
    cycle_trans_costs[i, ] <-
      (trans_c_matrix[treatment = i, ] %*% trans[, , treatment = i]) * 1/(1 + cDr)^(1:n_cycles - 2)
    
    cycle_costs[i, ] <- cycle_state_costs[i, ] + cycle_trans_costs[i, ]
    
    LE[i, ] <- c(1,1,0) %*% pop[, , treatment = i]
    
    LYs[i, ] <- LE[i, ] * 1/(1 + oDr)^(1:n_cycles - 1)
    
    cycle_QALE[i, ] <-
     state_q_matrix[treatment = i, ] %*%  pop[, , treatment = i]
    
    cycle_QALYs[i, ] <- cycle_QALE[i, ] * 1/(1 + oDr)^(1:n_cycles - 1)
    
    
    total_costs[i] <- sum(cycle_costs[treatment = i, -1])
    total_QALYs[i] <- sum(cycle_QALYs[treatment = i, -1])
  }
  
  list(pop = pop,
       cycle_costs = cycle_costs,
       cycle_QALYs = cycle_QALYs,
       total_costs = total_costs,
       total_QALYs = total_QALYs)
}

# Define cost and QALYs as functions

state_c_matrix <- function() {
  matrix(c(rgamma(1, cAsymp, 1), rgamma(1, cProg, 1), 0,
           rgamma(1, cAsymp + cDrug, 1), rgamma(1, cProg, 1), 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))
}

state_q_matrix <- function() {
  matrix(c(runif(1, uAsymp - 0.1, uAsymp + 0.1), runif(1, uProg - 0.1, uProg + 0.1), 0,
           runif(1, uAsymp - 0.1, uAsymp + 0.1), runif(1, uProg - 0.1, uProg + 0.1), 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))
}

trans_c_matrix <- function() {
  matrix(c(0, 0, cDeath,
           0, 0, cDeath),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))
}



## Run PSA analysis ----

costs <- NULL
qalys <- NULL

for (i in 1:100) {
  ce_res <- ce_markov(start_pop = c(1000, 0, 0),
                      p_matrix,
                      state_c_matrix(),
                      trans_c_matrix(),
                      state_q_matrix(),
                      is_trans_cost)
  
  costs <- rbind(costs, ce_res$total_costs)
  qalys <- rbind(qalys, ce_res$total_QALYs)
}


## Plot results ----

# Incremental costs and QALYs of with_drug vs to without_drug
c_incr_psa <- costs[, 2] - costs[, 1]
q_incr_psa <- qalys[, 2] - qalys[, 1]


plot(x = q_incr_psa, y = c_incr_psa,
     xlim = c(0, 1100),
     ylim = c(0, 10e6),
     pch = 16, cex = 1.5,
     xlab = "QALY difference",
     ylab = "Cost difference (£)",
     frame.plot = FALSE)
abline(a = 0, b = 30000) # Willingness-to-pay threshold

points(x = q_incr, y = c_incr, col = "red",
     pch = 16, cex = 1.5)
     