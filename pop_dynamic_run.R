library(dplyr)
library(readr)
library(DirichletReg)

run_simulation <- function(intervention_prob = 1,
                           repetitions = 1,
                           antibiotic_pressure = 0.006,
                           study = "NORM",
                           generations = 120,
                           fa = 0,
                           m = 0.01,
                           sigma_f = 0.1,
                           v = 0,
                           seed1 = -99,
                           seed2 = -99){


  if(seed1 != -99) set.seed(seed1)

  if(study == "NORM"){
    strain_amr_data <- study_parse_strain_amr_file("~/data/mmc2_full.csv", "ST", "CTX-M", 20, "NORM")
  }
  if(study == "BSAC_NORM"){
    strain_amr_data_NORM <- study_parse_strain_amr_file("~/data/mmc2_full.csv", "ST", "CTX-M", 20, "NORM")
    strain_amr_data <- study_parse_strain_amr_file("~/data/bsac_ctxm_ed.csv", "ST", "CTX-M", 20, "BSAC")
  }
  if(study == "BSAC"){
    strain_amr_data <- study_parse_strain_amr_file("~/data/bsac_ctxm_ed.csv", "ST", "CTX-M", 20, "BSAC")
  }

  genotype <- names(strain_amr_data$strains)

  pop_names <- factor("NORM")
  pop_size <- 50000
  kappa <- pop_size

  migration_pressure <- 0

  populations <- length(pop_size)
  loci_names <- c("CTXM", paste0("L", 1:20))
  IFloci <- loci_names
  loci_length <- length(loci_names)
  genotype_names <- genotype
  genotype_length <- length(genotype_names)

  e <- strain_amr_data$strains
  e <- e/sum(e)

  e_pre <- matrix(0, nrow=genotype_length, ncol=loci_length)
  rownames(e_pre) <- genotype_names
  colnames(e_pre) <- loci_names
  
  alpha_vec <- 10000*rt(loci_length * genotype_length, df=2) - 100
  alphas <- matrix(pmax(alpha_vec, 0.001), nrow=genotype_length, ncol=loci_length)

  e_pre_1 <- e_pre
  e_pre_2 <- e_pre
  for(i in 1:genotype_length){
    e_pre_1[i,] <- DirichletReg::rdirichlet(1, c(alphas[i,]))
    e_pre_2[i,] <- e_pre_1[i,]
  }
  e_pre_1[,1] <- strain_amr_data$amr[,2] + 0.01
  e_pre_2[,1] <- (strain_amr_data$amr[,2] + 0.01) * 0.1

  # Initialize state
  if(study == "BSAC_NORM"){
    e_NORM <- strain_amr_data_NORM$strains
    genotypes_NORM <- names(strain_amr_data_NORM$strains)
    e_mixed <- e
    e_mixed[genotype == "131_A"] <- e_NORM[genotypes_NORM == "131_A"]
    e_mixed[genotype == "131_B"] <- e_NORM[genotypes_NORM == "131_B"]
    e_mixed[genotype == "131_C1"] <- e_NORM[genotypes_NORM == "131_C1"]
    e_mixed[genotype == "131_C2"] <- e_NORM[genotypes_NORM == "131_C2"]
    e_mixed <- e_mixed / sum(e_mixed)

    state <- init_population(pop_names, pop_size, genotype_names, loci_names, kappa, e_mixed, e_pre_1)
  }else{
    state <- init_population(pop_names, pop_size, genotype_names, loci_names, kappa, e, e_pre_1)
 }

  state_stable <- state


  e_pre <- colSums(state %>% dplyr::select(IFloci)) / pop_size

  state <- init_population(pop_names, pop_size, genotype_names, loci_names, kappa, e, e_pre_2)

  print("State Initialised")
  e_init <- e
  e_init[genotype_names == "69" | genotype_names == "131_A" | genotype_names == "131_B" | genotype_names == "131_C1" | genotype_names == "131_C2"] <- 0.1 * e_init[genotype_names == "69" | genotype_names == "131_A" | genotype_names == "131_B" | genotype_names == "131_C1" | genotype_names == "131_C2"]
  e_init <- e_init / sum(e_init)
  
  genotype_index = c()
  for( i in genotype){
    genotype_index <- c(genotype_index, sample(which(state$genotype == i), size=round(e_init[i] * pop_size), replace = TRUE))
  }
  
  state <- state[genotype_index,]

  state_init <- state
  total_rows <- populations * genotype_length * generations * repetitions
  total_rows_pop <- populations * generations * repetitions
  summaries_list <- summarize_population(state, loci_names, pop_names, genotype_names, igeneration = 1, irepetition = 1, ab_pressure = antibiotic_pressure[1])
  summarized_state_genotype_init <- summaries_list[[1]]

  
  print("Start simulation ...")

  if(seed2 != -99) set.seed(seed2)
  for(irepetition in 1:repetitions){
    print(paste("Repetition :", irepetition))
    for(ab_pressure in antibiotic_pressure){

      interventions <- generate_interventions(generations,
                                              loci_length,
                                              intervention_prob,
                                              loci_names,
                                              pop_names)
  
      for(igeneration in 1:generations){

        if(igeneration %% 10 == 0) print(round(igeneration / generations, 2))
        
        if(igeneration == 1){
  
          state <- state_init
          if(irepetition == 1 & ab_pressure == antibiotic_pressure[1]){
            summarized_state_genotype <- summarized_state_genotype_init

          } else {
            summarized_state_genotype_init$repetition <- irepetition
            summarized_state_genotype_init$antibiotic_pressure <- ab_pressure
            summarized_state_genotype <- dplyr::bind_rows(summarized_state_genotype, summarized_state_genotype_init)
          }
        } else {
  
          state <- evolve_generation(state, dplyr::filter(interventions, generation == igeneration),
                                     e_pre, kappa, m, v, sigma_f,
                                     loci_names, igeneration, pop_names, genotype_names,
                                     migration_pressure, ab_pressure, state_stable, fa)
          summaries_list <- summarize_population(state, loci_names, pop_names, genotype_names, igeneration, irepetition, ab_pressure)
          summarized_state_genotype <- dplyr::bind_rows(summarized_state_genotype, summaries_list[[1]])
        }
      }
    }
  }
  return(list(summarized_state_genotype=summarized_state_genotype))
}

###########

study_parse_strain_amr_file <- function(filename, strain_column_name, amr_column_name, n_top, study){

  if(study == "NORM"){
    strain_data <- readr::read_delim(filename, delim=";")
    strain_data <- strain_data %>% dplyr::filter(year==2017)
    strain_data$CC131_clades[is.na(strain_data$CC131_clades)] <- ""
    strain_data$subglades <- paste0(strain_data$ST, strain_data$CC131_clades)
    
    strain_data$subglades[strain_data$subglades == "131A"] <- "131_A"
    strain_data$subglades[strain_data$subglades == "131B"] <- "131_B"
    strain_data$subglades[strain_data$subglades == "131C1"] <- "131_C1"
    strain_data$subglades[strain_data$subglades == "131C2"] <- "131_C2"
    strain_column_name <- "subglades"
    
    strain_table <- table(strain_data[strain_column_name])
    
    strain_sorted_tmp <- sort.int(strain_table, decreasing = TRUE, index.return= TRUE)
    strain_sorted <- strain_sorted_tmp$x
    sorted_index <- strain_sorted_tmp$ix
    
    n_top_strains <- strain_sorted[1:n_top]
    n_top_strain_names <- names(n_top_strains)
    
    other_strains <- strain_sorted[(n_top+1):length(strain_sorted)]
    other_strains_total <- sum(other_strains)
    other_strain_names <- names(other_strains)
    n_top_strains['STX'] <- other_strains_total
    
    strain_table_freq <- n_top_strains / sum(n_top_strains)
    
    strain_amr <- strain_data %>%
      dplyr::select(amr_column_name) %>%
      mutate_all(~ replace(., . == "CTX-M", 1)) %>%
      mutate_all(~ replace(., is.na(.), 0)) %>%
      tibble::deframe()
    strains <- strain_data[,strain_column_name] %>% tibble::deframe()
    
    amr_table <- table(strains, strain_amr)
    
    amr_table_top <- amr_table[n_top_strain_names,]
    
    amr_table_other <- amr_table[(row.names(amr_table) %in% other_strain_names),]
    
    amr_table_other_amr <- colSums(amr_table_other)
    amr_table_top <- rbind(amr_table_top, amr_table_other_amr)
    
    rownames(amr_table_top)[n_top+1] <- "STX"
    amr_table_top <- amr_table_top / rowSums(amr_table_top)
    
    return(list(strains=strain_table_freq, amr=amr_table_top))
  }
  if(study == "BSAC"){
    strain_data <- read.csv(filename, sep=";")
    strain_data$`CTX-M` <- (strain_data$blaCTX.M_any == "yes") * 1
    strain_data <- strain_data %>% dplyr::filter(year==2015 | year==2016 | year==2017)
    strain_data$ST131_clade[is.na(strain_data$X131_clade)] <- ""
    strain_data$subglades <- paste0(strain_data$ST, strain_data$X131_clade)
    
    strain_data$subglades[strain_data$subglades == "131A"] <- "131_A"
    strain_data$subglades[strain_data$subglades == "131B"] <- "131_B"
    strain_data$subglades[strain_data$subglades == "131C1"] <- "131_C1"
    strain_data$subglades[strain_data$subglades == "131C2"] <- "131_C2"
    strain_column_name <- "subglades"
    
    strain_table <- table(strain_data[strain_column_name])
    
    strain_sorted_tmp <- sort.int(strain_table, decreasing = TRUE, index.return= TRUE)
    strain_sorted <- strain_sorted_tmp$x
    sorted_index <- strain_sorted_tmp$ix
    
    n_top_strains <- strain_sorted[1:n_top]
    n_top_strain_names <- names(n_top_strains)
    
    other_strains <- strain_sorted[(n_top+1):length(strain_sorted)]
    other_strains_total <- sum(other_strains)
    other_strain_names <- names(other_strains)
    n_top_strains['STX'] <- other_strains_total
    
    strain_table_freq <- n_top_strains / sum(n_top_strains)
    
    strain_amr <- strain_data %>%
      dplyr::select(amr_column_name) %>%
      tibble::deframe()
    
    strains <- tibble(strain_data)[,strain_column_name] %>% tibble::deframe()
    
    amr_table <- table(strains, strain_amr)
    
    amr_table_top <- amr_table[n_top_strain_names,]
    
    amr_table_other <- amr_table[(row.names(amr_table) %in% other_strain_names),]
    
    amr_table_other_amr <- colSums(amr_table_other)
    amr_table_top <- rbind(amr_table_top, amr_table_other_amr)
    
    rownames(amr_table_top)[n_top+1] <- "STX"
    amr_table_top <- amr_table_top / rowSums(amr_table_top)
    
    return(list(strains=strain_table_freq, amr=amr_table_top))
    
  }

}


summarize_population <- function(state, IFloci, population, genotype, igeneration, irepetition, ab_pressure){
  summary_table_genotype_empty <- tidyr::expand_grid(population,genotype)
  pop_count <- state %>% dplyr::group_by(population)  %>% dplyr::tally() %>% dplyr::ungroup()
  pop_and_genotype_count <- state %>% dplyr::group_by(population, genotype)  %>% dplyr::tally() %>% dplyr::ungroup()

    pop_ab_prop <- state %>% dplyr::group_by(population) %>% dplyr::summarize_at(IFloci, mean) %>% dplyr::ungroup()
    pop_and_genotype_ab_prop <- state %>% dplyr::group_by(population, genotype) %>% dplyr::summarize_at(IFloci, mean) %>% dplyr::ungroup()

  summary_table_genotype <- pop_and_genotype_count %>%
    dplyr::inner_join(pop_count, by = "population", name = "total", suffix = c("", "_pop")) %>%
    dplyr::mutate(genotype_prop = n / n_pop) %>%
    dplyr::inner_join(pop_and_genotype_ab_prop, by = c("population", "genotype")) %>%
    dplyr::select(-c("n", "n_pop"))  %>%
    dplyr::inner_join(summary_table_genotype_empty, by = c("population", "genotype"))

  summary_table_genotype <-summary_table_genotype_empty %>%
    dplyr::left_join(summary_table_genotype, by = c("population", "genotype")) %>%
    dplyr::mutate(generation = igeneration, repetition = irepetition, antibiotic_pressure = ab_pressure) # %>%

  summary_table_population <- pop_count %>%
    dplyr::inner_join(pop_ab_prop, by = c("population"), suffix = c("_genotype", "_pop")) %>%
    dplyr::mutate(generation = igeneration, repetition = irepetition, antibiotic_pressure = ab_pressure)

  summary_table = list(summary_table_genotype)
  return(summary_table)

}


###########

evolve_generation <- function(state,
                              interventions,
                              e_pre,
                              kappa,
                              m,
                              v,
                              sigma_f,
                              IFloci,
                              igeneration,
                              pop_names,
                              genotypes,
                              migration_pressure,
                              antibiotic_pressure,
                              state_stable, 
                              fa) {


  N_t <- nrow(state)

  # Antibiotic intervention
  pop_counts <- table(state$population)
  interventions_long <- interventions %>% dplyr::slice(rep(1:dplyr::n(), times = as.numeric(pop_counts), each=TRUE)) %>% dplyr::select("CTXM")
  potential_destroy <-  (state %>% dplyr::select("CTXM") == 0)
  destroy <- rep(FALSE, N_t)
  for(j in pop_names){
    destroy[state$population == j] <- (rbinom(as.numeric(pop_counts[j]), 1, antibiotic_pressure) > 0) & potential_destroy[state$population == j]
  }
  state <- state %>% dplyr::filter(!destroy)

  # NFDS:
  new_state_index <- c()
  N_t_pop <- 0
  pop_counts <- table(state$population)

  for(j in pop_names){
    N_t <- pop_counts[j]

    f_t <- colSums(state %>% dplyr::filter(population == j) %>% dplyr::select(IFloci)) / N_t
    diff_vector <- e_pre - f_t
    present <- state %>% dplyr::filter(population == j) %>% dplyr::select(IFloci)
    pi_it <- rep(0, N_t)
    for(k in 1:length(IFloci)){
      pi_it[present[,k] == 1] = pi_it[present[,k] == 1] + diff_vector[k]
    }
    r_i <- ((state$genotype == "69") | (state$genotype == "131_A") | (state$genotype == "131_B") | (state$genotype == "131_C")) * fa

    lambda <- (kappa / N_t) * (1 + r_i) * (1 - m) * (1 - v) * (1 + sigma_f) ^ pi_it

    X_it <- rpois(N_t, lambda)
    index = which(state$population == j)

    new_state_index <- c(new_state_index, rep(index, X_it))
  }

  state <- state[new_state_index,]

  N_stable <- nrow(state_stable)
  N_state <- nrow(state)
  state <- rbind(state, state_stable[sample(kappa, m * kappa, replace=TRUE), ])
  return(state)

}

###########

init_population <- function(pop_names, pop_size, genotype_names, IFloci, kappa, e, e_pre){
  # Distribution of genotypes in the sample population
  # TODO : Fix the simple solution
  X_pre <- c(rmultinom(1, pop_size[1], c(e)))
  if(length(pop_size) > 1){
    for(i in 2:length(pop_size)){
      X_pre <- c(X_pre,c(rmultinom(1, pop_size[i], c(e))))
    }
  }
  # Long format for genotypes
  X_t <- rep(rep(genotype_names, length(pop_size)), times = X_pre)

  AMR_unif <- matrix(runif(sum(pop_size) * length(IFloci)),
                     nrow = sum(pop_size),
                     ncol = length(IFloci),
                     dimnames = list(c(), IFloci))
  e_pre_matrix <- e_pre[X_t,]

  # Antibiotic resistance indicator matrix
  AMR_it <- dplyr::as_tibble((e_pre_matrix - AMR_unif) > 0) * 1

  state <- dplyr::bind_cols(dplyr::tibble(population = factor(rep(pop_names, times = pop_size), levels = pop_names),
                                          genotype = as.factor(X_t), AMR_it))

  return(state)
}



generate_interventions <- function(generations,
                                   loci_length,
                                   intervention_prob,
                                   loci_names,
                                   pop_names){
  
  interventions <- dplyr::as_tibble(matrix(rbinom(generations * loci_length,
                                                  size = 1,
                                                  prob = intervention_prob[1]),
                                           nrow = generations,
                                           ncol = loci_length, dimnames = list(c(), loci_names)))
  
  interventions$population <- pop_names[1]
  interventions$generation <- 1:generations
  
  if(length(intervention_prob) > 1){
    for(i in 2:length(intervention_prob)){
      
      tmp <- dplyr::as_tibble(matrix(rbinom(generations * loci_length,
                                            size = 1,
                                            prob = intervention_prob[i]),
                                     nrow = generations,
                                     ncol = loci_length, dimnames = list(c(), loci_names)))
      tmp$population <- pop_names[i]
      tmp$generation <- 1:generations
      interventions <- dplyr::bind_rows(interventions, tmp)
    }
  }
  
  return(interventions)
}

myargs <- commandArgs(trailingOnly=TRUE)

sigma_f <- as.numeric(myargs[1])
repetition_id <- myargs[2]
fa <- as.numeric(myargs[3])
m <- as.numeric(myargs[4])
ab_pressure = as.numeric(myargs[5])
study = myargs[6]
out <- run_simulation(antibiotic_pressure = ab_pressure,
                      study = study,
                      generations=200,
                      repetitions=1,
                      sigma_f=sigma_f,
                      fa=fa,
                      m=m)

filename <- paste0("~/results/summarised_state_",
                   sigma_f, "_",
                   fa, "_",
                   m, "_",
                   repetition_id, "_",
                   ab_pressure, "_",
                   study, ".RData")

saveRDS(out, file=filename)
