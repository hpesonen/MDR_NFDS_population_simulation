# library(tidyverse)
# library(DirichletReg)
library(dplyr)
# library(tidyverse)
library(readr)
library(DirichletReg)

# run_simulation <- function(intervention_prob = 1,
#                            repetitions = 1,
#                            antibiotic_pressure = 0,
#                            generations = 120,
#                            fa = 0,
#                            m = 0.01,
#                            sigma_f = 0.1,
#                            v = 0,
#                            within_str_trans = 0.1 / 10000,
#                            between_str_trans = 0.1 / 10000,
#                            by_mutation = 0.001 / 10000,
#                            alpha = c(1,1),
#                            beta = c(10,10),
#                            transmissions = matrix(c(20,2,1), nrow=1, ncol=3),
#                            seed1 = -99,
#                            seed2 = -99){

run_simulation <- function(intervention_prob = 1,
                           repetitions = 1,
                           antibiotic_pressure = 0.006,
			   study = "NORM",
                           generations = 120,
                           fa = 0,
                           m = 0.01,
                           sigma_f = 0.1,
                           v = 0,
                           within_str_trans = 0.1 / 10000,
                           between_str_trans = 0.1 / 10000,
                           by_mutation = 0.001 / 10000,
                           alpha = c(1,1),
                           beta = c(10,10),
                           transmissions = matrix(c(20,2,1), nrow=1, ncol=3),
                           seed1 = -99,
                           seed2 = -99){


  if(seed1 != -99) set.seed(seed1)

  print(antibiotic_pressure)
  # strain_amr_data <- parse_strain_amr_file("~/Documents/research/rprojects/ONEHEARTBEAT/data/mmc2_full.csv", "ST", "CTX-M", 30)
  # print("here)
  
  if(study == "NORM"){
    strain_amr_data <- study_parse_strain_amr_file("~/mmc2_full.csv", "ST", "CTX-M", 20, "NORM")
  }
  if(study == "BSAC_NORM"){
    strain_amr_data_NORM <- study_parse_strain_amr_file("~/mmc2_full.csv", "ST", "CTX-M", 20, "NORM")
    strain_amr_data <- study_parse_strain_amr_file("~/bsac_ctxm_ed.csv", "ST", "CTX-M", 20, "BSAC")
  }
  if(study == "BSAC"){
    # strain_amr_data_NORM <- study_parse_strain_amr_file("~/mmc2_full.csv", "ST", "CTX-M", 20, "NORM")
    strain_amr_data <- study_parse_strain_amr_file("~/bsac_ctxm_ed.csv", "ST", "CTX-M", 20, "BSAC")
  }

  # strain_amr_data <- study_parse_strain_amr_file("~/mmc2_full.csv", "ST", "CTX-M", 20)
  #strain_amr_data <- parse_strain_amr_file("~/mmc2_full.csv", "ST", "CTX-M", 20)
  # return(strain_amr_data)
  genotype <- names(strain_amr_data$strains)

  # TODO separate data-parsing with relevant error-checks into separate function

  pop_names <- factor("NORM")
  pop_size <- 50000
  kappa <- pop_size

  migration_pressure <- 0

  populations <- length(pop_size)
  # loci_names <- ifelse(!is.null(IFloci)*rep(TRUE, length(IFloci$locus_name)), IFloci$locus_name, NA)
  # loci_names = c("CTX-M", "AB1", "AB2", "AB3", "AB4", "AB5", "AB6", "AB7")
  # loci_names = c("AB1", "AB2", "AB3", "AB4", "AB5", "AB6", "AB7", "AB8", "AB9", "AB10", "AB11" )
  loci_names <- c("CTXM", paste0("L", 1:20))
  IFloci <- loci_names
  # print(loci_names)
  loci_length <- length(loci_names)
  genotype_names <- genotype
  # genotype_names <- factor(genotype$genotype_names, level=genotype$genotype_names)
  genotype_length <- length(genotype_names)

  # Notice these!!!
  e <- strain_amr_data$strains
  # e[length(e)] <- 0.0041
  e <- e/sum(e)
  # print(e)
  # return()
  # e <- c(0.13521127, 0.10704225, 0.09577465, 0.08732394, 0.04225352, 0.03661972, 0.02535211, 0.02253521, 0.02253521, 0.01971831, 0.040563380)
  # e <- e/sum(e)

  # How large is the proportionality of AMR initially
  # Assumed the same for all populations
  # e_pre <- diag(c(rep(1, genotype_length)))
  # e_pre <- diag(e)
  # print(e_pre)
  e_pre <- matrix(0, nrow=genotype_length, ncol=loci_length)
  rownames(e_pre) <- genotype_names
  colnames(e_pre) <- loci_names
  
  # alphas <- diag(rep(100000, loci_length)) + 0.0001
  # alphas <-diag(rep(100000, loci_length-1)) + 0.0001
  # alphas <- cbind(0.001, alphas)
  # print(alphas)
    # alphas <- matrix(1000*runif(loci_length * genotype_length), nrow=genotype_length, ncol=loci_length)
  
  alphas <- 10000*rt(loci_length * genotype_length, df=2) - 100
  alphas <- pmax(alphas, 0.001)
  alphas <- matrix(alphas, nrow=genotype_length, ncol=loci_length)
  # alphas[,2] <- c(0, 100000, rep(0, genotype_length - 2)) + 0.001
  # print(intervention_prob)
  # print(alphas)
  e_pre_1 <- e_pre
  e_pre_2 <- e_pre
  for(i in 1:genotype_length){
    e_pre_1[i,] <- DirichletReg::rdirichlet(1, c(alphas[i,]))
    e_pre_2[i,] <- e_pre_1[i,]
    # e_pre[i,] <- rbeta(loci_length, alpha[i], beta[i])  
  }
  e_pre_1[,1] <- strain_amr_data$amr[,2] + 0.01
  e_pre_2[,1] <- (strain_amr_data$amr[,2] + 0.01) * 0.1
  # e_pre[,1] <- (strain_amr_data$amr[,2] + 0.01) * 0.1

  # e_pre[,1] <- 0.001
  # print(e_pre)
  # return()
  # rownames(e_pre) <- genotype_names
  # colnames(e_pre) <- loci_names
  # print(round(e_pre,2))
  # return()
  # matrix(rbeta(genotype_length * loci_length, alpha[i], beta[i]), nrow=genotype_length, ncol=loci_length)  
    
  # print(strain_amr_data$amr)
  # print(as_tibble(strain_amr_data$amr) %>% dplyr::select("CTX-M"))
  # e_pre <- as_tibble(strain_amr_data$amr) %>% dplyr::select("CTX-M")
  # rownames(e_pre) <- genotype
  # return()
  
  # e_pre <- as_tibble(strain_amr_data$amr) %>% dplyr::select("CTX-M")
  # e_pre <- pmax(strain_amr_data$amr[,2], 0.1)

  # print(e_pre)
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

  # state <- init_population(pop_names, pop_size, genotype_names, loci_names, kappa, e, e_pre)
  # state_stable <- state
  print("Initialised")
  # return()
  # print(table(state$genotype))
  e_init <- e
  e_init[genotype_names == "69" | genotype_names == "131_A" | genotype_names == "131_B" | genotype_names == "131_C1" | genotype_names == "131_C2"] <- 0.1 * e_init[genotype_names == "69" | genotype_names == "131_A" | genotype_names == "131_B" | genotype_names == "131_C1" | genotype_names == "131_C2"]
  # e_init[2] <- 0.1*e_init[2]
  e_init <- e_init / sum(e_init)
  
  genotype_index = c()
  for( i in genotype){
    genotype_index <- c(genotype_index, sample(which(state$genotype == i), size=round(e_init[i] * pop_size), replace = TRUE))
  }
  
  state <- state[genotype_index,]
  
  # print("Resampled")
  # print(table(state$genotype))
  # e_pre <- sum(strain_amr_data$amr[,2] * strain_amr_data$strains)
  # compared <- strain_amr_data$strains %*% e_pre
  # e_pre <- colSums(state %>% dplyr::select(IFloci)) / pop_size

  # e_pre <- e %*% e_pre_1
  
  state_init <- state
  # genotype_names <- state %>% dplyr::select(genotype) %>% unique() %>% dplyr::pull(genotype)
  total_rows <- populations * genotype_length * generations * repetitions
  total_rows_pop <- populations * generations * repetitions
  print("Summarising")
  summaries_list <- summarize_population(state, loci_names, pop_names, genotype_names, igeneration = 1, irepetition = 1, ab_pressure = antibiotic_pressure[1])
  
  summarized_state_genotype_init <- summaries_list[[1]]
  # summarized_state_population_init <- summaries_list[[2]]
  print("Start simulation")
  # return()
  if(seed2 != -99) set.seed(seed2)
  for(irepetition in 1:repetitions){
    print(paste("Repetition :", irepetition))
    for(ab_pressure in antibiotic_pressure){
      # print("H")
      interventions <- generate_interventions(generations,
                                              loci_length,
                                              intervention_prob,
                                              loci_names,
                                              pop_names)
  
      for(igeneration in 1:generations){
        if(igeneration %% 10 == 0) print(round(igeneration / generations, 2))
        
        # print(igeneration)
        if(igeneration == 1){
  
          # print(paste("AB resistance :", round(antibiotic_pressure[irepetition], 2)))
          print(paste("AB resistance :", round(ab_pressure, 2)))
          state <- state_init
          if(irepetition == 1 & ab_pressure == antibiotic_pressure[1]){
            summarized_state_genotype <- summarized_state_genotype_init
            # # summarized_state_population <- summarized_state_population_init
  
          } else {
            summarized_state_genotype_init$repetition <- irepetition
            # # summarized_state_population_init$repetition <- irepetition
            summarized_state_genotype_init$antibiotic_pressure <- ab_pressure
          # #   summarized_state_population_init$antibiotic_pressure <- ab_pressure

            summarized_state_genotype <- dplyr::bind_rows(summarized_state_genotype, summarized_state_genotype_init)
            # # # summarized_state_population <- dplyr::bind_rows(summarized_state_population, summarized_state_population_init)
          }
        } else {
  
          state <- evolve_generation(state, dplyr::filter(interventions, generation == igeneration),
                                     within_str_trans, between_str_trans,
                                     by_mutation, e_pre, kappa, m, v, sigma_f,
                                     loci_names, transmissions, igeneration, pop_names, genotype_names,
                                     migration_pressure, ab_pressure, state_stable, fa)
          # print("K2")
          summaries_list <- summarize_population(state, loci_names, pop_names, genotype_names, igeneration, irepetition, ab_pressure)
          # if(igeneration %% 10 == 0)
          summarized_state_genotype <- dplyr::bind_rows(summarized_state_genotype, summaries_list[[1]])
          # # # summarized_state_population <- dplyr::bind_rows(summarized_state_population, summaries_list[[2]])
          # print("K3")
        }
      }
    }
  }
  return(list(summarized_state_genotype=summarized_state_genotype))
  # return(list(summarized_state_genotype=summarized_state_genotype,
  #             summarized_state_population=summarized_state_population))
}



###########

# study_parse_strain_amr_file("~/mmc2_full.csv", "ST", "CTX-M", 20, "NORM")
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
    # strain_data <- read.csv(filename)
    # strain_data$`CTX-M` <- (strain_data$blaCTX.M_any == "yes") * 1
    # strain_data <- strain_data %>% dplyr::filter(year==2015 | year==2016 | year==2017)
    # strain_data$ST131_clade[is.na(strain_data$ST131_clade)] <- ""
    # strain_data$subglades <- paste0(strain_data$ST, strain_data$ST131_clade)

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
  pop_count             <- state %>% dplyr::group_by(population)  %>% dplyr::tally() %>% dplyr::ungroup()
  pop_and_genotype_count   <- state %>% dplyr::group_by(population, genotype)  %>% dplyr::tally() %>% dplyr::ungroup()

    pop_ab_prop           <- state %>% dplyr::group_by(population) %>% dplyr::summarize_at(IFloci, mean) %>% dplyr::ungroup()
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

  # summary_table = list(summary_table_genotype, summary_table_population)
  summary_table = list(summary_table_genotype)
  return(summary_table)

}


parse_strain_amr_file_ <- function(filename, strain_column_name, amr_column_name, n_top, i_year){
  # print(strain_column_name)
  strain_data <- readr::read_delim(filename, delim=";")
  strain_data_all <- strain_data
  
  strain_data <- strain_data %>% dplyr::filter(year==2017)
  
  # print("this 1?")
  # print(strain_data$ST)
  # print(strain_data$CC131_clades)
  strain_data$CC131_clades[is.na(strain_data$CC131_clades)] <- ""
  strain_data$subglades <- paste0(strain_data$ST, strain_data$CC131_clades)
  # print(strain_data$subglades)
  # print("this 2?")
  strain_data$subglades[strain_data$subglades == "131A"] <- "131_A"
  strain_data$subglades[strain_data$subglades == "131B"] <- "131_B"
  strain_data$subglades[strain_data$subglades == "131C1"] <- "131_C"
  strain_data$subglades[strain_data$subglades == "131C2"] <- "131_C"
  strain_column_name <- "subglades"
  # print(strain_column_name)
  
  strain_table <- table(strain_data[strain_column_name])
  
  # print(strain_table)
  
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
  years <- 2002:2017
  for(iyear in years){
    strain_data <- strain_data_all %>% dplyr::filter(year==iyear)
    # print(year)
    # print(strain_data_all)
    # print(strain_data)
    # print(dim(strain_data_all))
    # print(dim(strain_data))
    strain_data$CC131_clades[is.na(strain_data$CC131_clades)] <- ""
    strain_data$subglades <- paste0(strain_data$ST, strain_data$CC131_clades)
    # print(strain_data$subglades)
    # print("this 2?")
    strain_data$subglades[strain_data$subglades == "131A"] <- "131_A"
    strain_data$subglades[strain_data$subglades == "131B"] <- "131_B"
    strain_data$subglades[strain_data$subglades == "131C1"] <- "131_C"
    strain_data$subglades[strain_data$subglades == "131C2"] <- "131_C"
    strain_column_name <- "subglades"

    strain_amr <- strain_data %>%
      dplyr::select(amr_column_name) %>%
      mutate_all(~ replace(., . == "CTX-M", 1)) %>%
      mutate_all(~ replace(., is.na(.), 0)) %>%
      tibble::deframe()

    strains <- strain_data[,strain_column_name] %>% tibble::deframe()
    strains <- c(strains, n_top_strain_names, n_top_strain_names)
    strain_amr <- c(strain_amr, rep("0", length(n_top_strain_names)), rep("1", length(n_top_strain_names)))
    amr_table <- table(strains, strain_amr)
    amr_table_top <- amr_table[n_top_strain_names,]

    amr_table_other <- amr_table[(row.names(amr_table) %in% other_strain_names),]

    amr_table_other_amr <- colSums(amr_table_other)

    amr_table_top <- rbind(amr_table_top, amr_table_other_amr)
    
    rownames(amr_table_top)[n_top+1] <- "STX"
    # amr_table_top <- amr_table_top / rowSums(amr_table_top)
    if(iyear == years[1]){
      average_amr <- amr_table_top
    }else{
      average_amr <- average_amr + amr_table_top
    }
  }
  average_amr  <- average_amr / pmax(rowSums(average_amr), 1)

  return(list(strains=strain_table_freq, amr=average_amr / length(years)))
  # return(list(strains=strain_table_freq, amr=average_amr))
  
}


parse_strain_amr_file <- function(filename, strain_column_name, amr_column_name, n_top, i_year){
  # print(strain_column_name)
  strain_data <- readr::read_delim(filename, delim=";")
  strain_data <- strain_data %>% dplyr::filter(year==2017)
  
  # print("this 1?")
  # print(strain_data$ST)
  # print(strain_data$CC131_clades)
  strain_data$CC131_clades[is.na(strain_data$CC131_clades)] <- "" 
  strain_data$subglades <- paste0(strain_data$ST, strain_data$CC131_clades)
  # print(strain_data$subglades)
  # print("this 2?")
  strain_data$subglades[strain_data$subglades == "131A"] <- "131_A"
  strain_data$subglades[strain_data$subglades == "131B"] <- "131_B"
  strain_data$subglades[strain_data$subglades == "131C1"] <- "131_C"
  strain_data$subglades[strain_data$subglades == "131C2"] <- "131_C"
  strain_column_name <- "subglades"
  # print(strain_column_name)
  
  strain_table <- table(strain_data[strain_column_name])
  
  # print(strain_table)
  
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

  
  # strain_data_2002 <- readr::read_delim(filename, delim=";")

  # strain_data_2002 <- left_join(strain_data_2017_strains, strain_data_2002_incomplete, by="ST")

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


###########

parse_strain_amr_file_old <- function(filename, strain_column_name, amr_column_name, n_top, i_year){

  strain_data <- readr::read_delim(filename, delim=";")
  strain_data <- strain_data %>% dplyr::filter(year==2017)
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

  
  # strain_data_2002 <- readr::read_delim(filename, delim=";")

  # strain_data_2002 <- left_join(strain_data_2017_strains, strain_data_2002_incomplete, by="ST")

  strain_amr <- strain_data %>% 
    dplyr::select(amr_column_name) %>% 
    mutate_all(~ replace(., . == "CTX-M", 1)) %>% 
    mutate_all(~ replace(., is.na(.), 0)) %>%
    tibble::deframe()

  strains <- strain_data[,strain_column_name] %>% tibble::deframe()

  amr_table <- table(strains, strain_amr)

  amr_table_top <- amr_table[n_top_strain_names,]

  amr_table_other <- amr_table[!(row.names(amr_table) %in% other_strain_names),]

  amr_table_other_amr <- colSums(amr_table_other)
  amr_table_top <- rbind(amr_table_top, amr_table_other_amr)
  
  rownames(amr_table_top)[n_top+1] <- "STX"
  amr_table_top <- amr_table_top / rowSums(amr_table_top)
  
  return(list(strains=strain_table_freq, amr=amr_table_top))
}



###########


evolve_generation <- function(state,
                              interventions,
                              within_str_trans,
                              between_str_trans,
                              by_mutation,
                              e_pre,
                              kappa,
                              m,
                              v,
                              sigma_f,
                              IFloci,
                              transmissions,
                              igeneration,
                              pop_names,
                              genotypes,
                              migration_pressure,
                              antibiotic_pressure,
                              state_stable, 
                              fa) {


  N_t <- nrow(state)

  # Antibiotic intervention
  # It's assumed that 10% of the population is exposed to the antibiotic
  pop_counts <- table(state$population)
  # interventions_long <- interventions %>% dplyr::slice(rep(1:dplyr::n(), times = as.numeric(pop_counts), each=TRUE)) %>% dplyr::select(IFloci)
  # potential_destroy <-  rowSums((state %>% dplyr::select(IFloci) - interventions_long) < 0)
  interventions_long <- interventions %>% dplyr::slice(rep(1:dplyr::n(), times = as.numeric(pop_counts), each=TRUE)) %>% dplyr::select("CTXM")
  # print(table(interventions_long))
  # potential_destroy <-  rowSums((state %>% dplyr::select("CTXM") - interventions_long) < 0)
  potential_destroy <-  (state %>% dplyr::select("CTXM") == 0)
  destroy <- rep(FALSE, N_t)
  # testme <- rbinom(as.numeric(pop_counts[j]), potential_destroy, antiobiotic_pressure[tmp]) > 0
  for(j in pop_names){
    # destroy[state$population == j] <- rbinom(as.numeric(pop_counts[j]), potential_destroy[state$population == j], antibiotic_pressure) > 0
    destroy[state$population == j] <- (rbinom(as.numeric(pop_counts[j]), 1, antibiotic_pressure) > 0) & potential_destroy[state$population == j]
  }
  # print(table(destroy))
  # print(table(state$genotype, potential_destroy))
  state <- state %>% dplyr::filter(!destroy)

  # state <- rbind(state, state_stable[sample(kappa, m * kappa, replace=TRUE), ])

  # NFDS:
  new_state_index <- c()
  N_t_pop <- 0
  pop_counts <- table(state$population)

  for(j in pop_names){
    N_t <- pop_counts[j]

    # N_t <- state %>% dplyr::filter(population == j) %>% nrow()
    f_t <- colSums(state %>% dplyr::filter(population == j) %>% dplyr::select(IFloci)) / N_t
    diff_vector <- e_pre - f_t
    # print(diff_vector)
    # pi_it <- rowSums(as.matrix((state %>% dplyr::filter(population == j) %>%
    #                    dplyr::select(IFloci))) %*% matrix(rep(diff_vector, length(genotypes)), length(genotypes), length(IFloci), byrow=TRUE))
    present <- state %>% dplyr::filter(population == j) %>% dplyr::select(IFloci)
    pi_it <- rep(0, N_t)
    for(k in 1:length(IFloci)){
    # for(k in 1:1){
      pi_it[present[,k] == 1] = pi_it[present[,k] == 1] + diff_vector[k]
    }
    r_i <- ((state$genotype == "69") | (state$genotype == "131_A") | (state$genotype == "131_B") | (state$genotype == "131_C")) * fa
    # r_i = ((state$genotype == "69") | (state$genotype == "131")) * fa
    lambda <- (kappa / N_t) * (1 + r_i) * (1 - m) * (1 - v) * (1 + sigma_f) ^ pi_it
    # lambda <- (kappa / N_t) * (1 - m) * (1 - v) * (1 + sigma_f) ^ pi_it
    X_it <- rpois(N_t, lambda)
    index = which(state$population == j)

    new_state_index <- c(new_state_index, rep(index, X_it))
  }

  state <- state[new_state_index,]

  N_stable <- nrow(state_stable)
  N_state <- nrow(state)
  # if(N_stable - N_state > 0){
  #    state <- rbind(state, state_stable[sample(N_stable, N_stable - N_state, replace=TRUE), ])
  # }
  state <- rbind(state, state_stable[sample(kappa, m * kappa, replace=TRUE), ])
  # state <- state[-sample(kappa, m * kappa, replace=TRUE),]
  # state <- rbind(state, state_stable[sample(kappa, m * kappa, replace=TRUE), ])
  
  return(state)

}


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
  # state <- dplyr::bind_cols(dplyr::tibble(population = factor(rep(pop_names, times = pop_size), levels = pop_names),
  #                                           genotype = factor(X_t, levels = genotype_names)), AMR_it)
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


# print("start")
myargs <- commandArgs(trailingOnly=TRUE)
# print("myargs")
sigma_f <- as.numeric(myargs[1])
repetition <- myargs[2]
fa <- as.numeric(myargs[3])

m <- as.numeric(myargs[4])

ab_pressure = as.numeric(myargs[5])

study = myargs[6]

# sigma_f <- 0.1
# repetition <- 1

# out <- run_simulation(antibiotic_pressure = c(seq(0.0, 0.02, 0.002)), generations=500, repetitions=3)
# out <- run_simulation(antibiotic_pressure = c(0.001, 0.005, 0.01, 0.03), generations=500, repetitions=50)
# print("start")
# out <- run_simulation(antibiotic_pressure = c(0.0004, 0.001, 0.002, 0.006, 0.03), generations=500, repetitions=1, sigma_f=sigma_f, fa=fa, m=m)
out <- run_simulation(antibiotic_pressure = ab_pressure, study = study, generations=500, repetitions=1, sigma_f=sigma_f, fa=fa, m=m)
# print("end")

saveRDS(out, file=paste0("~/lancet_BSAC_NORM_3op/summarised_state_", sigma_f, "_", fa, "_", m, "_", repetition, "_", ab_pressure, "_", study, ".RData"))
