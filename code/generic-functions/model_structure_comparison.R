# function to compare models with different intercept structures

mod_structure_comp <- function(simp_mods, ran_mods, fix_mods) {
  
  # empty lists
  simp_sum <- vector(mode = "list", length = length(simp_mods))
  ran_sum <- vector(mode = "list", length = length(ran_mods))
  fix_sum <- vector(mode = "list", length = length(fix_mods))
  
  # models with one intercept
  for(i in 1:length(simp_mods)){
    
    simp_sum[[i]] <- tibble(coef(simp_mods[[i]])) # extract model coefficients
    colnames(simp_sum[[i]]) <- names(simp_mods)[i] # rename column to model name
    simp_sum[[i]] <- simp_sum[[i]] %>%
      mutate(coefficients = names(coef(simp_mods[[i]])) %>%
               str_replace_all("PercCovered", "PAC") %>%
               str_replace_all("RecentTreatment", "Treated") %>%
               str_replace_all("PAC:Treated", "Interaction")) # add coefficient names
    
  }
  
  # mixed-effects models
  for(i in 1:length(ran_mods)){
    
    ran_sum[[i]] <- tibble(fixef(ran_mods[[i]])$cond) # extract model coefficients
    colnames(ran_sum[[i]]) <- names(ran_mods)[i] # rename column to model name
    ran_sum[[i]] <- ran_sum[[i]] %>%
      mutate(coefficients = names(fixef(ran_mods[[i]])$cond) %>%
               str_replace_all("PercCovered", "PAC") %>%
               str_replace_all("RecentTreatment", "Treated") %>%
               str_replace_all("PAC:Treated", "Interaction")) # add coefficient names
    
  }
  
  # panel data models
  for(i in 1:length(fix_mods)){
    
    fix_sum[[i]] <- tibble(coef(fix_mods[[i]])) # extract model coefficients
    colnames(fix_sum[[i]]) <- names(fix_mods)[i] # rename column to model name
    fix_sum[[i]] <- fix_sum[[i]] %>%
      mutate(coefficients = names(coef(fix_mods[[i]])) %>%
               str_replace_all("PercCovered", "PAC") %>%
               str_replace_all("RecentTreatment", "Treated") %>%
               str_replace_all("PAC:Treated", "Interaction")) # add coefficient names
    
  }

  # combine dataframes
  mod_comp <- reduce(simp_sum, full_join) %>%
    full_join(reduce(ran_sum, full_join)) %>%
    full_join(reduce(fix_sum, full_join)) %>%
    relocate(coefficients)
  
  # output
  return(mod_comp)
  
}
