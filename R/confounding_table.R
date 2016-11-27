#  Okay, I have three types of change:
#  (1) compounding
#  (2) confounding
#  (3) contradicting
#
#  * In compounding delta_C & LU should have the same net sign relative to the 
#  underlying veg data.
#  * In confounding the LU is in the middle of the distribution and plays little 
#  to no role
#  * In contradicting we should see that the LU has a different sign relative
#  to the clim change.

taxa <- unique(vegclim_table$taxon)
clim_var <- unique(vegclim_table$climate)

conf_test <- function(taxa, clim_var, climtable) {

  require(dplyr)
  require(SpatialPack)
  
  coords <- climtable %>% select(cell, x, y) %>% distinct
  
  # Get all values, join to coords so we have the same number of cells across all datasets
  # then make sure we have only distinct sets.
  
  VHCH <- climtable %>% filter(c.ref == "PLSS" & 
                               climate == clim_var & 
                               taxon == taxa &
                               base == "PLSS") %>% 
    select(cell, x, y, clim, data) %>% 
    left_join(coords, .)
  
  VHCH$clim[VHCH$data == 0] <- NA
  
  VHCM <- climtable %>% filter(c.ref == "FIA" & 
                                 climate == clim_var & 
                                 taxon == taxa &
                                 base == "PLSS") %>% 
    select(cell, x, y, clim, data) %>% 
    left_join(coords, .) %>% distinct(cell, .keep_all = TRUE)
  
  VHCM$clim[VHCM$data == 0] <- NA
  
  VMCH <- climtable %>% filter(c.ref == "PLSS" & 
                                 climate == clim_var & 
                                 taxon == taxa &
                                 base == "FIA") %>% 
    select(cell, x, y, clim) %>% 
    left_join(coords, .) %>% distinct(cell, .keep_all = TRUE)
  
  VMCH$clim[VMCH$data == 0] <- NA
  
  all_cell <- coords$cell[coords$cell %in% VHCH$cell &
                           coords$cell %in% VHCM$cell &
                           coords$cell %in% VMCH$cell]
  
  # Get the t estimates, we're going to use the re-calculated df from the
  # modified t-test, which is pairwise (we don't want it pairwise).
  clim_t  <- t.test((VHCH %>% filter(cell %in% all_cell))$clim, 
                    (VHCM %>% filter(cell %in% all_cell))$clim)
  
  lu_t    <- t.test((VHCH %>% filter(cell %in% all_cell))$clim, 
                    (VMCH %>% filter(cell %in% all_cell))$clim)
  
  clim_mt <- modified.ttest((VHCH %>% filter(cell %in% all_cell))$clim, 
                            (VHCM %>% filter(cell %in% all_cell))$clim, 
                            coords = (coords %>% filter(cell %in% all_cell))[,-1])
  
  clim_df <- summary(clim_mt)$dof
  
  lu_mt <- modified.ttest((VHCH %>% filter(cell %in% all_cell))$clim, 
                            (VMCH %>% filter(cell %in% all_cell))$clim, 
                            coords = (coords %>% filter(cell %in% all_cell))[,-1])
  
  lu_df <- summary(lu_mt)$dof
  
  # We want to return the t-test p value & estimate
  data.frame(taxa = taxa, clim = clim_var,
             clim_est = diff(clim_t$estimate), 
             clim_p = 2 * pt(abs(clim_t$statistic), clim_df, lower.tail = FALSE),
             lu_est = diff(lu_t$estimate),   
             lu_p   = 2 * pt(abs(lu_t$statistic),     lu_df, lower.tail = FALSE))
  
}

tests <- do.call(rbind.data.frame,
        lapply(taxa, function(x) {
          do.call(rbind.data.frame,lapply(clim_var, function(y) {conf_test(x, y, climtable = vegclim_table)}))}))

tests_pls <- do.call(rbind.data.frame,
                           lapply(taxa, function(x) {
                             do.call(rbind.data.frame,lapply(clim_var, function(y) {conf_test(x, y, climtable = newveg)}))}))

tests_fia <- do.call(rbind.data.frame,
                           lapply(taxa, function(x) {
                             do.call(rbind.data.frame,lapply(clim_var, function(y) {conf_test(x, y, climtable = topveg)}))}))

tests_all <- do.call(rbind.data.frame,
                     lapply(taxa, function(x) {
                       do.call(rbind.data.frame,lapply(clim_var, function(y) {conf_test(x, y, climtable = allveg)}))}))

checker <- function(x) {
  # x will be the row:

  
  if (p.adjust(as.numeric(x[4]), method = "bonferroni", n = 120) < 0.05 & 
      p.adjust(as.numeric(x[6]), method = "bonferroni", n = 120) < 0.05) {
    
    # both p values are significant:
    if ((as.numeric(x[3]) * as.numeric(x[5])) > 0) {
      # A simple test to see if they're both positive
      return(ifelse(as.numeric(x[3]) > 0, '$\\Box$ ($+_c$, $+_v$)', '$\\Box$ ($-_c$,$-_v$)'))
      
    } else {
      # They are in different directions:
      return(ifelse(as.numeric(x[3]) > 0, '$\\bigcirc$ ($+_c$,$-_v$)', '$\\bigcirc$ ($-_c$,$+_v$)'))
    }
    
  } else {
    # If one or the either is not significant:
    if (as.numeric(x[4]) > (0.05 / 120) & as.numeric(x[6]) > (0.05 / 120)) {
      # If neither element shows significant change
      return('n/a ($._c$, $._v$)')
    } else {
      if (as.numeric(x[4]) > (0.05 / 120)) {
        # If climate change is significant, but the land use change is not significant.
        if (as.numeric(x[3]) > 0) {
          return('n/a ($+_c$, $._v$)')
        } else {
          return('n/a ($-_c$, $._v$)')
        }
      } else {
        if (as.numeric(x[5]) > 0) {
          return('n/a ($._c$, $+_v$)')
        } else {
          return('n/a ($._c$, $-_v$)')
        }
      }
    }
  }
}

get_conf <- function(x) {
  x$result <- apply(x, 1, checker)
  dcast(x, taxa ~ clim, value.var = 'result')
}

conf_table <- get_conf(tests)
conf_pls   <- get_conf(tests_pls)
conf_fia   <- get_conf(tests_fia)
conf_all   <- get_conf(tests_all)