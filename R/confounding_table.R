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

conf_test <- function(taxa, clim_var, climtable, overlap = FALSE) {

  require(dplyr)
  require(SpatialPack)
  
  coords <- climtable %>% select(cell, x, y) %>% distinct
  
  big_join <- coords %>% full_join(climtable, by = c("cell", "x", "y"))
  
  # Get all values, join to coords so we have the same number of cells across all datasets
  # then make sure we have only distinct sets.
  
  VHCH <- climtable %>% 
    filter(c.ref == "PLSS" & 
             climate == clim_var & 
             taxon == taxa &
             base == "PLSS" & data > 0) %>% 
    select(cell, x, y, clim, data) %>% right_join(coords, by = c("cell", "x", "y"))
  
  VHCH$clim[VHCH$data == 0] <- NA
  
  VHCM <- climtable %>% 
    filter(c.ref == "FIA" & 
             climate == clim_var & 
             taxon == taxa &
             base == "PLSS" & data > 0) %>% 
    select(cell, x, y, clim, data) %>% right_join(coords, by = c("cell", "x", "y"))
  
  VHCM$clim[VHCM$data == 0] <- NA
  
  VMCH <- climtable %>% 
    filter(c.ref == "PLSS" & 
             climate == clim_var & 
             taxon == taxa &
             base == "FIA" & data > 0) %>% 
    select(cell, x, y, clim, data) %>% right_join(coords, by = c("cell", "x", "y"))
  
  VMCH$clim[VMCH$data == 0] <- NA
  
  all_cell <- coords$cell[coords$cell %in% VHCH$cell &
                            coords$cell %in% VHCM$cell &
                            coords$cell %in% VMCH$cell]
  
  # Get the t estimates, we're going to use the re-calculated df from the
  # modified t-test, which is pairwise (we don't want it pairwise).
  if (overlap == FALSE) {
    clim_t  <- t.test(VHCH$clim, 
                      VHCM$clim)
    lu_t    <- t.test(VHCH$clim, 
                      VMCH$clim)
    
    clim_mt <- modified.ttest(VHCH$clim, 
                              VHCM$clim, 
                              coords = coords[,-1])
    
    lu_mt <- modified.ttest(VHCH$clim, 
                            VMCH$clim, 
                            coords = coords[,-1])
    
    
  } else {
    
    clim_t  <- t.test((VHCH %>% filter(cell %in% all_cell))$clim, 
                      (VHCM %>% filter(cell %in% all_cell))$clim)
    
    lu_t    <- t.test((VHCH %>% filter(cell %in% all_cell))$clim, 
                      (VMCH %>% filter(cell %in% all_cell))$clim)

    clim_mt <- modified.ttest((VHCH %>% filter(cell %in% all_cell))$clim, 
                              (VHCM %>% filter(cell %in% all_cell))$clim, 
                              coords = (coords %>% filter(cell %in% all_cell))[,-1])
  
    lu_mt <- modified.ttest((VHCH %>% filter(cell %in% all_cell))$clim, 
                            (VMCH %>% filter(cell %in% all_cell))$clim, 
                            coords = (coords %>% filter(cell %in% all_cell))[,-1])
    
  }
  
  clim_df <- summary(clim_mt)$dof
  lu_df <- summary(lu_mt)$dof
  
  # We want to return the t-test p value & estimate
  data.frame(taxa = taxa, clim = clim_var,
             clim_est = diff(clim_t$estimate), 
             clim_p = 2 * pt(abs(clim_t$statistic), clim_df, lower.tail = FALSE),
             lu_est = diff(lu_t$estimate),
             lu_p   = 2 * pt(abs(lu_t$statistic), lu_df, lower.tail = FALSE),
             clim_unp = clim_t$p.value,
             lu_unp   = lu_t$p.value)
  
}

tests <- do.call(rbind.data.frame,
                 lapply(unique(vegclim_table$taxon), function(x) {
                   do.call(rbind.data.frame,lapply(unique(vegclim_table$climate), 
                                                   function(y) {conf_test(x, y, climtable = vegclim_table)}))}))

tests_pls <- do.call(rbind.data.frame,
                     lapply(taxa, function(x) {
                       do.call(rbind.data.frame,
                               lapply(clim_var, function(y) {conf_test(x, y, climtable = newveg)}))}))

tests_fia <- do.call(rbind.data.frame,
                     lapply(taxa, function(x) {
                       do.call(rbind.data.frame,
                               lapply(clim_var, function(y) {conf_test(x, y, climtable = topveg)}))}))

tests_all <- do.call(rbind.data.frame,
                     lapply(taxa, function(x) {
                       do.call(rbind.data.frame,
                               lapply(clim_var, function(y) {conf_test(x, y, climtable = allveg)}))}))

checker <- function(x) {
  # x will be the row:
  
  csig <- p.adjust(as.numeric(x[4]), method = "bonferroni", n = 120) < 0.05
  
  if (csig == TRUE & as.numeric(x[3]) > 0) { clim <- "$+_c$" }
  if (csig == TRUE & as.numeric(x[3]) < 0) { clim <- "$-_c$" }
  if (csig == FALSE) { clim <- "$._c$" }
  
  veg  <- ifelse(p.adjust(as.numeric(x[6]), method = "bonferroni", n = 120) < 0.05, 
                 ifelse(as.numeric(x[3]) > 0, "$+_v$", "$-_v$"),
                 "$._v$")
  
  if (clim == "$._c$") {
    # There's not change in climate.
    symb <- "$-$"
  } else {
    if (substr(clim, 1, 3) == substr(veg, 1, 3)) {
      # They're equivalent and significant.
      symb <- "$\\Box$"
    } else {
      if (veg == "$._v$") {
        # They're not equivalent, but the veg isn't significant
        symb <- "$-$"
      } else {
        # They're not equivalent, and veg is significant
        symb <- "$\\bigcirc$"
      }
    }
  }

  return(paste0(symb, " (", clim, ",", veg, ")"))
}

get_conf <- function(x) {
  result <- apply(x, 1, checker)
  x$result <- as.character(result)
  out <- reshape2::dcast(x, taxa ~ clim, value.var = 'result')
}

conf_table <- get_conf(tests)
conf_pls   <- get_conf(tests_pls); saveRDS(conf_pls, file = 'data/output/conf_pls.RDS')
conf_fia   <- get_conf(tests_fia); saveRDS(conf_fia, file = 'data/output/conf_fia.RDS')
conf_all   <- get_conf(tests_all); saveRDS(conf_all, file = 'data/output/conf_all.RDS')
