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
  clim = t.test(subset(climtable, taxon == taxa & climate == clim_var & base == 'PLSS' & c.ref == 'PLSS' & data > 0)$clim,
         subset(climtable, taxon == taxa & climate == clim_var & base == 'PLSS' & c.ref == 'FIA' & data > 0)$clim)

  lu   = t.test(subset(climtable, taxon == taxa & climate == clim_var & base == 'PLSS' & c.ref == 'PLSS' & data > 0)$clim,
         subset(climtable, taxon == taxa & climate == clim_var & base == 'FIA' & c.ref == 'PLSS' & data > 0)$clim)

  # We want to return the t-test p value & estimate
  data.frame(taxa = taxa, clim = clim_var,
             clim_est = diff(clim$estimate), clim_p = clim$p.value,
             lu_est   = diff(lu$estimate),   lu_p   = lu$p.value)
  
}

tests <- do.call(rbind.data.frame,
        lapply(taxa, function(x) {
          do.call(rbind.data.frame,lapply(clim_var, function(y) {conf_test(x, y, climtable = vegclim_table)}))}))

tests_corrected <- do.call(rbind.data.frame,
                           lapply(taxa, function(x) {
                             do.call(rbind.data.frame,lapply(clim_var, function(y) {conf_test(x, y, climtable = newveg)}))}))

checker <- function(x) {
  # x will be the row:
  # compounding:
  # both ps are sig & both ests are in the same direction:
  if (p.adjust(as.numeric(x[4]), method = "bonferroni", n = 120) < 0.05 & 
      p.adjust(as.numeric(x[6]), method = "bonferroni", n = 120) < 0.05) {
    if ((as.numeric(x[3]) * as.numeric(x[5])) > 0) {
      return('$\\Box$')
    } else {
      return('$\\bigcirc$')
    }
  } else {
    if (as.numeric(x[4]) > (0.05 / 120)) {
      return('n/a')
    } else {
      return('-')
    }
  }
}

tests$result <- apply(tests, 1, checker)
conf_table <- dcast(tests, taxa ~ clim, value.var = 'result')

tests_corrected$result <- apply(tests_corrected, 1, checker)
conf_table_corr <- dcast(tests_corrected, taxa ~ clim, value.var = 'result')


conf_table <- conf_table[match(taxa, conf_table$taxa),]

mean_y <- function(taxa, climtable) {
  mean(subset(climtable, taxon == taxa & base == 'PLSS' & data > 0)$y)
}

mean_x <- function(taxa, climtable) {
  mean(subset(climtable, taxon == taxa & base == 'PLSS' & data > 0)$x)
}

conf_table$mean_y <- sapply(taxa, mean_y, climtable = vegclim_table)
conf_table$mean_x <- sapply(taxa, mean_x, climtable = vegclim_table)

conf_table <- conf_table[,1:5]

conf_cor <- conf_table[,2:5]
conf_cor[conf_cor == '$\\Box$'] <- 1
conf_cor[conf_cor == '-'] <- 0
conf_cor[conf_cor == '$\\bigcirc$'] <- -1

rownames(conf_cor) <- conf_table[,1]
cluster <- hclust(dist(conf_cor))

# Sort the table by similarity:
conf_table <- conf_table[cluster$order, ]
rownames(conf_table) <- NULL
