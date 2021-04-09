################################################################################
################################## run mcmc ####################################
################################################################################

set.seed(73762920)

# e coli
ec_mcmc_out <- post_prob(
  file_path = "G:\\My Drive\\Data\\transposon data\\gene_level_data\\E_coli_insdens.csv",
  insdens_thresh = 0.025,
  init_theta = 0.05,
  niter = 22000,
  burn_in = 2000,
  acc_window = 1000,
  ae_sd = 0.25,
  be_sd = 0.30,
  an_sd = 0.04,
  bn_sd = 0.05,
  print_prop_sd = T,
  print_it = T,
  thin = 1)
ec_ppe <- ec_mcmc_out[[1]]

par(mfrow = c(3,2))
lab_mag <- 1.5
lab_annot <- 1.4
par(mar =  c(4, 5, 2, 1) + 0.1)
plot(ec_mcmc_out[[2]]$alpha_e, xlab = "Iteration", ylab = expression(alpha[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(ec_mcmc_out[[2]]$alpha_n, xlab = "Iteration", ylab = expression(alpha[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(ec_mcmc_out[[2]]$beta_e, xlab = "Iteration", ylab = expression(beta[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(ec_mcmc_out[[2]]$beta_n, xlab = "Iteration", ylab = expression(beta[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(ec_mcmc_out[[2]]$theta, xlab = "Iteration", ylab = expression(theta),
     cex.lab = lab_mag, cex.axis = lab_annot)


#S aureus
sa_mcmc_out <- post_prob(file_path = "G:\\My Drive\\Data\\transposon data\\gene_level_data\\S_aureus_insdens.csv",
                         insdens_thresh = 0.025,
                         init_theta = 0.05,
                         niter = 22000,
                         burn_in = 2000,
                         acc_window = 1000,
                         ae_sd = 0.25,
                         be_sd = 0.30,
                         an_sd = 0.04,
                         bn_sd = 0.05,
                         print_prop_sd = T,
                         print_it = T,
                         thin = 1)
sa_ppe <- sa_mcmc_out[[1]]

par(mfrow = c(3,2))
lab_mag <- 1.5
lab_annot <- 1.4
par(mar =  c(4, 5, 2, 1) + 0.1)
plot(sa_mcmc_out[[2]]$alpha_e, xlab = "Iteration", ylab = expression(alpha[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(sa_mcmc_out[[2]]$alpha_n, xlab = "Iteration", ylab = expression(alpha[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(sa_mcmc_out[[2]]$beta_e, xlab = "Iteration", ylab = expression(beta[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(sa_mcmc_out[[2]]$beta_n, xlab = "Iteration", ylab = expression(beta[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(sa_mcmc_out[[2]]$theta, xlab = "Iteration", ylab = expression(theta),
     cex.lab = lab_mag, cex.axis = lab_annot)



# S Typhimurium
stmm_mcmc_out <- post_prob(file_path = "G:\\My Drive\\Data\\transposon data\\gene_level_data\\S_Typhimurium_insdens.csv",
                      insdens_thresh = 0.025,
                      init_theta = 0.05,
                      niter = 22000,
                      burn_in = 2000,
                      acc_window = 1000,
                      ae_sd = 0.25,
                      be_sd = 0.30,
                      an_sd = 0.04,
                      bn_sd = 0.05,
                      print_prop_sd = T,
                      print_it = T,
                      thin = 1)
stmm_ppe <- stmm_mcmc_out[[1]]

par(mfrow = c(3,2))
lab_mag <- 1.5
lab_annot <- 1.4
par(mar =  c(4, 5, 2, 1) + 0.1)
plot(stmm_mcmc_out[[2]]$alpha_e, xlab = "Iteration", ylab = expression(alpha[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(stmm_mcmc_out[[2]]$alpha_n, xlab = "Iteration", ylab = expression(alpha[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(stmm_mcmc_out[[2]]$beta_e, xlab = "Iteration", ylab = expression(beta[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(stmm_mcmc_out[[2]]$beta_n, xlab = "Iteration", ylab = expression(beta[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(stmm_mcmc_out[[2]]$theta, xlab = "Iteration", ylab = expression(theta),
     cex.lab = lab_mag, cex.axis = lab_annot)


################################################################################
################################## plot ROC ####################################
################################################################################

#s.aureus
load("G:\\My Drive\\Data\\transposon data\\S_aureus_ess_status.RData")
sa_ppe_and_status <- merge(sa_ppe, sa_true_ess_status)
pred <- ROCR::prediction(sa_ppe_and_status$post_prob,
                         labels = sa_ppe_and_status$status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf, xlim = c(0, 1), ylim = c(0, 1),
           xlab = (expression(bold("False positive rate"))),
           ylab = (expression(bold("True positive rate"))),
           lty = 2, col = "black", cex.lab =1.2, cex.axis = 1.8,
           lwd = 3, xaxis.cex.axis = 1.5, yaxis.cex.axis = 1.5,
           xaxis.lwd = 2, yaxis.lwd = 3)
auc <- unlist(ROCR::performance(pred,"auc")@y.values)
sa_auc <- round(auc, 3)

# ecoli
load("G:\\My Drive\\Data\\transposon data\\E_coli_ess_status.RData")
ec_ppe_and_status <- merge(ec_ppe, ec_true_ess_status)
pred <- ROCR::prediction(ec_ppe$post_prob, labels = ec_true_ess_status$status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf,avg="vertical", add = T, lty = 1,col = "black",lwd=3)
auc <- unlist(ROCR::performance(pred,"auc")@y.values)
ec_auc <-round(auc, 3)

# S typhimurium
load("G:\\My Drive\\Data\\transposon data\\S_Typhimurium_ess_status.RData")
stmm_ppe_and_status <- merge(stmm_ppe, stmm_true_ess_status)
pred <- ROCR::prediction(stmm_ppe$post_prob, labels = stmm_true_ess_status$status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf,avg="vertical", add = T, lty = 3,col = "black",lwd=3)
auc <- unlist(ROCR::performance(pred,"auc")@y.values)
stmm_auc <-format(round(auc, 3), nsmall = 3)

legend(x = 0.25, y= 0.32,
       legend = c(paste(" S.aureus, AUC = ", sa_auc, sep=""),
                  paste(" E.coli, AUC = ", ec_auc, sep=""),
                  paste(" S.Typhimurium, AUC = ", stmm_auc, sep="")),
       col = c("black", "black", "black"),
       lty = c(2, 1, 3),
       pch = 10,
       pt.cex = 0,
       #text.font = 10,
       cex = 1.6,
       horiz = FALSE,
       bty = "n",
       lwd = 3,
       x.intersp = 0,
       y.intersp = 1.4)

save(list = c("sa_ppe_and_status", "ec_ppe_and_status", "stmm_ppe_and_status"),
     file = paste0("G:\\My Drive\\Data\\transposon data\\gene_level_data\\",
                   "data_and_code_for_insdens_paper\\ppe_and_status_real_data.RData"))


