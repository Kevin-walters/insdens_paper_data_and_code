set.seed(9828666)

####################################### HH #####################################
hh_out <- post_prob(file_path = "G:\\My Drive\\Data\\transposon data\\gene_level_data\\sim_data\\sim_insdens_HH.csv",
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
hh_roc <- hh_out[[1]]
hh_post_par_vals <- hh_out[[2]]
par(mfrow = c(3,2))
lab_mag <- 1.5
lab_annot <- 1.4
par(mar =  c(4, 5, 2, 1) + 0.1)
plot(hh_out[[2]]$alpha_e, xlab = "Iteration", ylab = expression(alpha[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(hh_out[[2]]$alpha_n, xlab = "Iteration", ylab = expression(alpha[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(hh_out[[2]]$beta_e, xlab = "Iteration", ylab = expression(beta[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(hh_out[[2]]$beta_n, xlab = "Iteration", ylab = expression(beta[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(hh_out[[2]]$theta, xlab = "Iteration", ylab = expression(theta),
     cex.lab = lab_mag, cex.axis = lab_annot)


####################################### HL #####################################

hl_out <- post_prob(file_path = "G:\\My Drive\\Data\\transposon data\\gene_level_data\\sim_data\\sim_insdens_HL.csv",
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
hl_roc <- hl_out[[1]]
hl_post_par_vals <- hl_out[[2]]
par(mfrow = c(3,2))
lab_mag <- 1.5
lab_annot <- 1.4
par(mar =  c(4, 5, 2, 1) + 0.1)
plot(hl_out[[2]]$alpha_e, xlab = "Iteration", ylab = expression(alpha[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(hl_out[[2]]$alpha_n, xlab = "Iteration", ylab = expression(alpha[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(hl_out[[2]]$beta_e, xlab = "Iteration", ylab = expression(beta[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(hl_out[[2]]$beta_n, xlab = "Iteration", ylab = expression(beta[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(hl_out[[2]]$theta, xlab = "Iteration", ylab = expression(theta),
     cex.lab = lab_mag, cex.axis = lab_annot)



####################################### LH #####################################

lh_out <- post_prob(file_path = "G:\\My Drive\\Data\\transposon data\\gene_level_data\\sim_data\\sim_insdens_LH.csv",
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
lh_roc <- lh_out[[1]]
lh_post_par_vals <- lh_out[[2]]
par(mfrow = c(3,2))
lab_mag <- 1.5
lab_annot <- 1.4
par(mar =  c(4, 5, 2, 1) + 0.1)
plot(lh_out[[2]]$alpha_e, xlab = "Iteration", ylab = expression(alpha[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(lh_out[[2]]$alpha_n, xlab = "Iteration", ylab = expression(alpha[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(lh_out[[2]]$beta_e, xlab = "Iteration", ylab = expression(beta[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(lh_out[[2]]$beta_n, xlab = "Iteration", ylab = expression(beta[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(lh_out[[2]]$theta, xlab = "Iteration", ylab = expression(theta),
     cex.lab = lab_mag, cex.axis = lab_annot)



####################################### LL #####################################

ll_out <- post_prob(file_path = "G:\\My Drive\\Data\\transposon data\\gene_level_data\\sim_data\\sim_insdens_LL.csv",
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
ll_roc <- ll_out[[1]]
ll_post_par_vals <- ll_out[[2]]
par(mfrow = c(3,2))
lab_mag <- 1.5
lab_annot <- 1.4
par(mar =  c(4, 5, 2, 1) + 0.1)
plot(ll_out[[2]]$alpha_e, xlab = "Iteration", ylab = expression(alpha[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(ll_out[[2]]$alpha_n, xlab = "Iteration", ylab = expression(alpha[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(ll_out[[2]]$beta_e, xlab = "Iteration", ylab = expression(beta[e]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(ll_out[[2]]$beta_n, xlab = "Iteration", ylab = expression(beta[n]),
     cex.lab = lab_mag, cex.axis = lab_annot)
plot(ll_out[[2]]$theta, xlab = "Iteration", ylab = expression(theta),
     cex.lab = lab_mag, cex.axis = lab_annot)


################################################################################
################################## plot ROC ####################################
################################################################################
roc_lty <- c(2,1,4,3)
par(mfrow = c(1, 1))
#HH
hh_all_data <- read.csv("G:\\My Drive\\Data\\transposon data\\gene_level_data\\sim_data\\sim_insdens_HH.csv")
hh_status <- HH_ess_status$status
pred <- ROCR::prediction(hh_roc$post_prob, labels = hh_status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf, xlim = c(0, 1), ylim = c(0, 1),
           xlab = (expression(bold("False positive rate"))),
           ylab = (expression(bold("True positive rate"))),
           lty = roc_lty[1], col = "black", cex.lab =1.2, cex.axis = 1.8,
           lwd = 3, xaxis.cex.axis = 1.5, yaxis.cex.axis = 1.5,
           xaxis.lwd = 2, yaxis.lwd = 3)
auc <- unlist(ROCR::performance(pred,"auc")@y.values)
hh_auc <-round(auc, 3)

# LH
lh_all_data <- read.csv("G:\\My Drive\\Data\\transposon data\\gene_level_data\\sim_data\\sim_insdens_LH.csv")
lh_status <- LH_ess_status$status
pred <- ROCR::prediction(lh_roc$post_prob, labels = lh_status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf,avg="vertical", add = T, lty = roc_lty[2],col = "black",lwd=3)
auc <- unlist(ROCR::performance(pred,"auc")@y.values)
lh_auc <-round(auc, 3)

# HL
hl_all_data <- read.csv("G:\\My Drive\\Data\\transposon data\\gene_level_data\\sim_data\\sim_insdens_HL.csv")
hl_status <- HL_ess_status$status
pred <- ROCR::prediction(hl_roc$post_prob, labels = hl_status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf,avg="vertical", add = T, lty = roc_lty[3],col = "black",lwd=3)
auc <- unlist(ROCR::performance(pred,"auc")@y.values)
hl_auc <-round(auc, 3)

# LL
ll_all_data <- read.csv("G:\\My Drive\\Data\\transposon data\\gene_level_data\\sim_data\\sim_insdens_LL.csv")
ll_status <- LL_ess_status$status
pred <- ROCR::prediction(ll_roc$post_prob, labels = ll_status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf,avg="vertical", add = T, lty = roc_lty[4],col = "black",lwd=3)
auc <- unlist(ROCR::performance(pred,"auc")@y.values)
ll_auc <-round(auc, 3)

legend(x = 0.45, y= 0.42,
       legend = c(paste(" HH, AUC = ", hh_auc, sep=""),
                  paste(" LH, AUC = ", lh_auc, sep=""),
                  paste(" HL, AUC = ", hl_auc, sep=""),
                  paste(" LL, AUC = ", ll_auc, sep="")),
       col = rep("black", 4),
       lty = roc_lty,
       pch = 10,
       pt.cex = 0,
       #text.font = 10,
       cex = 1.6,
       horiz = FALSE,
       bty = "n",
       lwd = 3,
       x.intersp = 0,
       y.intersp = 1.4)

