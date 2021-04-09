set.seed(2457975)

load(paste0("G:\\My Drive\\Data\\transposon data\\gene_level_data\\",
            "data_and_code_for_insdens_paper\\ppe_and_status_real_data.RData"))

point_in_roc_space <- function(ppes, status, ca = c(99, 60, 30 , 5),
                               cb = c(1, 40, 70, 95)){
  r <- ca / (ca + cb)
  fprs <- sapply(r, function(x) sum(ppes[status == 0] > x) / sum(status == 0))
  tprs <- sapply(r, function(x) sum(ppes[status == 1] > x) / sum(status == 1))
  fpr_and_tpr <- cbind(fprs, tprs)
  return(fpr_and_tpr)
}

# e coli
pred <- ROCR::prediction(ec_ppe_and_status$post_prob,
                         labels = ec_ppe_and_status$status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf, xlim = c(0, 0.06), ylim = c(0, 1),
           xlab = (expression(bold("False positive rate"))),
           ylab = (expression(bold("True positive rate"))),
           lty = 1, cex.lab =1.2, cex.axis = 1.8,
           lwd = 1, xaxis.cex.axis = 1.5, yaxis.cex.axis = 1.5)

out <- point_in_roc_space(ppes = ec_ppe_and_status$post_prob,
                          status = ec_ppe_and_status$status)
cex <- 1.5
text(out[1, 1], out[1, 2], "99", cex = cex)
text(out[2, 1], out[2, 2], "60", cex = cex)
text(out[3, 1], out[3, 2], "30", cex = cex)
text(out[4, 1], out[4, 2], "5", cex = cex)
points(0.027, 0.91, pch = 16, cex = 2)

# s typhimurium
pred <- ROCR::prediction(stmm_ppe_and_status$post_prob,
                         labels = stmm_ppe_and_status$status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf, xlim = c(0, 0.02), ylim = c(0, 1),
           xlab = (expression(bold("False positive rate"))),
           ylab = (expression(bold("True positive rate"))),
           lty = 1, cex.lab =1.2, cex.axis = 1.8,
           lwd = 1, xaxis.cex.axis = 1.5, yaxis.cex.axis = 1.5)

out <- point_in_roc_space(ppes = stmm_ppe_and_status$post_prob,
                          status = stmm_ppe_and_status$status)

text(out[1, 1], out[1, 2], "99", cex = cex)
text(out[2, 1], out[2, 2], "60", cex = cex)
text(out[3, 1], out[3, 2], "30", cex = cex)
text(out[4, 1], out[4, 2], "5", cex = cex)
points(0.012, 0.45, pch = 16, cex = 2)

# s aureus
pred <- ROCR::prediction(sa_ppe_and_status$post_prob,
                         labels = sa_ppe_and_status$status)
perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
ROCR::plot(perf, xlim = c(0, 0.065), ylim = c(0, 1),
           xlab = (expression(bold("False positive rate"))),
           ylab = (expression(bold("True positive rate"))),
           lty = 1, cex.lab =1.2, cex.axis = 1.8,
           lwd = 1, xaxis.cex.axis = 1.5, yaxis.cex.axis = 1.5)

out <- point_in_roc_space(ppes = sa_ppe_and_status$post_prob,
                          status = sa_ppe_and_status$status)
text(out[1, 1], out[1, 2], "99", cex = cex)
text(out[2, 1], out[2, 2], "60", cex = cex)
text(out[3, 1], out[3, 2], "30", cex = cex)
text(out[4, 1], out[4, 2], "5", cex = cex)
points(0.048, 0.978, pch = 16, cex = 1.8)

