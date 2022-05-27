#' @title Latin Square Design
#'
#' @description Test for Significant Treatment effects and Differential Effects of Treatment Contrasts for a LSD
#'
#' @param resp The response variable vector
#' @param trt The treatment vector
#' @param rows A vector containing the row indices of a LSD
#' @param cols A vector containing the column indices of a LSD
#' @param los Level of significance (Default is 0.05)
#'
#' @return A list containing - ANOVA Table, Decision Table, Rejected Pairs, Mean Square Error, Critical Point
#'
#' @details In experimental designs if test of differential effect gets rejected we might be interested in analyzing which pair of treatments is behind the rejection.
#' For LSD we need to calculate the estimate of the treatment contrast (tau_k - tau_k') and test H0: tau_{k} = tau_{k'} for plausible rejection.
#' The estimate of (tau_k - tau_k') is (y_{00k} - y_{00k'}) which follows N(tau_j - tau_j', 2*sigma^2/v).
#' Under H_0, the test statistic (y_{00k} - y_{00k'})/sqrt(2*MSE/v) follows t_{(v-1)(v-2)}.
#'
#' @author Anik Chakraborty
#' @section Special Thanks: Professor Surupa Chakraborty and Professor Debjit Sengupta for helping me in building the concepts of Design of Experiments.
#' Professor Madhura Dasgupta for guiding me in R programming.
#
#' @seealso For CRD \code{\link{diff_eff.crd}}, for RBD \code{\link{diff_eff.rbd}}
#'
#' @export diff_eff.lsd

diff_eff.lsd = function(resp, trt, rows, cols, los = 0.05)
{
     # Latin Square Design

     if (!is.numeric(resp))
          stop("The response variable must be a numeric vector")

     l = c(length(resp), length(trt), length(rows), length(cols))

     if (any(diff(l) != numeric(3)))
          stop("Response, Treatments, Rows and Column vector must have same lenght")

     lsd = data.frame(Response = resp,
                       Treatments = factor(trt),
                       Rows = factor(rows),
                       Columns = factor(cols))

     n = nrow(lsd)                # Total observations
     mu_hat = mean(lsd$Response)  # Grand Mean

     model = summary(stats::aov(Response ~ Rows + Columns +
                                     Treatments, lsd))

     # Treatment wise Data
     trt = dplyr::summarise(dplyr::group_by(lsd, Treatments),
                            Means = mean(Response))

     y_bar = trt$Means; v = nrow(trt)

     # MSE: Estimate of variance in the model
     MSE = model[[1]]$`Mean Sq`[4]

     # Critical value
     crit = stats::qt(los/2, (v-1)*(v-2), lower.tail = F)

     # Initializing the decision table
     eff = data.frame(Pairs = 0, T_obs = 0, Decision = 0)
     k = 1

     for (i in seq_len(v-1))
          for (j in (i+1):v)
          {
               # Test Statistic
               Test_stat = (y_bar[i] - y_bar[j])/sqrt(2 * MSE/v)

               if (abs(Test_stat) > crit)
               {
                    dec = "***"
               } else dec = "-"

               eff[k,] = c(paste0("(", i, ",", j, ")"),
                           signif(Test_stat, 3), dec)

               k = k+1
          }

     # Rejected Pairs
     rej_pair = dplyr::select(dplyr::filter(eff, Decision == "***"), Pairs)

     # Output table
     output = list(`ANOVA Table` = model,
                   Means = as.data.frame(trt),
                   Critical_Value = paste("The critical value for the pairwise test:", signif(crit, 3)),
                   Decision_Table = eff,
                   Rejected_pairs = rej_pair,
                   `No. of rejected pairs` = nrow(rej_pair))

     return(output)
}
