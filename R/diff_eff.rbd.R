#' @title Randomized Block Design
#'
#' @description Test for Significant Treatment effects and Differential Effects of Treatment Contrasts for an RBD
#'
#' @param resp The response variable vector
#' @param trt The treatment vector
#' @param blocks The vector containing Block/Replicate indices
#' @param los Level of significance (Default is 0.05)
#'
#' @return A list containing - ANOVA Table, Decision Table, Rejected Pairs, Mean Square Error, Critical Point
#'
#' @details In experimental designs if test of differential effect gets rejected we might be interested in analyzing which pair of treatments is behind the rejection.
#' For RBD we need to calculate the estimate of the treatment contrast (tau_j - tau_j') and test H0: tau_{j} = tau_{j'} for plausible rejection.
#' The estimate of (tau_j - tau_j') is (y_{0j} - y_{0j'}) which follows N(tau_j - tau_j', 2*sigma^2/b).
#' Under H_0, the test statistic (y_{0j} - y_{0j'})/sqrt(2*MSE/b) follows t_{(b-1)(v-1)}.
#'
#' @author Anik Chakraborty
#' @section Special Thanks: Professor Surupa Chakraborty and Professor Debjit Sengupta for helping me in building the concepts of Design of Experiments.
#' Professor Madhura Dasgupta for guiding me in R programming.
#'
#' @seealso For CRD \code{\link{diff_eff.crd}}, for LSD \code{\link{diff_eff.lsd}}
#'
#' @export diff_eff.rbd

diff_eff.rbd = function(resp, trt, blocks, los = 0.05)
{
     # Randomized Block Design

     if (!is.numeric(resp))
          stop("The response variable must be a numeric vector")

     l = c(length(resp), length(trt), length(blocks))

     if (any(diff(l) != numeric(2)))
          stop("Response, Treatments and Block/Replicate vector must have same lenght")

     rbd = data.frame(Response = resp,
                       Treatments = factor(trt),
                       Blocks = factor(blocks))

     n = nrow(rbd)                # Total observations
     mu_hat = mean(rbd$Response)  # Grand Mean

     model = summary(stats::aov(Response ~ Blocks + Treatments,
                                rbd))

     # Treatment wise Data
     trt = dplyr::summarise(dplyr::group_by(rbd, Treatments),
                            Means = mean(Response))

     y_bar = trt$Means

     b = length(unique(rbd$Blocks)); v = nrow(trt)

     # MSE: Estimate of variance in the model
     MSE = model[[1]]$`Mean Sq`[3]

     # Critical value
     crit = stats::qt(los/2, (b-1)*(v-1), lower.tail = F)

     # Initializing the decision table
     eff = data.frame(Pairs = 0, T_obs = 0, Decision = 0)
     k = 1

     for (i in seq_len(v-1))
          for (j in (i+1):v)
          {
               # Test Statistic
               Test_stat = (y_bar[i] - y_bar[j])/sqrt(2 * MSE/b)

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
