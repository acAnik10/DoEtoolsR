#' @title Completely Randomized Design
#'
#' @description Test for Significant Treatment effects and Differential Effects of Treatment Contrasts for a CRD
#'
#' @param resp The response variable vector
#' @param trt The treatment vector
#' @param los Level of significance (Default is 0.05)
#'
#' @return A list containing - ANOVA Table, Decision Table, Rejected Pairs, Mean Square Error, Critical Point
#'
#' @details In experimental designs if test of differential effect gets rejected we might be interested in analyzing which pair of treatments is behind the rejection.
#' For CRD we need to calculate the estimate of the treatment contrast (tau_i - tau_i') and test H0: tau_{i} = tau_{i'} for plausible rejection.
#' The estimate of (tau_i - tau_i') is (y_{i0} - y_{i'0}) which follows N(tau_i - tau_i', sigma^2 * (1/n_i + 1/n_i')).
#' Under H_0, the test statistic (y_{i0} - y_{i'0})/sqrt(MSE*(1/n_i + 1/n_i')) follows t_{n-v}.
#'
#' @author Anik Chakraborty
#' @section Special Thanks: Professor Surupa Chakraborty and Professor Debjit Sengupta for helping me in building the concepts of Design of Experiments.
#' Professor Madhura Dasgupta for guiding me in R programming.
#'
#' @seealso For RBD \code{\link{diff_eff.rbd}}, for LSD \code{\link{diff_eff.lsd}}
#'
#' @export diff_eff.crd

diff_eff.crd = function(resp, trt, los = 0.05)
{
     # Completely Randomized Design

     if (!is.numeric(resp))
          stop("The response variable must be a numeric vector")

     if (length(resp) != length(trt))
          stop("Reponse and Treatment vectors must be of same length")

     crd = data.frame(Response = resp,
                      Treatments = factor(trt))

     n = nrow(crd)                # Total observations
     mu_hat = mean(crd$Response)  # Grand Mean

     model = summary(stats::aov(Response ~ Treatments, crd))

     # Treatment wise Data
     trt =  dplyr::summarise(dplyr::group_by(crd, Treatments),
                             Means = mean(Response),
                             size = dplyr::n())

     y_bar = trt$Means; n_ = trt$size; v = nrow(trt)

     # MSE: Estimate of variance in the model
     MSE = model[[1]]$`Mean Sq`[2]; MSE

     # Critical value
     crit = stats::qt(los/2, n-v, lower.tail = F)

     # Initializing the decision table
     eff = data.frame(Pairs = 0, T_obs = 0, Decision = 0)
     k = 1

     for (i in seq_len(v-1))
          for (j in (i+1):v)
          {
               # Test Statistic
               Test_stat = (y_bar[i] - y_bar[j])/sqrt(MSE * (1/n_[i] + 1/n_[j]))

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
