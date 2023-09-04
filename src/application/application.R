# Application of lpreba package

##################
#    Packages    #
##################

library(haven) # Reading Stata files
library(magrittr) # Provision of the pipe operator %>%
library(dplyr) # Grammar of data manipulation
library(Hmisc) # Provision of cutting a numeric variable into intervals
library(lpreba)

##############
#    Data    #
##############

# Source: MHE Data Archive (https://economics.mit.edu/faculty/angrist/data1/mhe)

setwd("./src/application")
data <- "./data/lee_2008.dta"

# Relevant variables for our analysis:
#
# difdemshare: Democratic vote share margin of victory, election t
# mdemwinnext: Probability of victory, election t+1

##########################
#    Data preparation    #
##########################

# Read in data and keep only relevant variables
data_lee_2008 <- read_stata(data, col_select = c(use, difdemshare, mdemwinnext)) %>%
    filter(difdemshare >= -0.25 & difdemshare <= 0.25 & use == 1)

# Create bins of width = 0.005
data_lee_2008$interval <- cut2(data_lee_2008$difdemshare, seq(-0.25, 0.25, 0.005))
data_lee_2008 <- data_lee_2008[!duplicated(data_lee_2008$interval), ]
data_lee_2008 <- data_lee_2008[order(data_lee_2008$interval), ]
data_lee_2008$difdemshare_binned <- seq(-0.25, 0.245, 0.005)

difdemshare_binned <- data_lee_2008$difdemshare_binned
mdemwinnext <- data_lee_2008$mdemwinnext

################################################################################
#                                    Figures                                   #
################################################################################

##################
#    Figure 1    #
##################

# Note:
#
# The code replicates Fig. 5a of Lee (2008), without the parametric (logit) fit.

pdf("./figures/application_figure_01.pdf", width = 14, height = 8)

plot(difdemshare_binned, mdemwinnext,
     xaxt = "n", xlab = "Democratic Vote Share Margin of Victory, Election t",
     yaxt = "n", ylab = "Probability of Victory, Election t+1",
     pch = 16, cex.lab = 1.25)
axis(1, at = seq(-0.25, 0.25, 0.05), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1.25)
abline(v = 0, lty = "dashed")
legend("topleft", legend = c("Local average (binsize = 0.005)"), pch = 16, bty = "n", cex = 1.25)

dev.off()

######################
#    Figures 3, 4    #
######################

bw_grid <- seq(0.05, 0.25, length.out = 200)

# Cross-validated bandwidths for each method (Nadaraya-Watson,
# local linear regression, boundary-adjusted Nadaraya-Watson) and side of the threshold 0
# (including plots of CV error over bandwidth grid).

pdf("./figures/application_figure_03a.pdf", width = 14, height = 8)
bw_CV_NW_lower <- bw_CV_fun(X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
                            bw_grid = bw_grid, degree = 0L)
dev.off()
pdf("./figures/application_figure_04a.pdf", width = 14, height = 8)
bw_CV_NW_upper <- bw_CV_fun(X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
                            bw_grid = bw_grid, degree = 0L)
dev.off()

pdf("./figures/application_figure_03b.pdf", width = 14, height = 8)
bw_CV_LL_lower <- bw_CV_fun(X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
                            bw_grid = bw_grid, degree = 1L)
dev.off()
pdf("./figures/application_figure_04b.pdf", width = 14, height = 8)
bw_CV_LL_upper <- bw_CV_fun(X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
                            bw_grid = bw_grid, degree = 1L)
dev.off()

# bw_CV_NW_adjusted_lower <- bw_CV_fun(X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
#                                      bw_grid = bw_grid, degree = 0L,
#                                      boundary_left = -0.25, boundary_right = -0.005,
#                                      boundary_adjustment = TRUE)
# bw_CV_NW_adjusted_upper <- bw_CV_fun(X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
#                                      bw_grid = bw_grid, degree = 0L,
#                                      boundary_left = 0, boundary_right = 0.245,
#                                      boundary_adjustment = TRUE)

# Estimates for each method and side of the threshold 0

NW_lower <- LP(x = difdemshare_binned[1:50], X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
               bw = bw_CV_NW_lower, degree = 0L)$estimates
NW_upper <- LP(x = difdemshare_binned[51:100], X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
               bw = bw_CV_NW_upper, degree = 0L)$estimates

LL_lower <- LP(x = difdemshare_binned[1:50], X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
               bw = bw_CV_LL_lower, degree = 1L)$estimates
LL_upper <- LP(x = difdemshare_binned[51:100], X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
               bw = bw_CV_LL_upper, degree = 1L)$estimates

NW_adjusted_lower <- NW_boundary(x = difdemshare_binned[1:50], X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
                                 bw = bw_CV_NW_lower,
                                 boundary_left = -0.25, boundary_right = -0.005)$estimates
NW_adjusted_upper <- NW_boundary(x = difdemshare_binned[51:100], X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
                                 bw = bw_CV_NW_upper,
                                 boundary_left = 0, boundary_right = 0.245)$estimates

##################
#    Figure 2    #
##################

pdf("./figures/application_figure_02.pdf", width = 14, height = 8)

plot(difdemshare_binned, mdemwinnext,
     xaxt = "n", xlab = "Democratic Vote Share Margin of Victory, Election t",
     yaxt = "n", ylab = "Probability of Victory, Election t+1",
     pch = 16, cex.lab = 1.25)
axis(1, at = seq(-0.25, 0.25, 0.05), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1.25)
abline(v = 0, lty = "dashed")
lines(difdemshare_binned[1:50], NW_lower, lwd = 1.5, col = "darkorange")
lines(difdemshare_binned[51:100], NW_upper, lwd = 1.5, col = "darkorange")
lines(difdemshare_binned[1:50], LL_lower, lwd = 1.5, col = "forestgreen")
lines(difdemshare_binned[51:100], LL_upper, lwd = 1.5, col = "forestgreen")
lines(difdemshare_binned[1:50], NW_adjusted_lower, lwd = 1.5, col = "chocolate4")
lines(difdemshare_binned[51:100], NW_adjusted_upper, lwd = 1.5, col = "chocolate4")
legend("topleft",
       legend = c("Local average (binsize = 0.005)", "Nadaraya-Watson", "Local linear", "NW boundary-adjusted"),
       pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 1), lwd = 1.5,
       col = c("black", "darkorange", "forestgreen", "chocolate4"),
       bty = "n", cex = 1.25)

dev.off()

##################
#    Figure 5    #
##################

# Confidence intervals (95% asymptotic)

confidence_intervals_left <- confidence_intervals_LP(x = difdemshare_binned[1:50], X = difdemshare_binned[1:50],
                                                     Y =  mdemwinnext[1:50], bw = bw_CV_LL_lower)
confidence_intervals_left_lower <- confidence_intervals_left$confidence_intervals_lower
confidence_intervals_left_upper <- confidence_intervals_left$confidence_intervals_upper

confidence_intervals_right <- confidence_intervals_LP(x = difdemshare_binned[51:100], X = difdemshare_binned[51:100],
                                                      Y = mdemwinnext[51:100], bw = bw_CV_LL_upper)
confidence_intervals_right_lower <- confidence_intervals_right$confidence_intervals_lower
confidence_intervals_right_upper <- confidence_intervals_right$confidence_intervals_upper

pdf("./figures/application_figure_05.pdf", width = 14, height = 8)

plot(difdemshare_binned, mdemwinnext,
     xaxt = "n", xlab = "Democratic Vote Share Margin of Victory, Election t",
     yaxt = "n", ylab = "Probability of Victory, Election t+1",
     pch = 16, cex.lab = 1.25)
axis(1, at = seq(-0.25, 0.25, 0.05), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1.25)
abline(v = 0, lty = "dashed")
lines(difdemshare_binned[1:50], LL_lower, lwd = 1.5, col = "forestgreen")
lines(difdemshare_binned[51:100], LL_upper, lwd = 1.5, col = "forestgreen")
polygon(x = c(difdemshare_binned[1:50], rev(difdemshare_binned[1:50])),
        y = c(confidence_intervals_left_upper, rev(confidence_intervals_left_lower)),
        col = adjustcolor("forestgreen", alpha.f = 0.1),
        border = NA)
polygon(x = c(difdemshare_binned[51:100], rev(difdemshare_binned[51:100])),
        y = c(confidence_intervals_right_upper, rev(confidence_intervals_right_lower)),
        col = adjustcolor("forestgreen", alpha.f = 0.1),
        border = NA)
legend("topleft",
       legend = c("Local average (binsize = 0.005)", "Local linear", "95% asymptotic confidence bands"),
       pch = c(16, NA, 15), lty = c(NA, 1, NA), lwd = 1.5,
       col = c("black", "forestgreen", adjustcolor("forestgreen", alpha.f = 0.2)),
       bty = "n", cex = 1.25)

dev.off()