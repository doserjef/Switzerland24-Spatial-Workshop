
# Another example of how really evil patterns in detection probability
# can screw the inferences from a detection-naive model (as shown
# in the slides in the HM intro on Monday of the course)


simNpC2 <- function (T = 20, expN = c(100, 75), dp = c(0.5, 0.5)) 
{
    T <- round(T[1])
    lambda <- seq(expN[1], expN[2], length.out = T)
    p <- seq(dp[1], dp[2], length.out = T)
    N <- rpois(T, lambda)
    C <- rbinom(T, N, p)
    return(list(T = T, expN = expN, dp = dp, lambda = lambda, 
        p = p, N = N, C = C))
}

set.seed(1)
str(dat <- simNpC2( dp = c(0.1, 0.5)))

# Quick and dirty example of CR type of analysis
Nhat <- dat$C/dat$p
x <- 1:dat$T 
fm <- lm(Nhat ~ x)
xpred <- seq(1, 20, length.out = 100)
fv <- predict(fm, newdata = data.frame (x = xpred), se = TRUE)
ci <- cbind(fv$fit - 2 * fv$se.fit,fv$fit + 2 * fv$se.fit)

# Plot 1
xlim <- c(0, 21)
par(mfrow = c(1, 3), mar = c(6, 6, 6, 4), cex.main = 2, cex.axis = 2, cex.lab = 2)
# The truth
plot(1:dat$T, dat$lambda, xlab = "Year", ylab = "", main = "True abundance (N)", xlim = xlim,
     ylim = c(0, 120), type = "l", lwd = 4, col = rgb(1,0,0,0.5), frame = FALSE, las = 1)
points(1:dat$T, dat$N, pch = 16, col = rgb(1,0,0,0.8), cex = 2)

# The observation process: "more and better people see more"
plot(1:dat$T, dat$p, xlab = "Year", ylab = "", main = "Detectability (p):\n 'more people see more'", 
     ylim = c(0, 1), type = "l", lwd = 5, col = rgb(0,0,0,0.5), frame = FALSE, las = 1, xlim = xlim,)

# The observed data: function of true N and of observation process
plot(1:dat$T, dat$C, xlab = "Year", ylab = "", main = "Observed counts (C)", xlim = xlim,
     ylim = c(0, 120), pch = 16, col = rgb(0,0,0, 0.8), cex = 2, frame = FALSE)
abline(lm(dat$C ~ I(1:dat$T)), lwd = 4, col = rgb(0,0,0,0.5)) 
lines(1:dat$T, dat$lambda, type = "l", lty = 1, lwd = 3, col = rgb(1,0,0,0.5))


# Plot 2
xlim <- c(0, 21)
par(mfrow = c(1, 3), mar = c(6, 6, 6, 4), cex.main = 2, cex.axis = 2, cex.lab = 2)
# The truth
plot(1:dat$T, dat$lambda, xlab = "Year", ylab = "", main = "True abundance (N)", xlim = xlim,
     ylim = c(0, 120), type = "l", lwd = 4, col = rgb(1,0,0,0.5), frame = FALSE, las = 1)
points(1:dat$T, dat$N, pch = 16, col = rgb(1,0,0,0.8), cex = 2)

# The observation process: "more and better people see more"
plot(1:dat$T, dat$p, xlab = "Year", ylab = "", main = "Detectability (p):\n 'more people see more'", 
     ylim = c(0, 1), type = "l", lwd = 5, col = rgb(0,0,0,0.5), frame = FALSE, las = 1, xlim = xlim,)

# The observed data: function of true N and of observation process
plot(1:dat$T, dat$C, xlab = "Year", ylab = "", main = "Observed counts (C)\n ... and estimates", xlim = xlim,
     ylim = c(0, 120), pch = 16, col = rgb(0,0,0, 0.8), cex = 2, frame = FALSE)
abline(lm(dat$C ~ I(1:dat$T)), lwd = 4, col = rgb(0,0,0,0.5)) 
lines(1:dat$T, dat$lambda, type = "l", lty = 1, lwd = 3, col = rgb(1,0,0,0.5))

# Add estimates
lines(xpred, fv$fit, col = rgb(0,0,1, 0.7), lwd = 5)
polygon(c(xpred, rev(xpred)), c(ci[,1], rev(ci[,2])), col = rgb(0,0,1, 0.2), border = NA )












