#Plot Graph for -log10(p-value) with Bonferroni
setwd("~/ResearchCode/SNPsR")
infile_path <- "result/skat/conclusion-p-value/conclusion-old-group.csv"
p_data <- read.csv(infile_path, header = TRUE)
x <- p_data[, 8]
y <- p_data[, 9]

jpeg(file = "graphic/plotminuslog10andBon2.jpeg")
plot(x,
     y, 
     main = "Plot -log10(p-value) with Bonferroni = 2.00", 
     xlab = "Bp position", 
     ylab = "-log10(p-value)", 
     ylim = c(0,3))

abline(h=2, col = "red", lwd = 2)
abline(h=1.3, col = "blue", lwd = 2)
dev.off()

#Plot Nonparametrix Regression

library(splines)
ss <- seq(x[1]+1, x[length(y)]-1,by=50)
stry = lm(y~ns (x,knots ==ss))
lines(x,stry$fitted, col ="purple")



