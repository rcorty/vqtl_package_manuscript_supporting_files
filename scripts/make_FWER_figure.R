library(tidyr); library(dplyr); library(vqtl)

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
	if (missing(col)) { stop("Please provide a vector of colours.") }
	apply(sapply(col, col2rgb)/255, 2,
				function(x)
					rgb(x[1], x[2], x[3], alpha = alpha))
}

set.seed(2)

my.cross <- sim.cross(map = sim.map(len = rep(100, 4), n.mar = 30, eq.spacing = TRUE, include.x = FALSE),
											n.ind = 200,
											type = 'f2')
my.cross$pheno$sex <- rep(x = c(0, 1), each = 100)
my.cross <- calc.genoprob(my.cross)

my.cross$pheno$phenotype1 <- rnorm(n = nind(my.cross))
my.cross$pheno$phenotype2 <- rnorm(n = nind(my.cross), mean = 0.8*my.cross$geno$`1`$data[,15])
my.cross$pheno$phenotype3 <- rnorm(n = nind(my.cross), sd = my.cross$geno$`2`$data[,15])
my.cross$pheno$phenotype4 <- rnorm(n = nind(my.cross), mean = my.cross$geno$`3`$data[,15], sd = my.cross$geno$`3`$data[,15])


a1 <- scanonevar.perm(cross = my.cross,
											mean.formula = 'phenotype2 ~ sex + mean.QTL.add + mean.QTL.dom',
											var.formula = '~sex + var.QTL.add + var.QTL.dom',
											n.perms = 5,
											verbose.return = TRUE)
a1max <- a1 %>% group_by(perm) %>% summarise(max.full.lod = max(full.lod),
																						 max.mean.lod = max(mean.lod),
																						 max.var.lod = max(var.lod))

a2 <- list()
for (i in 1:5) {
	my.cross$pheno$phenotype2 <- sample(my.cross$pheno$phenotype2)
	a2[[i]] <- scanone(cross = my.cross,
										 pheno.col = 'phenotype2')
	a2[[i]][['perm']] <- i
}
a2 <- rbind_all(a2)







# pdf(file = '../2016_G3_PackageVQTL_Corty/images/FWER.pdf', width = 6, height = 8)
par(mar = c(2.1, 3.1, 3.1, 0.1), xpd = FALSE)
layout(mat = matrix(data = c(1:20), ncol = 5), widths = c(7, 3/4, 3/4, 3/4))
ylim <- c(0, 5.2)

p1 <- a1 %>% filter(perm == 1)
class(p1) <- c('scanonevar', 'tbl_df', 'tbl', 'data.frame')
p2 <- a2 %>% filter(perm == 1)
plot(x =  p1, y = p2,
		 show.equations = FALSE, legend.pos = NA, ylim = ylim, title.cex = 1.2, title = 'randomization 1 of phenotype 2', suppress.chromosome = TRUE)
abline(h = c(max(p1$full.lod), max(p1$mean.lod), max(p1$var.lod), max(p2$lod)), col = c('black', 'blue', 'red', 'forestgreen'), lwd = 2)


p3 <- a1 %>% filter(perm == 2)
class(p3) <- c('scanonevar', 'tbl_df', 'tbl', 'data.frame')
p4 <- a2 %>% filter(perm == 2)
plot(x =  p3, y = p4,
		 show.equations = FALSE, legend.pos = NA, ylim = ylim, title.cex = 1.2, title = 'randomization 2 of phenotype 2', suppress.chromosome = TRUE)
abline(h = c(max(p3$full.lod), max(p3$mean.lod), max(p3$var.lod), max(p4$lod)), col = c('black', 'blue', 'red', 'forestgreen'), lwd = 2)


p5 <- a1 %>% filter(perm == 3)
class(p5) <- c('scanonevar', 'tbl_df', 'tbl', 'data.frame')
p6 <- a2 %>% filter(perm == 3)
plot(x =  p5, y = p6,
		 show.equations = FALSE, legend.pos = NA, ylim = ylim, title.cex = 1.2, title = 'randomization 3 of phenotype 2', suppress.chromosome = TRUE)
abline(h = c(max(p5$full.lod), max(p5$mean.lod), max(p5$var.lod), max(p6$lod)), col = c('black', 'blue', 'red', 'forestgreen'), lwd = 2)

# three dots in bottom left corner to indicate "more perms follow..."
par(mar = c(2, 0, 0, 0))
plot(x = rep(0.5, 3), y = c(0.1, 0.5, 0.9), axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = 16, cex = 1.5)




# first column of right side, drops down into DGLM-joint histogram
par(mar = c(2.1, 0, 3.1, 0), xpd = NA)
plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.4, x1 = 0.4,
				 y0 = max(p1$full.lod), y1 = -23, col = add.alpha(col = 'black', alpha = 0.6))
segments(x0 = rep(0, 4), x1 = c(0.4, 1, 1, 1),
				 y0 = c(max(p1$full.lod), max(p1$mean.lod), max(p1$var.lod), max(p2$lod)),
				 y1 = c(max(p1$full.lod), max(p1$mean.lod), max(p1$var.lod), max(p2$lod)),
				 col = add.alpha(col = c('black', 'blue', 'red', 'forestgreen'), alpha = 0.6))


plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.5, x1 = 0.5,
				 y0 = max(p3$full.lod), y1 = -13, col = add.alpha(col = 'black', alpha = 0.6))
segments(x0 = rep(0, 4), x1 = c(0.5, 1, 1, 1),
				 y0 = c(max(p3$full.lod), max(p3$mean.lod), max(p3$var.lod), max(p4$lod)),
				 y1 = c(max(p3$full.lod), max(p3$mean.lod), max(p3$var.lod), max(p4$lod)),
				 col = add.alpha(col = c('black', 'blue', 'red', 'forestgreen'), alpha = 0.6))

plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.6, x1 = 0.6,
				 y0 = max(p5$full.lod), y1 = -3, col = add.alpha(col = 'black', alpha = 0.6))
segments(x0 = rep(0, 4), x1 = c(0.6, 1, 1, 1),
				 y0 = c(max(p5$full.lod), max(p5$mean.lod), max(p5$var.lod), max(p6$lod)),
				 y1 = c(max(p5$full.lod), max(p5$mean.lod), max(p5$var.lod), max(p6$lod)),
				 col = add.alpha(col = c('black', 'blue', 'red', 'forestgreen'), alpha = 0.6))

hist(x = rnorm(n = 1000), main = NA, axes = FALSE, ylab = NA, xlab = NA, col = add.alpha(col = 'black', alpha = 0.6))
mtext(text = 'mvQTL', side = 1, cex = 0.8)




# second column of right side, drops down into DGLM-mean histogram
plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.4, x1 = 0.4,
				 y0 = max(p1$mean.lod), y1 = -23, col = add.alpha(col = 'blue', alpha = 0.6))
segments(x0 = rep(0, 3), x1 = c(0.4, 1, 1),
				 y0 = c(max(p1$mean.lod), max(p1$var.lod), max(p2$lod)),
				 y1 = c(max(p1$mean.lod), max(p1$var.lod), max(p2$lod)),
				 col = add.alpha(col = c('blue', 'red', 'forestgreen'), alpha = 0.6))

plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.5, x1 = 0.5,
				 y0 = max(p3$mean.lod), y1 = -13, col = add.alpha(col = 'blue', alpha = 0.6))
segments(x0 = rep(0, 3), x1 = c(0.5, 1, 1),
				 y0 = c(max(p3$mean.lod), max(p3$var.lod), max(p4$lod)),
				 y1 = c(max(p3$mean.lod), max(p3$var.lod), max(p4$lod)),
				 col = add.alpha(col = c('blue', 'red', 'forestgreen'), alpha = 0.6))

plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.6, x1 = 0.6,
				 y0 = max(p5$mean.lod), y1 = -3, col = add.alpha(col = 'blue', alpha = 0.6))
segments(x0 = rep(0, 3), x1 = c(0.6, 1, 1),
				 y0 = c(max(p5$mean.lod), max(p5$var.lod), max(p6$lod)),
				 y1 = c(max(p5$mean.lod), max(p5$var.lod), max(p6$lod)),
				 col = add.alpha(col = c('blue', 'red', 'forestgreen'), alpha = 0.6))

hist(x = rnorm(n = 1000), main = NA, axes = FALSE, ylab = NA, xlab = NA, col = add.alpha(col = 'blue', alpha = 0.6))
mtext(text = 'mQTL', side = 1, cex = 0.8)





# third column of right side, drops down into DGLM-var histogram
plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.4, x1 = 0.4,
				 y0 = max(p1$var.lod), y1 = -23, col = add.alpha(col = 'red', alpha = 0.6))
segments(x0 = rep(0, 2), x1 = c(0.4, 1),
				 y0 = c(max(p1$var.lod), max(p2$lod)),
				 y1 = c(max(p1$var.lod), max(p2$lod)),
				 col = add.alpha(col = c('red', 'forestgreen'), alpha = 0.6))

plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.5, x1 = 0.5,
				 y0 = max(p3$var.lod), y1 = -13, col = add.alpha(col = 'red', alpha = 0.6))
segments(x0 = rep(0, 2), x1 = c(0.5, 1),
				 y0 = c(max(p3$var.lod), max(p4$lod)),
				 y1 = c(max(p3$var.lod), max(p4$lod)),
				 col = add.alpha(col = c('red', 'forestgreen'), alpha = 0.6))

plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.6, x1 = 0.6,
				 y0 = max(p5$var.lod), y1 = -3, col = add.alpha(col = 'red', alpha = 0.6))
segments(x0 = rep(0, 2), x1 = c(0.6, 1),
				 y0 = c(max(p5$var.lod), max(p6$lod)),
				 y1 = c(max(p5$var.lod), max(p6$lod)),
				 col = add.alpha(col = c('red', 'forestgreen'), alpha = 0.6))

hist(x = rnorm(n = 1000), main = NA, axes = FALSE, ylab = NA, xlab = NA, col = add.alpha(col = 'red', alpha = 0.6))
mtext(text = 'vQTL', side = 1, cex = 0.8)








# third column of right side, drops down into LM histogram
plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.4, x1 = 0.4,
				 y0 = max(p2$lod), y1 = -23, col = add.alpha(col = 'forestgreen', alpha = 0.6))
segments(x0 = 0, x1 = 0.4,
				 y0 = c(max(p2$lod)),
				 y1 = c(max(p2$lod)),
				 col = add.alpha(col = c('forestgreen'), alpha = 0.6))

plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.5, x1 = 0.5,
				 y0 = max(p4$lod), y1 = -13, col = add.alpha(col = 'forestgreen', alpha = 0.6))
segments(x0 = 0, x1 = c(0.5),
				 y0 = c(max(p4$lod)),
				 y1 = c(max(p4$lod)),
				 col = add.alpha(col = c('forestgreen'), alpha = 0.6))

plot(x = NA, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = NA, ylab = NA)
segments(x0 = 0.6, x1 = 0.6,
				 y0 = max(p6$lod), y1 = -3, col = add.alpha(col = 'forestgreen', alpha = 0.6))
segments(x0 = 0, x1 = c(0.6),
				 y0 = c(max(p6$lod)),
				 y1 = c(max(p6$lod)),
				 col = add.alpha(col = c('forestgreen'), alpha = 0.6))

hist(x = rnorm(n = 1000), main = NA, axes = FALSE, ylab = NA, xlab = NA, col = add.alpha(col = 'forestgreen', alpha = 0.6))
mtext(text = 'QTL', side = 1, cex = 0.8)


# dev.off()
