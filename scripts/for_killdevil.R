args <- commandArgs(trailingOnly = TRUE)
phen <- args[1]

if (!any(phen == c('phen1', 'phen2', 'phen3', 'phen4'))) {
	stop("phen can only be phen1 - phen4")
}

library(tidyr); library(dplyr); library(qtl); library(vqtl);

set.seed(27599)

my.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(100, 3), n.mar = 11, eq.spacing = TRUE, include.x = FALSE, anchor.tel = TRUE),
													 n.ind = 400,
													 type = 'f2')
my.cross$pheno$sex <- rep(x = c(0, 1), each = 200)
my.cross <- qtl::calc.genoprob(cross = my.cross, step = 2)

my.cross$pheno$phenotype1x <- rnorm(n = qtl::nind(my.cross))

my.cross$pheno$phenotype2x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0.26*(my.cross$geno$`1`$data[,6] - 2),
																		sd = exp(my.cross$pheno$sex - 0.5))

my.cross$pheno$phenotype3x <- rnorm(n = qtl::nind(my.cross),
																		sd = exp(0.21*(my.cross$geno$`2`$data[,6] - 2) + (my.cross$pheno$sex - 0.5)))

my.cross$pheno$phenotype4x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0.2*(my.cross$geno$`3`$data[,6] - 2),
																		sd = exp(0.15*(my.cross$geno$`3`$data[,6] - 2) + (my.cross$pheno$sex - 0.5)))

n.perms <- 500
n.cores <- 10

if (phen == 'phen1') {
	# make fig 1 -- LOD score scans
	ymax <- 6
	a0 <- qtl::scanone(cross = my.cross,
										 pheno.col = 'phenotype1x')
	a1 <- scanonevar(cross = my.cross,
									 mean.formula = phenotype1x ~ sex + mean.QTL.add + mean.QTL.dom,
									 var.formula = ~sex + var.QTL.add + var.QTL.dom)
	plot(pa1 <- plot(x = a1, y = a0, ylim = c(0, ymax)))
	ggplot2::ggsave(plot = pa1, filename = 'images/LOD_scan_phen1x.pdf', height = 2.5, width = 9)


	# make fig 2 -- empirircal p-value scans
	seed <- 27599
	ymax <- 5
	system.time(a2 <- scanonevar.perm(sov = a1, n.perms = n.perms, random.seed = seed, n.cores = n.cores))
	plot(pa2 <- plot(x = a2, ylim = c(0, ymax)))
	ggplot2::ggsave(plot = pa2, filename = 'images/empir_p_scan_phen1x.pdf', height = 2.5, width = 9)
}


if (phen == 'phen2') {
	# make fig 1 -- LOD score scans
	ymax <- 6
	b0 <- qtl::scanone(cross = my.cross,
										 pheno.col = 'phenotype2x')
	b1 <- scanonevar(cross = my.cross,
									 mean.formula = phenotype2x ~ sex + mean.QTL.add + mean.QTL.dom,
									 var.formula = ~sex + var.QTL.add + var.QTL.dom)
	plot(pb1 <- plot(x = b1, y = b0, ylim = c(0, ymax)))
	ggplot2::ggsave(plot = pb1,	filename = 'images/LOD_scan_phen2x.pdf', height = 2.5, width = 9)


	# make fig 2 -- empirircal p-value scans
	seed <- 27599
	ymax <- 5
	system.time(b2 <- scanonevar.perm(sov = b1, n.perms = n.perms, random.seed = seed, n.cores = n.cores))
	plot(pb2 <- plot(x = b2, ylim = c(0, ymax)))
	ggplot2::ggsave(plot = pb2, filename = 'images/empir_p_scan_phen2x.pdf', height = 2.5, width = 9)
}

if (phen == 'phen3') {
	# make fig 1 -- LOD score scans
	ymax <- 6
	c0 <- qtl::scanone(cross = my.cross,
										 pheno.col = 'phenotype3x')
	c1 <- scanonevar(cross = my.cross,
									 mean.formula = phenotype3x ~ sex + mean.QTL.add + mean.QTL.dom,
									 var.formula = ~sex + var.QTL.add + var.QTL.dom)
	plot(pc1 <- plot(x = c1, y = c0, ylim = c(0, ymax)))
	ggplot2::ggsave(plot = pc1, filename = 'images/LOD_scan_phen3x.pdf', height = 2.5, width = 9)


	# make fig 2 -- empirircal p-value scans
	seed <- 27599
	ymax <- 5
	system.time(c2 <- scanonevar.perm(sov = c1, n.perms = n.perms, random.seed = seed, n.cores = n.cores))
	plot(pc2 <- plot(x = c2, ylim = c(0, ymax)))
	ggplot2::ggsave(plot = pc2, filename = 'images/empir_p_scan_phen3x.pdf', height = 2.5, width = 9)
}


if (phen == 'phen4') {
	# make fig 1 -- LOD score scans
	ymax <- 6
	d0 <- qtl::scanone(cross = my.cross,
										 pheno.col = 'phenotype4x')
	d1 <- scanonevar(cross = my.cross,
									 mean.formula = phenotype4x ~ sex + mean.QTL.add + mean.QTL.dom,
									 var.formula = ~sex + var.QTL.add + var.QTL.dom)
	plot(pd1 <- plot(x = d1, y = d0, ylim = c(0, ymax)))
	ggplot2::ggsave(plot = pd1, filename = 'images/LOD_scan_phen4x.pdf', height = 2.5, width = 9)


	# make fig 2 -- empirircal p-value scans
	seed <- 27599
	ymax <- 5
	system.time(d2 <- scanonevar.perm(sov = d1, n.perms = n.perms, random.seed = seed, n.cores = n.cores))
	plot(pd2 <- plot(x = d2, ylim = c(0, ymax)))
	ggplot2::ggsave(plot = pd2, filename = 'images/empir_p_scan_phen4x.pdf', height = 2.5, width = 9)
}
