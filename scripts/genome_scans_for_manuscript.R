args <- commandArgs(trailingOnly = TRUE)
focal.phen.name <- args[1]


library(tidyr); library(dplyr); library(qtl); library(vqtl);
set.seed(27599)

# simulate the cross
my.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(100, 3), n.mar = 11, eq.spacing = TRUE, include.x = FALSE, anchor.tel = TRUE),
													 n.ind = 400,
													 type = 'f2')
my.cross$pheno$sex <- rep(x = c(0, 1), each = 200)
my.cross$pheno$batch <- sample(x = 10, size = 400, replace = TRUE)
my.cross <- qtl::calc.genoprob(cross = my.cross, step = 2)
batch.effects <- runif(n = 10, min = -0.5, max = 0.5)

# simulate phenotype1 through phenotype4, for the body of the paper
my.cross$pheno$phenotype1 <- rnorm(n = qtl::nind(my.cross),
																	 mean = 0,
																	 sd = 1)
my.cross$pheno$phenotype2 <- rnorm(n = qtl::nind(my.cross),
																	 mean = 0.28*(my.cross$geno$`1`$data[,6] - 2),
																	 sd = 1)
my.cross$pheno$phenotype3 <- rnorm(n = qtl::nind(my.cross),
																	 mean = 0,
																	 sd = exp(0.23*(my.cross$geno$`2`$data[,6] - 2)))
my.cross$pheno$phenotype4 <- rnorm(n = qtl::nind(my.cross),
																	 mean = 0.24*(my.cross$geno$`3`$data[,6] - 2),
																	 sd = exp(0.16*(my.cross$geno$`3`$data[,6] - 2)))

# simulate phenotype1x through phenotype4x, for the appendix
my.cross$pheno$phenotype1x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0,
																		sd = exp(my.cross$pheno$sex - 0.5))
my.cross$pheno$phenotype2x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0.48*(my.cross$geno$`1`$data[,6] - 2),
																		sd = exp(0.5*batch.effects[my.cross$pheno$batch]))
my.cross$pheno$phenotype3x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0,
																		sd = exp(0.24*(my.cross$geno$`2`$data[,6] - 2) + 0.5*(batch.effects[my.cross$pheno$batch])))
my.cross$pheno$phenotype4x <- 
	rnorm(n = qtl::nind(my.cross),
				mean = 0.21*(my.cross$geno$`3`$data[,6] - 2),
				sd = exp(0.12*(my.cross$geno$`3`$data[,6] - 2) + 0.2*(batch.effects[my.cross$pheno$batch])))


# do a genome scan, then permutations to assess significance
so <- qtl::scanone(cross = my.cross,
									 pheno.col = focal.phen.name,
									 addcovar = my.cross$pheno$batch)
sov <- scanonevar(cross = my.cross,
									mean.formula = formula(paste(focal.phen.name, '~ batch + mean.QTL.add + mean.QTL.dom')),
									var.formula = ~batch + var.QTL.add + var.QTL.dom)
p1 <- plot(x = sov, y = so, ylim = c(0, 6))
ggplot2::ggsave(plot = p1, height = 2.5, width = 9,
								filename = paste0('images/LOD_scan_', focal.phen.name,'.pdf'))
message('Done with so and sov')


perms <-  qtl::scanone(cross = my.cross,
											 pheno.col = focal.phen.name,
											 addcovar = my.cross$pheno$batch,
											 n.perm = 1000,
											 verbose = FALSE)
the.evd <- evd::fgev(x = perms)
so$empir.p <- evd::pgev(q = so$lod,
												loc = fitted(the.evd)[1],
												scale = fitted(the.evd)[2],
												shape = fitted(the.evd)[3],
												lower.tail = FALSE)
message('Done with so perms')

sovp <- scanonevar.perm(sov = sov, n.perms = 1000, random.seed = 27599, n.cores = 40)
p2 <- plot(x = sovp, y = so, ylim = c(0, 3.5))
ggplot2::ggsave(plot = p2, height = 2.5, width = 9,
								filename = paste0('images/empir_p_scan_', focal.phen.name,'.pdf'))
message('Done with sov perms')

