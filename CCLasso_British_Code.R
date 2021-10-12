
# CCLasso Installation, Example, and British Species and Genus Runs
#
# Gancz 4/21/2021

############################################################################################################################

# Installation

############################################################################################################################

#install gtools 
#install.packages("gtools")
#install.packages("corrplot")
#install.packages("tidyverse")
library(gtools)
library(ggplot2)
library(dplyr)
library(corrplot)
library(readxl)
library(tidyverse)
library(corrplot)
library(hrbrthemes)
library(qgraph)


#Create SparrCC.Count Function

SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # dimension for w (latent variables)
  p <- ncol(x);
  n <- nrow(x);
  # posterior distribution (alpha)
  x <- x + 1;
  # store generate data
  y <- matrix(0, n, p);
  # store covariance/correlation matrix
  cov.w <- cor.w <- matrix(0, p, p);
  indLow <- lower.tri(cov.w, diag = T);
  # store covariance/correlation for several posterior samples
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax);
  for(i in 1:imax) {
    # generate fractions from posterior distribution
    y <- t(apply(x, 1, function(x) 
      gtools::rdirichlet(n = 1, alpha = x)));
    # estimate covariance/correlation
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin);
    # store variance/correlation only low triangle 
    covs[, i] <- cov_cor$cov.w[indLow];
    cors[, i] <- cov_cor$cor.w[indLow];
  }
  # calculate median for several posterior samples
  cov.w[indLow] <- apply(covs, 1, median); 
  cor.w[indLow] <- apply(cors, 1, median);
  #
  cov.w <- cov.w + t(cov.w);
  diag(cov.w) <- diag(cov.w) / 2;
  cor.w <- cor.w + t(cor.w);
  diag(cor.w) <- 1;
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}

# Create SparCC Frac Function 
#-------------------------------------------------------------------------------

SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # Log transformation
  x <- log(x);
  p <- ncol(x);
  # T0 = var(log(xi/xj)) variation matrix
  TT <- stats::var(x);
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT;
  # Variance and correlation coefficients for Basic SparCC  
  rowT0 <- rowSums(T0);
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2))/(p - 2);
  var.w[var.w < Vmin] <- Vmin;
  #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
  #  sqrt(outer(var.w, var.w, "*")) / 2;
  Is <- sqrt(1/var.w);
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5;
  # Truncated correlation in [-1, 1]
  cor.w[cor.w <= - 1] <- - 1; 
  cor.w[cor.w >= 1] <- 1;
  # Left matrix of estimation equation
  Lmat <- diag(rep(p - 2, p)) + 1; 
  # Remove pairs
  rp <- NULL;
  # Left components
  cp <- rep(TRUE, p);
  # Do loops until max iteration or only 3 components left
  k <- 0;  
  while(k < kmax && sum(cp) > 3) {
    # Left T0 = var(log(xi/xj)) after removing pairs
    T02 <- T0;
    # Store current correlation to find the strongest pair
    curr_cor.w <- cor.w;
    # Remove diagonal
    diag(curr_cor.w) <- 0;
    # Remove removed pairs
    if(!is.null(rp)) {
      curr_cor.w[rp] <- 0;
    }
    # Find the strongest pair in vector form
    n_rp <- which.max(abs(curr_cor.w));
    # Remove the pair if geater than alpha
    if(abs(curr_cor.w[n_rp]) >= alpha) {
      # Which pair in matrix form
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)));
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1;
      # Update remove pairs
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp);
      rp <- c(rp, n_rp);
      # Update T02
      T02[rp] <- 0;
      # Which component left
      cp <- (diag(Lmat) > 0);
      # Update variance and truncated lower by Vmin
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]));
      var.w[var.w <= Vmin] <- Vmin;
      # Update correlation matrix and truncated by [-1, 1]
      #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
      #  sqrt(outer(var.w, var.w, "*")) / 2;    
      Is <- sqrt(1/var.w);
      cor.w <- (var.w + rep(var.w, each = p) - T0) * 
        Is * rep(Is, each = p) * 0.5;
      # Truncated correlation in [-1, 1]
      cor.w[cor.w <= - 1] <- - 1;
      cor.w[cor.w >= 1] <- 1;
    }
    else {
      break;
    }
    # 
    k <- k + 1;
  }
  # Covariance
  Is <- sqrt(var.w);
  cov.w <- cor.w * Is * rep(Is, each = p);
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}

#Install CC Lasso and friends
#-------------------------------------------------------------------------------

cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
                    lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) {
  n <- nrow(x);
  p <- ncol(x);
  
  if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  x <- log(x);
  vx2 <- stats::var(x);
  
  # Diagonal weight for loss function
  rmean_vx2 <- rowMeans(vx2);
  wd <- 1/diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2));
  wd2 <- sqrt(wd);
  
  # Some global parameters for optimization with single lambda
  rho <- 1;
  u_f <- eigen(diag(p) - 1/p)$vectors;
  wd_u <- (t(u_f) %*% (wd * u_f))[-p, -p];
  wd_u_eig <- eigen(wd_u);
  d0_wd <- 1 / ( (rep(wd_u_eig$values, each = p-1) + wd_u_eig$values) / 
                   (2 * rho) + 1 );
  u0_wd <- wd_u_eig$vectors;
  
  # Golden section method for the selection of lambda (log10 scale)
  sigma <- vx2;
  lam_int2 <- log10(range(lam_int));
  a1 <- lam_int2[1]; 
  b1 <- lam_int2[2];
  # Store lambda and corresponding cross validation's loss
  lams <- NULL; 
  fvals <- NULL;
  # Two trial points in first 
  a2 <- a1 + 0.382 * (b1 - a1); 
  b2 <- a1 + 0.618 * (b1 - a1);
  fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
                         sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
                         wd2 = wd2);
  lams <- c(lams, b2); 
  fvals <- c(fvals, fb2$cv_loss);
  fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
                         sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
                         wd2 = wd2);
  lams <- c(lams, a2); 
  fvals <- c(fvals, fa2$cv_loss);
  # Error tolerance for convergence
  err_lam2 <- 1e-1 * max(1, lam_int2);
  err_fval <- 1e-4;
  
  err <- b1 - a1;
  k <- 0;
  while(err > err_lam2 && k < k_max) {
    fval_max <- max(fa2$cv_loss, fb2$cv_loss);
    
    if(fa2$cv_loss > fb2$cv_loss) {
      a1 <- a2;      
      a2 <- b2;
      fa2 <- fb2;
      b2 <- a1 + 0.618 * (b1 - a1);
      fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
                             sigma = fa2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
                             wd2 = wd2);
      
      lams <- c(lams, b2); 
      fvals <- c(fvals, fb2$cv_loss);
    } else {
      b1 <- b2;
      b2 <- a2;
      fb2 <- fa2;
      a2 <- a1 + 0.382 * (b1 - a1);
      fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
                             sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
                             wd2 = wd2);
      
      lams <- c(lams, a2); 
      fvals <- c(fvals, fa2$cv_loss);
    }
    fval_min <- min(fa2$cv_loss, fb2$cv_loss);      
    
    k <- k + 1;
    err <- b1 - a1;
    if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break;
    }
  }
  info_cv <- list(lams = lams, fvals = fvals, k = k + 2, 
                  lam_int = 10^c(a1, b1)); 
  if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
    cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",", 
        10^b1, "])\n", sep = "");
  }
  
  lambda <- 10^((a2 + b2)/2);
  # Bootstrap for cclasso
  lambda2 <- lambda / rho;
  info_boot <- boot_cclasso(x = x, sigma = fb2$sigma, lambda2 = lambda2, 
                            n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  
  return(list(var_w = info_boot$var_w, cor_w = info_boot$cor_w, 
              p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv));
}
#-------------------------------------------------------------------------------
# Bootstrap for cclasso
boot_cclasso <- function(x, sigma, lambda2, n_boot = 20, 
                         wd, u_f, u0_wd, d0_wd) {
  n <- nrow(x);
  p <- ncol(x);
  
  # Store the result of bootstrap
  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
  cors_mat <- matrix(0, p, p);
  ind_low <- lower.tri(cors_mat);
  
  # Bootstrap procedure
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
                     ncol = n_boot);
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k];
    sigma2 <- cclasso_sub(sigma = sigma, vx = var(x[ind_samp, ]), 
                          lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    
    vars_boot[, k] <- diag(sigma2);
    Is <- 1 / sqrt(vars_boot[, k]);
    cors_mat <- Is * sigma2 * rep(Is, each = p);
    cors_boot[, k] <- cors_mat[ind_low];
  }
  Is <- 1 / sqrt(diag(sigma));
  cors_mat <- sigma * Is * rep(Is, each = p);
  cors_boot[, n_boot + 1] <- cors_mat[ind_low];
  vars_boot[, n_boot + 1] <- diag(sigma);
  
  #----------------------------------------  
  # Variance estimation via bootstrap
  vars2 <- rowMeans(vars_boot);
  #----------------------------------------
  # Correlations' relationship for artificial null sample
  tol_cor <- 1e-3;
  sam_art0 <- matrix(rnorm(n * p), nrow = n) * rep(sqrt(vars2), each = n);
  cors_art0 <- cor(sam_art0)[ind_low];
  sam_art <- sam_art0 - log(rowSums(exp(sam_art0)));
  sigma_art <- cclasso_sub(sigma = sigma, vx = var(sam_art), 
                           lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  Is <- 1 / sqrt(diag(sigma_art));
  cors_mat <- Is * sigma_art * rep(Is, each = p);
  cors_art2 <- cors_mat[ind_low];
  # Bias of estimation between absolute data and cclasso of compositional data
  cors0m2 <- log( ((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) * 
                                                           (1 - cors_art0)) );
  tmp <- abs(cors_art2) >= tol_cor;
  bias02 <- ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0);
  # Modification of estimation for cclasso
  cors2 <- log( (1 + cors_boot) / (1 - cors_boot) );    
  cors2mod <- (cors_boot >= tol_cor) * (cors2 + bias02) + 
    (cors_boot <= -tol_cor) * (cors2 - bias02);
  cors2mod <- 1 - rowMeans(2 / (exp(cors2mod) + 1));
  cors2_mat <- diag(p);
  cors2_mat[ind_low] <- cors2mod;
  cors2_mat <- t(cors2_mat);
  cors2_mat[ind_low] <- cors2mod;
  # P-values with null distribution of correlation estimations of absolute data
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
  pval_mat <- diag(p);
  pval_mat[ind_low] <- p_vals;
  pval_mat <- t(pval_mat);
  pval_mat[ind_low] <- p_vals;
  #----------------------------------------
  
  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));    
}
#-------------------------------------------------------------------------------
# cross validation's loss of cclasso for single lambda
cv_loss_cclasso <- function(lambda2, x, k_cv, sigma, 
                            wd, u_f, u0_wd, d0_wd, wd2) {
  n <- nrow(x);
  p <- ncol(x);
  
  n_b <- floor(n / k_cv);
  cv_loss <- 0;
  for(k in 1:k_cv) {
    itest <- (n_b * (k - 1) + 1):(n_b * k);
    vxk <- stats::var(x[itest, ]);
    vx2k <- stats::var(x[-itest, ]);
    
    sigma <- cclasso_sub(sigma = sigma, vx = vx2k, lambda2 = lambda2, 
                         wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    
    dsig <- sigma - vxk;
    tmp <- rowMeans(dsig);
    dsig <- dsig - tmp - rep(tmp, each = p) + mean(tmp);
    cv_loss <- cv_loss + base::norm(wd2 * dsig, "F")^2;
  }
  
  return(list(cv_loss = cv_loss, sigma = sigma));
}
#-------------------------------------------------------------------------------
# cclasso for single lambda
cclasso_sub <- function(sigma, vx, lambda2, 
                        wd, u_f, u0_wd, d0_wd, 
                        k_max = 200, x_tol = 1e-4) {
  p <- ncol(sigma);
  sigma2 <- sigma;
  LAMBDA <- matrix(0, p, p);
  lambda2 <- matrix(lambda2, p, p);
  diag(lambda2) <- 0;
  
  k <- 0;
  err <- 1;
  while(err > x_tol && k < k_max) {
    # Update sigma
    x_sigma <- t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f;
    x_sigma[-p,-p] <- u0_wd %*% ((t(u0_wd) %*% x_sigma[-p, -p] %*% u0_wd) * 
                                   d0_wd) %*% t(u0_wd);
    sigma_new <- vx + u_f %*% x_sigma %*% t(u_f);
    # Update sigma2
    A <- LAMBDA + sigma_new;
    sigma2_new <- (A > lambda2) * (A - lambda2) + (A < -lambda2) * 
      (A + lambda2);
    # Update Lambda
    LAMBDA <- LAMBDA + (sigma_new - sigma2_new);
    
    err <- max( abs(sigma_new - sigma)/(abs(sigma) + 1), 
                abs(sigma2_new - sigma2)/(abs(sigma2) + 1) ); 
    k <- k + 1;
    sigma <- sigma_new;
    sigma2 <- sigma2_new;
  }
  
  if(k >= k_max) {
    cat("WARNING of cclasso_sub:\n", "\tMaximum Iteration:", k_max, 
        "&& Relative error:", err, "!\n");
  }
  
  return(sigma);
}

############################################################################################################################

# Demo Example from https://github.com/huayingfang/CCLasso

############################################################################################################################

n <- 192
p <- 747;
x <- matrix(rnorm(n * p), nrow = n); 
x.frac <- exp(x) / rowSums(exp((x)));
totCount <- round(runif(n = n,  min = 1000, max = 2000));
x.count <- x.frac * totCount;
# 2. run cclasso 
# using fraction
res_ccl_frac <- cclasso(x = x.frac, counts = F);
# using counts
res_ccl_count <- cclasso(x = x.count, counts = T);
# 3. run SparCC.count and SparCC.frac
res_spa_count <- SparCC.count(x = x.count);
res_spa_frac <- SparCC.frac(x = x.frac);
# 4. get the correlation matrix
{
  #cat("CCLasso using fraction data:\n");
  #print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(res_ccl_count$cor_w, 2));
  #cat("SparCC using fraction data:\n");
  #print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(res_spa_count$cor.w, 2));
}

############################################################################################################################

# British Genus Level Analysis 

############################################################################################################################

#importing excel of genus-level otu table that is filtered and rarefied
Genus3_Transposed <- read.delim("~/MoL_Analysis /1.6.2021/Analysis/CCLasso/0_Comparison_All_Genus_British_Modern_Filtered_Rarefy49000/Genus3_Transposed.txt", row.names=1)
Genus <- Genus3_Transposed

x <- data.matrix(Genus)

# 2. run cclasso ;
# using counts
Genus_res_ccl_count <- cclasso(x, counts = T);

# 3. run SparCC.count and SparCC.frac
Genus_res_spa_count <- SparCC.count(x);
#res_spa_frac <- SparCC.frac(x = x.frac);

# 4. get the correlation matrix
{
  #cat("CCLasso using fraction data:\n");
  #print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(Genus_res_ccl_count$cor_w, 2));
  #cat("SparCC using fraction data:\n");
  #print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(Genus_res_spa_count$cor.w, 2));
}

#renaming correlation matrix from CCLasso
Genus_Cor <- data.frame(Genus_res_ccl_count$cor_w)
names(Genus_Cor) <- names(Genus)
row.names(Genus_Cor) <- names(Genus)


#Exporting rescclcount cor2 matrix
#write.table(a, "cor_W")

#write.table(a2, "a2.txt")

############################################################################################################################

# British Species Level Analysis 

############################################################################################################################

#importing excel of species-level otu table that is filtered and rarefied
Species_Trans <- read.csv("~/MoL_Analysis /1.6.2021/Analysis/CCLasso/Species_Filtered_British_Modern_Filtered_Rarefy37000/Species_Trans.csv", row.names=1)
Species <- Species_Trans

x <- data.matrix(Species)

# 2. run cclasso ;
# using counts
Species_res_ccl_count <- cclasso(x, counts = T);

# 3. run SparCC.count and SparCC.frac
Species_res_spa_count <- SparCC.count(x);
#res_spa_frac <- SparCC.frac(x = x.frac);

# 4. get the correlation matrix
{
  #cat("CCLasso using fraction data:\n");
  #print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(Species_res_ccl_count$cor_w, 2));
  #cat("SparCC using fraction data:\n");
  #print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(Species_res_spa_count$cor.w, 2));
}

#renaming correlation matrix from CCLasso
Species_Cor <- data.frame(rSpecies_es_ccl_count$cor_w)
names(Species_Cor) <- names(Species_Trans)
row.names(Species_Cor) <- names(Species_Trans)


#Exporting rescclcount cor2 matrix
#write.table(a, "cor_W_genus")

############################################################################################################################

# British Genus Level Analysis: Strep 

############################################################################################################################

#importing excel of genus-level otu table that is filtered and rarefied
Genus_Strep_Trans <- read.delim("~/MoL_Analysis /1.6.2021/Analysis/CCLasso/0_Comparison_All_Genus_British_Modern_Filtered_Strep_Rarefy49000/Genus_Strep_Trans.txt", row.names=1)
Genus_strep <- Genus_Strep_Trans

x <- data.matrix(Genus_strep)

# 2. run cclasso ;
# using counts
Genus_strep_res_ccl_count <- cclasso(x, counts = T);

# 3. run SparCC.count and SparCC.frac
Genus_strep_res_spa_count <- SparCC.count(x);
#res_spa_frac <- SparCC.frac(x = x.frac);

# 4. get the correlation matrix
{
  #cat("CCLasso using fraction data:\n");
  #print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(Genus_res_ccl_count$cor_w, 2));
  #cat("SparCC using fraction data:\n");
  #print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(Genus_res_spa_count$cor.w, 2));
}

#renaming correlation matrix from CCLasso
Genus_Cor <- data.frame(Genus_res_ccl_count$cor_w)
names(Genus_Cor) <- names(Genus)
row.names(Genus_Cor) <- names(Genus)


#Exporting rescclcount cor2 matrix
#write.table(a, "cor_W")

#write.table(a2, "a2.txt")

############################################################################################################################

# British Genus Level Analysis: Meth

############################################################################################################################

#importing excel of genus-level otu table that is filtered and rarefied
Genus_Meth_Trans <- read.delim("~/MoL_Analysis /1.6.2021/Analysis/CCLasso/0_Comparison_All_Genus_British_Modern_Filtered_Meth_Rarefy49000/Genus_Meth_Trans.txt", row.names=1)
Genus_meth <- Genus_Meth_Trans

x <- data.matrix(Genus_meth)

# 2. run cclasso ;
# using counts
Genus_meth_res_ccl_count <- cclasso(x, counts = T);

# 3. run SparCC.count and SparCC.frac
Genus_meth_res_spa_count <- SparCC.count(x);
#res_spa_frac <- SparCC.frac(x = x.frac);

# 4. get the correlation matrix
{
  #cat("CCLasso using fraction data:\n");
  #print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(Genus_res_ccl_count$cor_w, 2));
  #cat("SparCC using fraction data:\n");
  #print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(Genus_res_spa_count$cor.w, 2));
}

#renaming correlation matrix from CCLasso
Genus_Cor <- data.frame(Genus_res_ccl_count$cor_w)
names(Genus_Cor) <- names(Genus)
row.names(Genus_Cor) <- names(Genus)


#Exporting rescclcount cor2 matrix
#write.table(a, "cor_W")

#write.table(a2, "a2.txt")

############################################################################################################################

# British Genus Level Analysis: Actino

############################################################################################################################

#importing excel of genus-level otu table that is filtered and rarefied
Genus_Actino_Trans <- read.delim("~/MoL_Analysis /1.6.2021/Analysis/CCLasso/0_Comparison_All_Genus_British_Modern_Filtered_Actino_Rarefy49000/Genus_Actino_Trans.txt", row.names=1)
Genus_Actino <- Genus_Actino_Trans

x <- data.matrix(Genus_Actino)

# 2. run cclasso ;
# using counts
Genus_Actino_res_ccl_count <- cclasso(x, counts = T);

# 3. run SparCC.count and SparCC.frac
Genus_Actino_res_spa_count <- SparCC.count(x);
#res_spa_frac <- SparCC.frac(x = x.frac);

# 4. get the correlation matrix
{
  #cat("CCLasso using fraction data:\n");
  #print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(Genus_res_ccl_count$cor_w, 2));
  #cat("SparCC using fraction data:\n");
  #print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(Genus_res_spa_count$cor.w, 2));
}
############################################################################################################################

# British Genus Level Analysis: Sample 70, min freq 50000
############################################################################################################################

#importing excel of genus-level otu table that is filtered and rarefied
genus_samp70_trans <- read.delim("~/MoL_Analysis /1.6.2021/Analysis/CCLasso/0_Comparison_All_Genus_British_Modern_Filtered_Rarefy49000_minsam70/genus_samp70_trans.txt", row.names=1)

x <- data.matrix(genus_samp70_trans)

# 2. run cclasso ;
# using counts
genus_samp70_res_ccl_count <- cclasso(x, counts = T);

# 3. run SparCC.count and SparCC.frac
genus_samp70_res_spa_count <- SparCC.count(x);
#res_spa_frac <- SparCC.frac(x = x.frac);

# 4. get the correlation matrix
{
  #cat("CCLasso using fraction data:\n");
  #print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(genus_samp70_res_ccl_count$cor_w, 2));
  #cat("SparCC using fraction data:\n");
  #print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(genus_samp70_res_spa_count$cor.w, 2));
}

