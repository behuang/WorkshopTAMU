nrec <- function(mpp, penalty=3) {
  if (!inherits(mpp, "mpprob")) stop("Must run mpprob first\n")
  output = list()  
  nfou=nrow(mpp$founders)
  mosaic = mpp$estfnd
  recevents = list()
  
  for (i in 1:length(mpp$map)) { 
  for (j in 1:nrow(mpp$finals)) {
      mat <- matrix(mpp$prob[[i]][j,], ncol=nfou, byrow=T)
      mosaic[[i]][j,] = optpath(mat, penalty)$ids
    }
    recevents[[i]] <- apply(mosaic[[i]], 1, function(x) return(sum(diff(x[!is.na(x)])!=0)))
  }
  ## Recombination events occurring in the genome
  output$mosaic <- mosaic
  ## Number of recombination events per chromosome
  output$nrec <- recevents
  ## Total number of recombinations per line
  output$totrec <- rowSums(do.call("cbind", output$nrec))
  output
}


#### C function for computing recombination events? 
library(inline)
library(RcppArmadillo)

optpath <- cxxfunction(signature(probs="numeric", penalty="numeric"), plugin="RcppArmadillo", body='
  arma::mat pr = Rcpp::as<arma::mat>(probs);
  int nr = pr.n_rows, nf=pr.n_cols;
  double mms = Rcpp::as<double>(penalty);
  arma::mat onevec = arma::ones<arma::mat>(nf, 1);
  arma::vec strains = arma::zeros<arma::vec>(nf);
  arma::mat pathscore = arma::zeros<arma::mat>(1, nf);
  arma::mat D = arma::ones<arma::mat>(nf, nf);
  arma::mat bestpath = arma::randu<arma::mat>(nr, nf);
  arma::mat y = arma::zeros<arma::mat>(1, nf);
  arma::mat X = arma::zeros<arma::mat>(nf, nf);
  arma::vec w = arma::zeros<arma::vec>(nf);
  arma::vec ids = arma::zeros<arma::vec>(nr);
  arma::vec steps = arma::zeros<arma::vec>(nr);
  arma::mat tmp = arma::zeros<arma::mat>(1, nf);
  int i, j;
  arma::uword m, k, step;

  for (i=0; i<nf; i++) strains(i) = i+1;
  D.diag() -= 1;  
  D = D*mms;

  for (i=0; i<nr; i++)
  {
    y = pr.row(i)+pathscore;
    X = (onevec * y)-D;
    for (j=0; j<nf; j++) {
      tmp = X.row(j);
      tmp.max(m, k);
      w(j) = k;
      pathscore(0,j) = X(j, w(j));
    }
    bestpath.row(i) = w.t();
  }

  pathscore.max(m, k);
  step = k;

  for (i=0; i<nr-1; i++) {
    steps(i) = bestpath(nr-2-i, step);
    step = steps(i);
  }
  steps(nr-1) = step;

  for (i=0; i<nr; i++) ids(i) = strains(steps(nr-i-1));

  return Rcpp::List::create(Rcpp::Named("ids")=ids);
' 
) 
