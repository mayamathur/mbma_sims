
# version that assumes only nonaffirmative results:
# https://github.com/mayamathur/phacking/blob/e335a054aedc448d9a1bc0f6e82300a0fb876583/Sherlock%20code/init_stan_model_MBMA.R

# # to test if this compiles in interactive session:
# # helper code
# path = "/home/groups/manishad/MBMA"
# setwd(path)
# source("init_stan_model_MBMA.R")


# 2022-5-3: SIMPLIFIED JEFFREYS PRIOR WRT TAU ----------------------


model.text <- "

functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){

		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
		matrix[2,2] fishinfo;
    real sigma;
		real LL;
		real UU;
		// will just be set to 1
		int n;

    real Si;
    real cz;
    real dnor;
    real pnor;
    real r;

		// this will be the TOTALS for all observations
		matrix[2,2] fishinfototal;
		fishinfototal[1,1] = 0;
  	fishinfototal[1,2] = 0;
  	fishinfototal[2,1] = 0;
  	fishinfototal[2,2] = 0;


		// build a Fisher info matrix for EACH observation
		for (i in 1:k) {

		  Si = sqrt( tau^2 + sei[i]^2 );
      cz = (sei[i] * tcrit[i] - mu) / Si;
      dnor = exp( normal_lpdf(cz | 0, 1 ) );
      pnor = exp( normal_lcdf(cz | 0, 1 ) );
      r = dnor/pnor;
      
      kmm = Si^(-2)*(cz*r + r^2 - 1);
      kms = tau*Si^(-3)*r*( cz^2 + cz*r + 1 );
      kss = ( tau^2 * Si^(-4) ) * ( cz^3*r + cz^2*r^2 + cz*r - 2 ); 

  		fishinfo[1,1] = -kmm;
      fishinfo[1,2] = -kms;
      fishinfo[2,1] = -kms;
      fishinfo[2,2] = -kss;

  		// add the new fisher info to the total one
  		fishinfototal = fishinfototal + fishinfo;
		}

		return sqrt(determinant(fishinfototal));
	}
}

data{
	int<lower=0> k;
  real sei[k];
  real tcrit[k];
  real affirm[k];
	real y[k];
}

parameters{
  real mu;
	real<lower=0> tau;
}


model{
  // this is to remove prior, as a sanity check:
  // target += 0;
  //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
	target += log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
	for(i in 1:k) {
      if ( affirm[i] == 0 ) {
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
      } else if ( affirm[i] == 1 ) {
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
      }
	}
}

// this chunk doesn't actually affect the model that's being fit to the data;
//  it's just re-calculating the prior, lkl, and post to return to user
// Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'

generated quantities{
  real log_lik = 0;
  real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
  real log_post;
  // this is just an intermediate quantity for log_lik
  // will be equal to UU or LL above, depending on affirm status
  real critScaled;

  // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
  //  to R functions for MAP, MLE, etc.
  real log_lik_sanity = 0;
  real log_prior_sanity = log( jeffreys_prior(2, 2, k, sei, tcrit, affirm) );

  for ( i in 1:k ){
    log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
    log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );

    critScaled = tcrit[i] * sei[i];

    // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
    // see 'Truncation with upper bounds in Stan' section
    // nonaffirm case:
    if ( y[i] <= critScaled ) {
    // from sanity checks in doParallel, I know this matches joint_nll_2
      log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
      log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );

    // affirm case:
    } else if ( y[i] > critScaled ) {
      log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
      log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
    }
  }
  log_post = log_prior + log_lik;
}
"

# necessary to prevent ReadRDS errors in which cores try to work with other cores' intermediate results
# https://groups.google.com/g/stan-users/c/8snqQTTfWVs?pli=1
options(mc.cores = parallel::detectCores())

# "isystem" arg is just a placeholder to avoid Stan's not understanding special characters
#  in getwd(), even though we don't actually use the dir at all
# note: removing the isystem arg does NOT fix the very sporadic "error reading from connection" on cluster
cat( paste("\n init_stan_model: about to call stan_model") )
stan.model <- stan_model(model_code = model.text,
                         isystem = "~/Desktop")


cat( paste("\n init_stan_model: done calling stan_model") )