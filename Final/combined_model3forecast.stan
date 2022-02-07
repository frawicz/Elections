data {
  // Data for dynamic component of the model
  int<lower=0> nParties; // Number of Parties
  int<lower=0> nPolls; // Number of published Polls
  int<lower=0> nPeriods; // Number of days the model should account for
  int<lower=0> y[nPolls,nParties]; // The poll results
  int<lower=0> nInst;
  int<lower=0> nforecast; 
  int<lower=0, upper = nInst> iid[nPolls];
  vector<lower=0>[nParties] a_pred ;
  int<lower=0, upper = nPeriods> date[nPolls]; // Date (as days from 1 to nPeriods) when the Poll was published
  int prior_house[nInst, nParties]; // prior for house effects;
}
parameters {
  simplex[nParties] vE;
  cholesky_factor_corr[nParties-1] S_cor; 
  cholesky_factor_corr[nParties-1] S_shock_cor; 
  vector<lower=0>[nParties-1] sigma_evo;
  vector<lower=0>[nParties-1] sigma_shock;
  matrix[nforecast,nParties-1] alphastarforecast;
  matrix[nPeriods-1,nParties-1] alphastar;
  matrix[nInst,nParties] house_effect_raw;
}
transformed parameters {
  matrix[nPeriods,nParties] alphastar_prior;
  matrix[nPeriods,nParties] testpar;
  matrix[nforecast,nParties] alphastarforecast_prior;
  matrix[nPeriods,nParties] beta;
  matrix[nPeriods,nParties] ea;
  matrix[nInst,nParties] house_effect;
  matrix[nforecast,nParties] ef;
  matrix[nforecast,nParties] forecast;
  


  alphastar_prior[1:(nPeriods-1),1:(nParties-1)] = alphastar;
  alphastarforecast_prior[1:(nforecast),1:(nParties-1)] = alphastarforecast;
  //print("alphastarforecast_prior");
  //print(alphastarforecast_prior);

  house_effect[1:nInst, 1:nParties] = house_effect_raw;
 


  // Transform structural forecast to log-ratio
  for (j in 1:(nParties)){
    alphastar_prior[nPeriods,j] = log(vE[j]/vE[nParties]);
    testpar[nPeriods,j] = alphastar_prior[nPeriods,j]; 	
}

  for (i in 1:(nPeriods-1)) {
    alphastar_prior[i,nParties] = 0;
    testpar[i,nParties] = alphastar_prior[i,nParties];   
}

  for(k in 1:nPolls) {
    testpar[date[k],] = alphastar_prior[date[k],];
    alphastar_prior[date[k],] = (alphastar_prior[date[k],]  + house_effect[iid[k],]);
}
for(fore in 1:nforecast){
  alphastarforecast_prior[fore,nParties] = 0;  
}
  // Transform values from log-ratio space back to vote shares
  for (i in 1:(nPeriods)) 
    for (j in 1:(nParties)) 
      ea[i,j] = exp(alphastar_prior[i,j]);

  for (i in 1:(nPeriods)) 
    for (j in 1:(nParties)) 
      beta[i,j] = ea[i,j]/sum(ea[i,]);   




for(fore in 1:nforecast){
  for (j in 1:(nParties))
    ef[fore,j] = exp(alphastarforecast_prior[fore,j]);
	//print("ef");
	//print(ef);  

  for (j in 1:(nParties))   
    forecast[fore,j] = ef[fore,j]/sum(ef[fore,1:nParties]);
	//print("forecast");
	//print(forecast);    
}
 

              
}
model {
  
  int pos; // Auxilliary counter
  pos = 1; // Pre-set counter to 1
	

 //for (i in 1:(nParties)){
 // vE[i] ~ dirichlet(a_pred[i]);
//	}
  vE ~ dirichlet(a_pred);


  // Dynamic component of the model

  // Evolution variance prior
  sigma_evo ~ normal(0, 0.02);

  // Forecast variance shock prior
  sigma_shock ~ normal(0, 0.02);

  // Evolution covariance prior
  S_cor ~ lkj_corr_cholesky(50);

  // Forecast covariance prior
  S_shock_cor ~ lkj_corr_cholesky(100);



  // Backward Random Walk of latent support (in log-ratio space)
 // for (i in (nPeriods-1):1){
//	print("Alphastar Prior: ",  alphastar_prior[(i+1),1:(nParties-1)]);
 //   alphastar[i,1:(nParties-1)] ~ multi_normal_cholesky(alphastar_prior[(i+1),1:(nParties-1)], diag_pre_multiply(sigma_evo, S_cor)); 
//	print("Alphastar: ",  alphastar[i,1:(nParties-1)]);
 //   }



  // Backward Random Walk of latent support (in log-ratio space)
  for (i in 2:(nPeriods-1)){
    alphastar[i,1:(nParties-1)] ~ multi_normal_cholesky(alphastar_prior[(i-1),1:(nParties-1)], diag_pre_multiply(sigma_evo, S_cor)); 
    }




  // Add forecast variance schock 
  alphastarforecast[1,1:(nParties-1)] ~ multi_normal_cholesky(alphastar_prior[nPeriods, 1:(nParties-1)], diag_pre_multiply(sigma_shock, S_shock_cor));

  for(fore in 2:nforecast){
   alphastarforecast[fore,1:(nParties-1)] ~ multi_normal_cholesky(alphastarforecast[fore-1, 1:(nParties-1)], diag_pre_multiply(sigma_shock, S_shock_cor));
   //print(alphastarforecast);
}

  // Model Polls based on latent state (in vote share space)
  for (k in 1:nPolls){
    y[k,1:nParties] ~ multinomial(beta[date[k],]'); 
    }
    
    
  // Estimate house effects
  for (j in 1:(nParties)) 
        for (c in 1:(nInst))
    //      house_effect_raw[c, j] ~ normal(prior_house[c, j]/1000.0, 0.02); // Prior 1 percent point sd
		house_effect_raw[c, j] ~ normal(0, 0.02); // Prior 1 percent point sd
  
}
