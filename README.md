In PG.m, particle Gibbs sampling is done for updating x_{n.t}. 
Particle Gibbs is supposed to be better in updating latent variables which varies with time. 
This has to be updated for each subject individually. 

%----------------------


In basis_phi.m the phi matrix has been calculated given a x.
In two other files, the sigma2epsilon, beta have been updated using full conditional densities. 
