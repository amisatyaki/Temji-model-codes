function beta = fcd_beta(data,sigma_epsilon,PHI,N,tau)
% data is NdT x 1 vector, containing all the data points
% sigma_epsilon is a vector containing sigma^2_epsilon1, ...,
% sigma^2_epsilond
% PHI is NdT x K matrix containg phi_1, ... phi_N
% N is the number of participants
% tau is the prior variance of components of beta
[NdT,K] = size(PHI);
d = size(epsilon,1);
T = NdT/(d*N);
Depsiloninv = kron(eye(T),diag(1./sigma_epsilon));
Sigmaepsiloninv = kron(eye(N),Depsiloninv);
Sigmabetainv = PHI'*Sigmaepsiloninv*PHI + tau^2*eye(K);
Sigmabeta = inv(Sigmabetainv);
mubeta = Sigmabeta*PHI'*Sigmaepsiloninv*data;
beta = mvnrnd(mubeta,Sigmabeta)';