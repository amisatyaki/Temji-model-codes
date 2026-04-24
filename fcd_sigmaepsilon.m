function sigmaepsilon = fcd_sigmaepsilon(data,PHI,beta,N,d,aepsilon,bepsilon)
% data is NdT x 1 vector, containing all the data points
% PHI is NdT x K matrix containg phi_1, ... phi_N
% N is the number of participants
% d is the number of nodes
% beta is the coefficient of phi's
% aepsilon, bepsilon are hypriors of sigma_epsion_i
[NdT,~] = size(PHI);
T = NdT/(d*N);
sigmaepsilon = zeros(d,1);
error = data - PHI*beta;
reshapederror = reshape(error,d,(N*T));
for i=1:d
    scaleinv = bepsilon + sum(reshapederror(i,:).^2);
    sigmaepsilon(i) = 1/gamma(aepsilon+(N*T)/2,1/scaleinv);
end



