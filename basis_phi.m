function basis=basis_phi(wvector,x,bvector,K)
% wvector contains w_k, for k=1... K
% bvector contains b_k for k=1 ... K
% x is either a vector or a scaler at which basis will be evaluated
% K is the total number of basis
L=size(x,1);
basis = zeros(L,K);
for i=1:L
    for k=1:K
        basis(i,k) = sqrt(2/K)*cos(wvector(k)*x(i)+bvector(k));
    end
end