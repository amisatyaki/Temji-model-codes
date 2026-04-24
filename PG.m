function [xn]=PG(data,theta,xn_0,x0_n,T,d,bvector,wvector,K)
%-----------------------------------------------
% This is particle Gibbs algorithm following Andrieu et al. (2010)
% data contains all the observations for nth subject, dT x 1 vector
% theta contains all the parameters other than latents, it is a list
% theta = {beta, sigma2epsilon, An, mn_1,... mn_d, sigma2eta, W1, ..., WT}
% xn_0 is the initial value of xn vector containing dT x 1 elements, that
% is xn_0 = (x_n.1 (0), ..., x_n.T(0))'
% x0_n is d x 1 vector of latent variable at time point 0
% T is the total number of time points
% d is the number of nodes
%-----------------------------------------
N = 10* T; % no of particles, it can be changed
beta = theta{1,1};
sigma2epsilon = theta{1,2};
D_epsiloninv = diag(1./sigma2epsilon);
An = theta{1,3};
mn = zeros(L,d);
xn_0 = reshape(xn_0,d,T);
data = reshape(data,d,T);
for i=(3+1):(3+d)
    mn(:,i-3) = theta{1,i};
end
sigma2eta = theta{1,3+d+1};
D_eta = diag(sigma2eta);
mundot1 = zeros(d,1);
xntemp = zeros(T,d,N);
w1=zeros(N,1);
a=zeros(T-1,N);
for t=1:T
    W = theta{1,3+d+1+t};
    if (t==1)
        for i=1:d
            u=zeros(d,1);
            for j=1:d
                u = u+ (An(i,j)*x0_n(j))*mn(:,j);
            end
            mundot1(i) = mn(:,i)'*W*u;
        end
        for r=1:(N-1)
            xntemp(t,:,r)=mvnrnd(mundot1,D_eta)';
            w1(r) = (-0.5)*(data(:,t) - basis_phi(wvector,xntemp(t,:,r),bvector,K)*beta)'*...
                D_epsiloninv*(data(:,t) - basis_phi(wvector,xntemp(t,:,r),bvector,K)*beta);
            w1(r) = exp(w1(r));
        end
        xntemp(t,:,N) = xn_0(:,t);
        w1(N) = (-0.5)*(data(:,t) - basis_phi(wvector,xntemp(:,N),bvector,K)*beta)'*...
                D_epsiloninv*(data(:,t) - basis_phi(wvector,xntemp(:,N),bvector,K)*beta);
        w1(N) = exp(w1(N));    
    end
    w1 = w1./(sum(w1));
    if(t>1)
        for r=1:(N-1)
            a(t-1,r)=randsample(N,1,true,w1);
            x0_n = xntemp(t-1,:,a);
            for i=1:d
                u=zeros(d,1);
                for j=1:d
                    u = u+ (An(i,j)*x0_n(j))*mn(:,j);
                end
                mundot1(i) = mn(:,i)'*W*u;
            end
            xntemp(t,:,r) = mvnrnd(mundot1,D_eta)';
            w1(r) = (-0.5)*(data(:,t) - basis_phi(wvector,xntemp(t,:,r),bvector,K)*beta)'*...
                D_epsiloninv*(data(:,t) - basis_phi(wvector,xntemp(t,:,r),bvector,K)*beta);
            w1(r) = exp(w1(r));
        end
        a(t-1,N) = N;
        xntemp(t,:,N) = xn_0(:,t);
        w1(N) = (-0.5)*(data(:,t) - basis_phi(wvector,xntemp(:,N),bvector,K)*beta)'*...
                D_epsiloninv*(data(:,t) - basis_phi(wvector,xntemp(:,N),bvector,K)*beta);
        w1(N) = exp(w1(N));
    end
end
b = randsample(N,1,true,w1);
xn=zeros(d,T);
for t=1:T
    b1=b;
    xn(:,T-t+1)=xntemp(T-t+1,:,b1);
    b=a(T-t,b1);
end
