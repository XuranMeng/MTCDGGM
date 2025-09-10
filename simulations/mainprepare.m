



mex mexfun/mexbwsolve.c
mex mexfun/mexfwsolve.c
mex mexfun/mexfz.c
mex mexfun/mexMatvec.c
mex mexfun/mexProjL2.c
mex mexfun/mexProxL2.c
mex mexfun/mextriang.c
mex mexfun/mexMatvecmultiple.c
mex mexfun/mexmatrixmul.c
mex mexfun/mexCondat.c
mex mexfun/mexsigma_update_Classic_Lasso_SSNAL.c
mex mexfun/mexscale.c
addpath(genpath(pwd));
savepath


p = 70; % dimension of response
q = 120; % dimension of covariates
nu=0.3;
Treptition=200;
dim=(p-1)*(q+1);


%Procedure 1
rng(1118)
generate(n,p,q,nu,Treptition)
