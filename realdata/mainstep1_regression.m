mex mexfun/mexbwsolve.c
mex mexfun/mexfwsolve.c
mex mexfun/mexfz.c
mex mexfun/mexMatvec.c
mex mexfun/mexProjL2.c
mex mexfun/mexProxL2.c
mex mexfun/mextriang.c
mex mexfun/mexMatvecmultiple.c
mex mexfun/mexmatrixmul.c

savepath
addpath(genpath(pwd));

pools=5
p=73

q=120
n=178
Z0=csvread('Z0emma.csv',1,0);
iU=csvread('iUemma.csv',1,0);
size(Z0)

[y_mul,hatbeta_mul]=mulGMMReg(Z0,iU,pools)
% find(hatbeta_mul~=0)'
Llist=find(hatbeta_mul~=0);
numel(Llist)

[jlist,jplist,hlist]=betaindex(p,q,Llist);
unames=readtable('unames.csv')
unamelist=unique(hlist)
unames.x(unamelist(2:10)-2)


mkdir('betaandM/')
% hatbeta_mul(abs(hatbeta_mul)>0.1)
csvwrite('betaandM/hatbeta_mul.csv',hatbeta_mul);
csvwrite('betaandM/y_mul.csv',y_mul);


