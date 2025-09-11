mex mexfun/mexbwsolve.c
mex mexfun/mexfwsolve.c
mex mexfun/mexfz.c
mex mexfun/mexMatvec.c
mex mexfun/mexProjL2.c
mex mexfun/mexProxL2.c
mex mexfun/mextriang.c
mex mexfun/mexMatvecmultiple.c
mex mexfun/mexmatrixmul.c
addpath(genpath(pwd));
savepath


Treptition=200;

startrep = 1;

p = 70;
q = 120;
% if n == 200
%     startrep = 1;
% else
%     startrep = 1;
% end

nu=0.3;
dim=(p-1)*(q+1);
if n==400
    pools=4;
elseif n==200
    pools = 7;
else
    pools=8;
end


for reptition=startrep:Treptition
    disp([n p q reptition]);
    path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
    Z0=csvread(append(path,'Z0.csv'));
    iU=csvread(append(path,'iU.csv'));
    disp(reptition);
    [y_mul,hatbeta_mul]=mulGMMReg(Z0,iU,pools);
    csvwrite(append(path,'hatbeta_mulcv.csv'),hatbeta_mul);
    csvwrite(append(path,'y_mulcv.csv'),y_mul);
    % 
end
% end
