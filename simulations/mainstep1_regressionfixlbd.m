
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

p=70;
q=120;
nu=0.3;
dim=(p-1)*(q+1);
Treptition=200;


path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
    
beta=csvread(append(path,'beta.csv'));


parpool('local', 5) 


parfor reptition= 1:Treptition
    lbde=0.3
    alpha=2
    lbdg=lbde/sqrt(alpha)
    path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
    Z0=csvread(append(path,'Z0.csv'));
    iU=csvread(append(path,'iU.csv'));
    disp(reptition);
    [y_mul,hatbeta_mul]=mulGMMRegfixlbd(Z0,iU,lbde,lbdg);
    % find(hatbeta_mul~=0)'
    file_name_beta_mul=append('hatbeta_mullbde',string(lbde),'alpha',string(alpha) ,'.csv');
    file_name_y_mul=append('y_mullbde',string(lbde),'alpha',string(alpha) ,'.csv');
    csvwrite(append(path,file_name_beta_mul),hatbeta_mul);
    csvwrite(append(path,file_name_y_mul),y_mul);

    % 
end

delete(gcp('nocreate'));
% a=squeeze(hatbeta_sepa(1,:));
% a=a(:);
% size(a)
% a(120,:)
% find(a~=0)

% reptition=3
% path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
% hatbeta_mul=csvread(append(path,'hatbeta_mul.csv'));
% hatbeta_sepa=csvread(append(path,'hatbeta_sepa.csv'));

% hatbeta_mul(hatbeta_mul~=0)';
% hatbeta_mul(abs(hatbeta_mul)>0.1)';

% hatbeta_sepa(hatbeta_sepa~=0)';
% hatbeta_sepa(abs(hatbeta_sepa)>0.1)';


% A = [1 2 3; 4 5 6; 7 8 9] % Example matrix
% reshape(A,[],1) % convert matrix to column vector
