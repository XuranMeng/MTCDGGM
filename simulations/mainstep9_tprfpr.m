addpath(pwd);

p = 70; % dimension of response
q = 120; % dimension of covariates



lbde1=0.3;

nu=0.3;
path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
beta=csvread(append(path,'beta.csv'));





for reptition=1:200
    disp(reptition);
    path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
    
    Z0=csvread(append(path,'Z0.csv'));
    iU=csvread(append(path,'iU.csv'));

    interM = [];
    for h=1:(q+1)
        interM=[interM Z0.*iU(:,h)]; %these are the interaction terms between Z0 and iU
    end
    if n==800
        hatbeta=csvread(append(path,'/hatbeta_mullbde',string(lbde1),'alpha2.csv'));
    else
        hatbeta=csvread(append(path,'/hatbeta_mulcv.csv'));
    end
    


    % y_mulcv=csvread(append(path,'/y_mulcv.csv'));
    % y_mullbde1=csvread(append(path,'/y_mullbde',string(lbde1),'alpha2.csv'));
    % y_mullbde2=csvread(append(path,'/y_mullbde',string(lbde2),'alpha2.csv'));


    
    j=1;
    beta_j=beta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
    hatbeta_j=hatbeta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
    Llist=find(hatbeta_j~=0);
    if isempty(Llist)
        csvwrite(fullfile(path, 'compare_ours.csv'), []);
        continue;   
    end
    hatbeta_j_sub=hatbeta_j(Llist,:);


    y1=Z0(:,j);
    y=y1-mean(y1);

    W_j=interM;
    W_j(:,((0:q)*p+j)) = [];

    Mj=csvread(append(path,'/Mcompare_',string(j),'.csv'));
    M_j=Mj;



    MW_j = mexmatrixmul(M_j', W_j');
    MW_jMW_jT = mexmatrixmul(MW_j, MW_j');  % This is MW_j * MW_j'
    adjustfactor = sqrt(diag(MW_jMW_jT) / n);

    

    % residual_mulcv = y - mexmatrixmul(W_j, hatbeta_mulcv_j);
    % hatsd_cv=std(residual_mulcv)*sqrt(size(y,1)-1)/sqrt(size(y,1)-sum(hatbeta_mulcv_j~=0));
    residual = y - mexmatrixmul(W_j, hatbeta_j);
    hatsd=std(residual)*sqrt(size(y,1)-1)/sqrt(size(y,1)-sum(hatbeta_j~=0));
    debiased_hatbeta_j = hatbeta_j_sub + mexmatrixmul(MW_j, residual) / n;

    
    disp(Llist)
    sd_sel=sqrt(n)./adjustfactor./hatsd;
    z0=debiased_hatbeta_j.*sd_sel;
    pval=2 * (1 - normcdf(abs(z0)));

    ResultMat = [Llist, debiased_hatbeta_j, sd_sel,z0,pval];
    csvwrite(append(path,'/compare_ours.csv'),ResultMat);



  
end


