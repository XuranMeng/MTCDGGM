addpath(pwd);


p = 70; % dimension of response
q = 120; % dimension of covariates






lbde1=0.3
lbde2=0.6









nu=0.3;
path=append("datarepounknown/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
beta=csvread(append(path,'beta.csv'));





pre_bias_mullbde1_list=[];
post_bias_mullbde1_list=[];
adjust_mullbde1_list=[];



pre_bias_mullbde2_list=[];
post_bias_mullbde2_list=[];
adjust_mullbde2_list=[];


AIL_list=[];

pre_bias_mulcv_list=[];
post_bias_mulcv_list=[];
% adjust_mulcv_list=[];

for reptition=1:200
    disp(reptition);
    path=append("datarepounknown/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
    
    X0=csvread(append(path,'Z0.csv'));
    iU=csvread(append(path,'iU.csv'));
    hatbGamma=csvread(append(path,'hatbGamma.csv'),0,0);
    Z0=X0-iU(:,2:end)*hatbGamma;

    interM = [];
    for h=1:(q+1)
        interM=[interM Z0.*iU(:,h)]; %these are the interaction terms between Z0 and iU
    end

    % hatbeta_mulcv=csvread(append(path,'/hatbeta_mulcv.csv'));
    hatbeta_mullbde1=csvread(append(path,'/hatbeta_mullbde',string(lbde1),'alpha2.csv'));
    hatbeta_mullbde2=csvread(append(path,'/hatbeta_mullbde',string(lbde2),'alpha2.csv'));


    % y_mulcv=csvread(append(path,'/y_mulcv.csv'));
    % y_mullbde1=csvread(append(path,'/y_mullbde',string(lbde1),'alpha2.csv'));
    % y_mullbde2=csvread(append(path,'/y_mullbde',string(lbde2),'alpha2.csv'));


    





    collect_pre_bias_mullbde1_list=[];
    collect_post_bias_mullbde1_list=[];
    collect_adjust_mullbde1_list=[];
    collect_pre_bias_mullbde2_list=[];
    collect_post_bias_mullbde2_list=[];
    collect_adjust_mullbde2_list=[];
    collect_pre_bias_mulcv_list=[];
    collect_post_bias_mulcv_list=[];
    collect_adjust_mulcv_list=[];
    collect_interval_length=[];
    for j=1:2
        beta_j=beta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        Llist=find(beta_j~=0);
        beta_j_sub=beta_j(Llist,:);
        % hatbeta_mulcv_j=hatbeta_mulcv((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        hatbeta_mullbde1_j=hatbeta_mullbde1((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        hatbeta_mullbde2_j=hatbeta_mullbde2((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        % y_mulcv_j=y_mulcv((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        % y_mullbde1_j=y_mullbde1((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        % y_mullbde2_j=y_mullbde2((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);


        % hatbeta_mulcv_j_sub=hatbeta_mulcv_j(Llist,:);
        hatbeta_mullbde1_j_sub=hatbeta_mullbde1_j(Llist,:);
        hatbeta_mullbde2_j_sub=hatbeta_mullbde2_j(Llist,:);

        y1=Z0(:,j);
        y=y1-mean(y1);

        W_j=interM;
        W_j(:,((0:q)*p+j)) = [];

        Mj=csvread(append(path,'/M_',string(j),'.csv'));
        M_j=Mj;



        MW_j = mexmatrixmul(M_j', W_j');
        MW_jMW_jT = mexmatrixmul(MW_j, MW_j');  % This is MW_j * MW_j'
        adjustfactor = sqrt(diag(MW_jMW_jT) / n);

       

        % residual_mulcv = y - mexmatrixmul(W_j, hatbeta_mulcv_j);
        % hatsd_cv=std(residual_mulcv)*sqrt(size(y,1)-1)/sqrt(size(y,1)-sum(hatbeta_mulcv_j~=0));
        residual_mullbde1 = y - mexmatrixmul(W_j, hatbeta_mullbde1_j);
        hatsd_lbd1=std(residual_mullbde1)*sqrt(size(y,1)-1)/sqrt(size(y,1)-sum(hatbeta_mullbde1_j~=0));
        residual_mullbde2 = y - mexmatrixmul(W_j, hatbeta_mullbde2_j);
        hatsd_lbd2=std(residual_mullbde2)*sqrt(size(y,1)-1)/sqrt(size(y,1)-sum(hatbeta_mullbde2_j~=0));


        % Compute debiased beta for both cases using precomputed residuals
        % debiased_hatbeta_mulcv_j = hatbeta_mulcv_j_sub + mexmatrixmul(MW_j, residual_mulcv) / n;
        debiased_hatbeta_mullbde1_j = hatbeta_mullbde1_j_sub + mexmatrixmul(MW_j, residual_mullbde1) / n;
        debiased_hatbeta_mullbde2_j = hatbeta_mullbde2_j_sub + mexmatrixmul(MW_j, residual_mullbde2) / n;



        

        
        collect_pre_bias_mullbde1_list=[collect_pre_bias_mullbde1_list,hatbeta_mullbde1_j_sub'+0.3];
        collect_post_bias_mullbde1_list=[collect_post_bias_mullbde1_list,debiased_hatbeta_mullbde1_j'+0.3];
        collect_adjust_mullbde1_list=[collect_adjust_mullbde1_list,(debiased_hatbeta_mullbde1_j'+0.3)/hatsd_lbd1*sqrt(n)./adjustfactor'];

        collect_pre_bias_mullbde2_list=[collect_pre_bias_mullbde2_list,hatbeta_mullbde2_j_sub'+0.3];
        collect_post_bias_mullbde2_list=[collect_post_bias_mullbde2_list,debiased_hatbeta_mullbde2_j'+0.3];
        collect_adjust_mullbde2_list=[collect_adjust_mullbde2_list,(debiased_hatbeta_mullbde2_j'+0.3)/hatsd_lbd2*sqrt(n)./adjustfactor'];
        collect_interval_length=[collect_interval_length,adjustfactor'./sqrt(n)];


        % collect_pre_bias_mulcv_list=[collect_pre_bias_mulcv_list,hatbeta_mulcv_j_sub'+0.3];
        % collect_post_bias_mulcv_list=[collect_post_bias_mulcv_list,debiased_hatbeta_mulcv_j'+0.3];
        % collect_adjust_mulcv_list=[collect_adjust_mulcv_list,(debiased_hatbeta_mulcv_j'+0.3)/hatsd_cv*sqrt(n)./adjustfactor'];
 


    end
    pre_bias_mullbde1_list=[pre_bias_mullbde1_list;collect_pre_bias_mullbde1_list];
    post_bias_mullbde1_list=[post_bias_mullbde1_list;collect_post_bias_mullbde1_list];
    adjust_mullbde1_list=[adjust_mullbde1_list;collect_adjust_mullbde1_list];


    pre_bias_mullbde2_list=[pre_bias_mullbde2_list;collect_pre_bias_mullbde2_list];
    post_bias_mullbde2_list=[post_bias_mullbde2_list;collect_post_bias_mullbde2_list];
    adjust_mullbde2_list=[adjust_mullbde2_list;collect_adjust_mullbde2_list];
    AIL_list=[AIL_list;collect_interval_length];
    



    % pre_bias_mulcv_list=[pre_bias_mulcv_list;collect_pre_bias_mulcv_list];
    % post_bias_mulcv_list=[post_bias_mulcv_list;collect_post_bias_mulcv_list];
    % adjust_mulcv_list=[adjust_mulcv_list;collect_adjust_mulcv_list];



end


path=append("datarepounknown/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
csvwrite(append(path,'pre_bias_mullbde1_list.csv'),pre_bias_mullbde1_list);
csvwrite(append(path,'post_bias_mullbde1_list.csv'),post_bias_mullbde1_list);
csvwrite(append(path,'adjust_mullbde1_list.csv'),adjust_mullbde1_list);



csvwrite(append(path,'pre_bias_mullbde2_list.csv'),pre_bias_mullbde2_list);
csvwrite(append(path,'post_bias_mullbde2_list.csv'),post_bias_mullbde2_list);
csvwrite(append(path,'adjust_mullbde2_list.csv'),adjust_mullbde2_list);


csvwrite(append(path,'AIL_list.csv'),AIL_list);
% csvwrite(append(path,'pre_bias_mulcv_list.csv'),pre_bias_mulcv_list);
% csvwrite(append(path,'post_bias_mulcv_list.csv'),post_bias_mulcv_list);
% csvwrite(append(path,'adjust_mulcv_list.csv'),adjust_mulcv_list);




