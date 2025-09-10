




n=800;
p = 70; % dimension of response
q = 120; % dimension of covariates







lbde=0.6



nu=0.3;

path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
beta=csvread(append(path,'beta.csv'));



debiase1list=[];
debiase2list=[];
debiase3list=[];
debiase4list=[];


result1list=[];
result2list=[];
result3list=[];
result4list=[];

for reptition=1:200
    disp(reptition);
    path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
    
    Z0=csvread(append(path,'Z0.csv'));
    iU=csvread(append(path,'iU.csv'));

    interM = [];
    for h=1:(q+1)
        interM=[interM Z0.*iU(:,h)]; %these are the interaction terms between Z0 and iU
    end


    hatbeta=csvread(append(path,'/hatbeta_mullbde',string(lbde),'alpha2.csv'));



    for j=1
        beta_j=beta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        Llist=find(beta_j~=0);
        Llist=[1;2;Llist];
        beta_j_sub=beta_j(Llist,:);
        

        hatbeta_j=hatbeta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        hatbeta_j_sub=hatbeta_j(Llist,:);

        y1=Z0(:,j);
        y=y1-mean(y1);


        W_j=interM;
        W_j(:,((0:q)*p+j)) = [];
        M_j=csvread(append(path,'/Mnew_',string(j),'.csv'));
        % M_j=Mj(:,Llist);
        % numel(find(Mj~=0))
        % find(M_j~=0)
        % M_j=csvread(append(path,'/Mmulgroup_',string(j),'.csv'));
        % M_j=csvread(append(path,'/Mrand_',string(j),'.csv'));
        % M_j=csvread(append(path,'/Msinglegroup_',string(j),'.csv'));

        MW_j = mexmatrixmul(M_j', W_j');
        MW_jMW_jT = mexmatrixmul(MW_j, MW_j');  % This is MW_j * MW_j'
        adjustfactor = sqrt(diag(MW_jMW_jT) / n);

        

        residual = y - mexmatrixmul(W_j, hatbeta_j);
        % Compute debiased beta for both cases using precomputed residuals
        debiased_hatbeta__j = hatbeta_j_sub + mexmatrixmul(MW_j, residual) / n;

        sigma1square=mexmatrixmul(mexmatrixmul([0,0,1,-1],MW_jMW_jT),[0;0;1;-1;])/n;
        sigma2square=mexmatrixmul(mexmatrixmul([1,0,-1,0],MW_jMW_jT),[1;0;-1;0;])/n;
        sigma3square=mexmatrixmul(mexmatrixmul([1,2,0,0],MW_jMW_jT),[1;2;0;0;])/n;
        sigma4square=mexmatrixmul(mexmatrixmul([2,0,-1,0;0,1,1,0],MW_jMW_jT),[2,0;0,1;-1,1;0,0;])/n;

        debiase1=(mexmatrixmul([0,0,1,-1],debiased_hatbeta__j)-0);
        debiase2=(mexmatrixmul([1,0,-1,0],debiased_hatbeta__j)-0.3);
        debiase3=(mexmatrixmul([1,2,0,0],debiased_hatbeta__j)-0.0);
        debiase4=(sigma4square^(-0.5)*(mexmatrixmul([2,0,-1,0;0,1,1,0],debiased_hatbeta__j)-[0.3;-0.3]));


        result1=(mexmatrixmul([0,0,1,-1],debiased_hatbeta__j)-0)*sqrt(n)/sqrt(sigma1square);
        result2=(mexmatrixmul([1,0,-1,0],debiased_hatbeta__j)-0.3)*sqrt(n)/sqrt(sigma2square);
        result3=(mexmatrixmul([1,2,0,0],debiased_hatbeta__j)-0.0)*sqrt(n)/sqrt(sigma3square);
        result4=(sigma4square^(-0.5)*(mexmatrixmul([2,0,-1,0;0,1,1,0],debiased_hatbeta__j)-[0.3;-0.3]))*sqrt(n);



    end

    debiase1list=[debiase1list;debiase1];
    debiase2list=[debiase2list;debiase2];
    debiase3list=[debiase3list;debiase3];
    debiase4list=[debiase4list;debiase4'];


    result1list=[result1list;result1];
    result2list=[result2list;result2];
    result3list=[result3list;result3];
    result4list=[result4list;result4'];
end

path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');


 csvwrite(append(path,'/debiase1list.csv'),debiase1list);
 csvwrite(append(path,'/debiase2list.csv'),debiase2list);
 csvwrite(append(path,'/debiase3list.csv'),debiase3list);
 csvwrite(append(path,'/debiase4list.csv'),debiase4list);



 csvwrite(append(path,'/result1list.csv'),result1list);
 csvwrite(append(path,'/result2list.csv'),result2list);
 csvwrite(append(path,'/result3list.csv'),result3list);
 csvwrite(append(path,'/result4list.csv'),result4list);
% size(pre_bias_sepalbde1_list)



