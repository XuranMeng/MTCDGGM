addpath(genpath(pwd));

n=178;
p = 73; % dimension of response
q = 120; % dimension of covariates


Z0=csvread('Z0.csv',1,0);
iU=csvread('iU.csv',1,0);

path='betaandM/'
hatbeta_mul=csvread(append(path,'hatbeta_mul.csv'));

% hatbeta_sepa=csvread(append(path,'hatbeta_sepa.csv'));

interM = [];
for j=1:(q+1)
interM=[interM Z0.*iU(:,j)]; %these are the interaction terms between Z0 and iU
end



adjustfactormullist=[]
betamuldebiaselist=[]
sigmamullist=[]
for j=1:p
    disp(j);
    betamul_j = hatbeta_mul((1:(p-1)*(q+1)) + (j-1)*(p-1)*(q+1), :);
    betamulzero=zeros(size(betamul_j));
    finalresultmulszero=zeros(size(betamul_j));
    adjustfactormulzero=zeros(size(betamul_j));
    sigma_jzero=zeros(size(betamul_j));
    Llistmul = find(betamul_j ~= 0);


    % betasepa_j = hatbeta_sepa((1:(p-1)*(q+1)) + (j-1)*(p-1)*(q+1), :);
    % betasepazero=zeros(size(betasepa_j));
    % finalresultsepaszero=zeros(size(betasepa_j));
    % adjustfactorsepazero=zeros(size(betasepa_j));
    
    % Llistsepa = find(betasepa_j ~= 0);


    Wj=interM;
    Wj(:,((0:q)*p+j)) = [];
    y1=Z0(:,j);


    if ~isempty(Llistmul)
        Mmul_j=csvread(append('betaandM/', 'Mmul_', string(j), '.csv'));
        MW_j = mexmatrixmul(Mmul_j', Wj');
        MW_jMW_jT = mexmatrixmul(MW_j, MW_j');  % This is MW_j * MW_j'
        a=betamul_j(Llistmul) + mexmatrixmul(MW_j, y1 - mexmatrixmul(Wj(:,Llistmul), betamul_j(Llistmul))) / n;
        Wjsub=Wj(:,Llistmul);
        sigma_j=norm(y1-Wjsub*(Wjsub'*Wjsub)^(-1)*Wjsub'*y1,2)/sqrt(n-numel(Llistmul))

        % hatsigma_jsquare=norm(y1 - mexmatrixmul(Wj(:,Llistmul), betamul_j(Llistmul)),2)^2/(n-numel(Llistmul));
        adjustfactor = sqrt(diag(MW_jMW_jT) / n);
        betatodo=betamulzero;
        adjustfactortodo=adjustfactormulzero;
        betatodo(Llistmul)=a;
        betamuldebiase_j=betatodo;
        adjustfactortodo(Llistmul)=adjustfactor;
        adjustfactormul_j=adjustfactortodo;
        sigmamul_j=sigma_jzero+sigma_j;
    else
        adjustfactormul_j=adjustfactormulzero;
        betamuldebiase_j=betamulzero;
        sigmamul_j=sigma_jzero;
    end

    % If Llistsepa is not empty, perform the sepa-related operations
    % if ~isempty(Llistsepa)
    %     Msepa_j=csvread(append('betaandM/', 'Msepaindex_', string(j), '.csv'))
    adjustfactormullist=[adjustfactormullist;adjustfactormul_j];
    betamuldebiaselist=[betamuldebiaselist;betamuldebiase_j];
    sigmamullist=[sigmamullist;sigmamul_j];
    


    % end
end
csvwrite('betaandM/adjustfactormullist.csv',adjustfactormullist)

csvwrite('betaandM/betamuldebiaselist.csv',betamuldebiaselist)

csvwrite('betaandM/sigmamullist.csv',sigmamullist)



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



