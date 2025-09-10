%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(pwd));

n=178;
p = 73; % dimension of response
q = 120; % dimension of covariates

path='betaandM/'
hatbeta_mul=csvread(append(path,'hatbeta_mul.csv'));
betamuldebiaselist=csvread(append(path,'betamuldebiaselist.csv'));
sigmamullist=csvread(append(path,'sigmamullist.csv'));
adjustfactormullist=csvread(append(path,'adjustfactormullist.csv'));

Llist=find(hatbeta_mul~=0);
Llistpostive=find(hatbeta_mul>0.000001);
Llistnegative=find(hatbeta_mul<-0.000001);



% [jlist,jplist,hlist]=betaindex(p,q,Llist);
% checklist1=jplist(find(jlist==1))
% checklist2=hlist(find(jlist==1))

% (p-1).*checklist2+checklist1-1

% Lcompare=find(betamuldebiaselist~=0);
% numel(Llist)

% sum(abs(Llist-Lcompare))




% hatbeta_sepa=csvread(append(path,'hatbeta_sepa.csv'));

% find(adjustfactormullist~=0)-find(hatbeta_mul~=0)

% csvwrite('betaandM/adjustfactormullist.csv',adjustfactormullist)

% csvwrite('betaandM/betamuldebiaselist.csv',betamuldebiaselist)


% adjustfactormullist=csvread('betaandM/adjustfactormullist.csv');

% betamuldebiaselist=csvread('betaandM/betamuldebiaselist.csv');
% betamuldebiaselist(betamuldebiaselist~=0)'
% hatbeta_mul(hatbeta_mul~=0)'


% sigmamullist=csvread('betaandM/sigmamullist.csv');




size(adjustfactormullist)


finalresults=betamuldebiaselist(Llist)./adjustfactormullist(Llist)*sqrt(n);
% finalresults=betamuldebiaselist(Llist)./adjustfactormullist(Llist)*sqrt(n)./sigmamullist(Llist);

indexorigin=Llist;
indexoriginpositive=find(hatbeta_mul>0.000001);
indexoriginnegative=find(hatbeta_mul<-0.000001);
numel(indexorigin)


index005=Llist(find(abs(finalresults)>norminv(1 - 0.05/2)));
index005positive=Llist(find(finalresults>norminv(1 - 0.05/2)));
index005negative=Llist(find(finalresults<-norminv(1 - 0.05/2)));
numel(index005)


index001=Llist(find(abs(finalresults)>norminv(1 - 0.01/2)));
index001positive=Llist(find((finalresults)>norminv(1 - 0.01/2)));
index001negative=Llist(find((finalresults)<-norminv(1 - 0.01/2)));



index0001=Llist(find(abs(finalresults)>norminv(1 - 0.001/2)));
index0001positive=Llist(find((finalresults)>norminv(1 - 0.001/2)));
index0001negative=Llist(find((finalresults)<-norminv(1 - 0.001/2)));

numel(index0001)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[jlist,jplist,hlist]=betaindex(p,q,Llist);
unamelist=unique(hlist)
unames.x(unamelist(2))






