n=800;
p = 70; % dimension of response
q = 120; % dimension of covariates

nu=0.3;




mean(randn(200,1))


path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');



debiase1list=csvread(append(path,'/debiase1list.csv'));
debiase2list=csvread(append(path,'/debiase2list.csv'));
debiase3list=csvread(append(path,'/debiase3list.csv'));
debiase4list=csvread(append(path,'/debiase4list.csv'));

result1list=csvread(append(path,'/result1list.csv'));
result2list=csvread(append(path,'/result2list.csv'));
result3list=csvread(append(path,'/result3list.csv'));
result4list=csvread(append(path,'/result4list.csv'));




mean(result1list)
mean(result2list)
mean(result3list)
mean(result4list)


std(result1list)
std(result2list)
std(result3list)
std(result4list) 


sum(abs(result1list)<1.96)/200
sum(abs(result2list)<1.96)/200
sum(abs(result3list)<1.96)/200

sum(sum(result4list.^2,2)<5.991)/200









