n=400;
p = 70; % dimension of response
q = 120; % dimension of covariates
lbde1=0.3;
lbde2=0.6;
nu=0.3;







path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');

pre_bias_mullbde1_list=csvread(append(path,'pre_bias_mullbde1_list.csv'));
post_bias_mullbde1_list=csvread(append(path,'post_bias_mullbde1_list.csv'));
adjust_mullbde1_list=csvread(append(path,'adjust_mullbde1_list.csv'));


% pre_bias_sepalbde1_list=csvread(append(path,'pre_bias_sepalbde1_list.csv'));
% post_bias_sepalbde1_list=csvread(append(path,'post_bias_sepalbde1_list.csv'));
% adjust_sepalbde1_list=csvread(append(path,'adjust_sepalbde1_list.csv'));


% pre_bias_sepalbde2_list=csvread(append(path,'pre_bias_sepalbde2_list.csv'));
% post_bias_sepalbde2_list=csvread(append(path,'post_bias_sepalbde2_list.csv'));
% adjust_sepalbde2_list=csvread(append(path,'adjust_sepalbde2_list.csv'));

pre_bias_mullbde2_list=csvread(append(path,'pre_bias_mullbde2_list.csv'));
post_bias_mullbde2_list=csvread(append(path,'post_bias_mullbde2_list.csv'));
adjust_mullbde2_list=csvread(append(path,'adjust_mullbde2_list.csv'));



pre_bias_mulcv_list=csvread(append(path,'pre_bias_mulcv_list.csv'));
% pre_bias_sepacv_list=csvread(append(path,'pre_bias_sepacv_list.csv'));
post_bias_mulcv_list=csvread(append(path,'post_bias_mulcv_list.csv'));
% post_bias_sepacv_list=csvread(append(path,'post_bias_sepacv_list.csv'));
adjust_mulcv_list=csvread(append(path,'adjust_mulcv_list.csv'));
% adjust_sepacv_list=csvread(append(path,'adjust_sepacv_list.csv'));

AIL_list=csvread(append(path,'AIL_list.csv'))*1.96*2;

%%%%%%%%




clc;
disp('lbde1-part-mul:   ');
disp('Pre-Bias:   '), disp(mean((pre_bias_mullbde1_list))), disp('Pre-Bias-sd:  '),disp(std((pre_bias_mullbde1_list)));
disp('Post-Bias:  '), disp(mean((post_bias_mullbde1_list))),disp('Post-Bias-sd:  '),disp(std((post_bias_mullbde1_list)));
disp('Emp-SD:     '), disp(std(adjust_mullbde1_list)); 
disp('Cov-Prob:   '); disp(sum(abs(adjust_mullbde1_list)<1.96)/size(adjust_mullbde1_list,1)*100.0);
disp('Rej-0:   '); disp(sum((post_bias_mullbde1_list-0.3)./post_bias_mullbde1_list.*adjust_mullbde1_list<-1.96)/size(post_bias_mullbde1_list,1)*100.0);





clc;
disp('lbde2-part-mul:   ');
disp('Pre-Bias:   '), disp(mean((pre_bias_mullbde2_list))), disp('Pre-Bias-sd:  '),disp(std((pre_bias_mullbde2_list)));
disp('Post-Bias:  '), disp(mean((post_bias_mullbde2_list))),disp('Post-Bias-sd:  '),disp(std((post_bias_mullbde2_list)));
disp('Emp-SD:     '), disp(std(adjust_mullbde2_list)); 
disp('Cov-Prob:   '); disp(sum(abs(adjust_mullbde2_list)<1.96)/size(adjust_mullbde2_list,1)*100.0);
disp('Rej-0:   '); disp(sum((post_bias_mullbde2_list-0.3)./post_bias_mullbde2_list.*adjust_mullbde2_list<-1.96)/size(post_bias_mullbde2_list,1)*100.0);





clc;
disp('cv-part-mul:   ');
disp('Pre-Bias:   '), disp(mean((pre_bias_mulcv_list))), disp('Pre-Bias-sd:  '),disp(std((pre_bias_mulcv_list)));
disp('Post-Bias:  '), disp(mean((post_bias_mulcv_list))),disp('Post-Bias-sd:  '),disp(std((post_bias_mulcv_list)));
disp('Emp-SD:     '), disp(std(adjust_mulcv_list)); 
disp('Cov-Prob:   '); disp(sum(abs(adjust_mulcv_list)<1.96)/size(adjust_mulcv_list,1)*100.0);
disp('Rej-0:   '); disp(sum((post_bias_mulcv_list-0.3)./post_bias_mulcv_list.*adjust_mulcv_list<-1.96)/size(post_bias_mulcv_list,1)*100.0);



disp('AIL:   '); disp(mean(AIL_list,1));
disp('AILstd:   '); disp(std(AIL_list,1));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%Orcale Part
n=800;
p = 70; % dimension of response
q = 120; % dimension of covariates


nu=0.3


path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');

oracle_bias=csvread(append(path,'oracle_bias.csv'));


clc
disp('oracle_Property:   ');
disp('Bias:   '), disp(mean((oracle_bias))), disp('Pre-Bias-sd:  '),disp(std((oracle_bias)));
disp('Emp-SD:     '), disp(std(oracle_bias)*sqrt(n)); 
disp('Cov-Prob:   '); disp(sum(abs(oracle_bias*sqrt(n))<1.96)/size(oracle_bias,1)*100.0);
disp('Rej-0:   '); disp(sum(((oracle_bias-0.3)*sqrt(n))<-1.96)/size(oracle_bias,1)*100.0);
% disp('AIL:   '); disp(sum(((oracle_bias-0.3)*sqrt(n))<-1.96)/size(oracle_bias,1)*100.0);

2*1.96/sqrt(n)













