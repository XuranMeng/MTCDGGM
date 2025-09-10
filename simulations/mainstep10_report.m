addpath(pwd);

n=800;

p = 70; % dimension of response
q = 120; % dimension of covariates
dim=(p-1)*(q+1);
lbde1=0.3;

nu=0.3;
path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
beta=csvread(append(path,'beta.csv'));







FPR_ours_list=[];
TPR_ours_list=[];
FPR_deco_list=[];
TPR_deco_list=[];

threshold=0.001;
for reptition = 1:200
    disp(reptition)
    path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
        

    filename_ours=append(path,'compare_ours.csv');
    fileInfo = dir(filename_ours);
    if fileInfo.bytes > 0
        compare_ours = csvread(filename_ours);
    else
        % Handle empty case
        compare_ours = [];
    end

    filename_deco=append(path,'compare_deco.csv');
    fileInfo = dir(filename_deco);
    if fileInfo.bytes > 0
        compare_deco = readmatrix(filename_deco, 'NumHeaderLines', 1);
    else
        % Handle empty case
        compare_deco = [];
    end

    
    count_70_ours = 0;
    count_139_ours = 0;
    count_else_ours = 0;



    count_70_deco = 0;
    count_139_deco = 0;
    count_else_deco = 0;


    if isempty(compare_ours)
        disp('skip this')
    else
        index_col = compare_ours(:, 1);
        pvalue_col = compare_ours(:, end);
        for i = 1:size(compare_ours, 1)
            idx = index_col(i);
            pval = pvalue_col(i);

            if (idx == 139 && pval < threshold)
                count_139_ours = count_139_ours + 1;
            elseif (idx == 70 && pval < threshold)
                count_70_ours = count_70_ours + 1;
            elseif (pval < threshold) 
                count_else_ours = count_else_ours + 1;
            end
        end
    end

    if isempty(compare_deco)
        disp('skip this')
    else
        index_col = compare_deco(:, 1);
        pvalue_col = compare_deco(:, 3);
        for i = 1:size(compare_deco, 1)
            idx = index_col(i);
            pval = pvalue_col(i);

            if (idx == 139 && pval < threshold)
                count_139_deco = count_139_deco + 1;
            elseif (idx == 70 && pval < threshold)
                count_70_deco = count_70_deco + 1;
            elseif (pval < threshold) 
                count_else_deco = count_else_deco + 1;
            end
        end
    end
    TPRours=(count_70_ours+count_139_ours)/2;
    FPRours=(count_else_ours)/(dim-2);
    TPRdeco=(count_70_deco+count_139_deco)/2;
    FPRdeco=(count_else_deco)/(dim-2);
    TPR_ours_list = [TPR_ours_list, TPRours];
    FPR_ours_list = [FPR_ours_list, FPRours];
    TPR_deco_list = [TPR_deco_list, TPRdeco];
    FPR_deco_list = [FPR_deco_list, FPRdeco];
end

clc;
mean(TPR_ours_list)
std(TPR_ours_list)

mean(FPR_ours_list*1000)
std(FPR_ours_list*1000)


mean(TPR_deco_list)
std(TPR_deco_list)


mean(FPR_deco_list*1000)
std(FPR_deco_list*1000)
