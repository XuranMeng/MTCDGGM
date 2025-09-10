% read all the data
% Data_iU : (sample size)*(q+1)
% Data_Z0:(sample size)*(p)
% Data_individualW:(j this is the index j of W_j)*(n)*(p-1 * q+1)
% Data_individualy:(j this is the index j of W_j)*(n)
function [Data_iU,Data_Z0,Data_individualW,Data_individualy]=readdata(n,p,q,nu,reptition)
 
    %this is the large design matrix, which is sparse
    interM = [];
    for j=1:(q+1)
        interM=[interM Z0.*iU(:,j)]; %these are the interaction terms between Z0 and iU
    end
    Data_iU(:,:)=iU;
    Data_Z0(:,:)=Z0;
    for j=1:p
        Xj=interM;
        Xj(:,((0:q)*p+j)) = [];
        y1=Z0(:,j);
        b=y1-mean(y1);
        Data_individualW(j,:,:)=Xj;
        Data_individualy(j,:)=b;
    end
    % Data_combiney(reptition,:,:)=b0;
    

end