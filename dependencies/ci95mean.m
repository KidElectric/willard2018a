function out=ci95mean(data,dim,flag)
% out=ci95mean(data,dim,flag)
% Calcululate 95% confidence interval of mean based on SEM.
% Note: Wrapper to correctly label out=stdErrMean(data,dim,flag)
% interval around mean. Updated to use ts = tinv([0.025  0.975],n-1)
% instead of 1.96*SE.
% Caculate mean (1) +/- 1.96*stdErr(2-3) of data (columns 1-3 of variable
% "out") =95% confidence interval. 2016
% Flag = 1 output m, and +/- 95% confidence interval;
%      = 0, m and m +/- 95% confidence interval  (default)
if ~exist('dim','var') || isempty(dim)
    if size(data,1)==1
        dim=2;
    else
        dim=1;
    end
end

m=nanmean(data,dim);
se=stdErr(data,dim); %NaN compatible
ts = tinv([0.025  0.975],size(data,dim)-1); 
if ~exist('flag','var') || flag==0
    out(:,1)=m;
    out(:,2)=(m+ts(2)*se)';
    out(:,3)=(m+ts(1)*se)';
else
    out(:,1)=m;
    out(:,2)=(ts(2)*se)';
    out(:,3)=(ts(1)*se)';
end