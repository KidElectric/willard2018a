function out=stdErrMean(data,dim,flag)
% out=stdErrMean(data,dim,flag) %calculates mean +/- SEM
% see: stdErr() for S.E.M. function. ci95Mean() for 95% CI of mean
% out = n x 3 matrix of: mean, + stdErr, -stdErr of data
% Flag = 1 output m, and +/- stdErr ;
%      = 0, m and m +/- StdErr  (default)
if ~exist('dim','var') || isempty(dim)
    if size(data,1)==1
        dim=2;
    else
        dim=1;
    end
end

m=nanmean(data,dim);
se=stdErr(data,dim); %NaN compatible
if ~exist('flag','var') || flag==0
    out(:,1)=m;
    out(:,2)=(m + se);
    out(:,3)=(m + -se);
else
    out(:,1)=m;
    out(:,2)=se;
    out(:,3)=-se;
end