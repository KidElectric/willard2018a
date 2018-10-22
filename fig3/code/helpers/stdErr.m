function se=stdErr(data,dim,n)
% se=stdErr(data,dim,n)
%NaNs are excluded from estimate.

% Calculate standard error of the mean
% if ~exist('dim','var') || isempty(dim)
%     dim=2;
% end

if ~exist('n','var') || isempty(n)
    n=sum(~isnan(data),dim); %Number of observations
end
se=nanstd(data,0,dim)./sqrt(n);
