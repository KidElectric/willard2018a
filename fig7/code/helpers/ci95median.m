function out=ci95median(data,dim,flag)
% ci=ci95median(data,dim,flag)
if ~exist('dim','var') || isempty(dim)
    if size(data,1)==1
        dim=2;
    else
        dim=1;
    end
end

m=nanmedian(data,dim);
se=bootstrapCI95(data,flag);
out(:,1)=m;
out(:,2)=se(:,2)';
out(:,3)=se(:,1)';
