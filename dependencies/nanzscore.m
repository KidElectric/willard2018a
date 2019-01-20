function z=nanzscore(data,flag,dim)

if dim>2
    error('Only dim 1 and 2 supported')
end
ms=bsxfun(@minus,data,nanmean(data,dim));
% if dim==2
    z=bsxfun(@rdivide,ms,nanstd(data,flag,dim));
% else
%     z=bsxfun(@rdivide,ms',nanstd(data,flag,dim)')';
% end