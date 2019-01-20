function ci=bootstrapCI95(X,flag)
% ci=bootstrapCI95(X) return 95% confidence interval from 1000 iteration
% bootstrap of median(X) (dataSample choose 200).
val=zeros(1,1000);
ci=zeros(size(X,2),2);
for j=1:size(X,2)
    for i=1:1000
        tt=datasample(X(:,j),100);
        val(i)=nanmedian(tt);
    end
    ci(j,:)=[ quantile(val,0.025), quantile(val,0.975)];
    if exist('flag','var') && flag==1
        ci(j,:)=ci(j,:)-nanmedian(val);
    end
end