function ind=equalize_class_quant(y,classIDs)
% ind=equalize_class_quant(y,classIDs)
% Selectively downsample other categories to match minimum category
%y= vector of all classes to equalize
%classIDs= vector of class ids,
% Example:
% y=datasample([1 2 3 4],100);
% aa=histc(y,1:4)
% aa =
%     28    22    36    14 <- classes have different numbers of
%     observations
% ind=equalize_class_quant(y,1:4);
% y1=y(ind);
% histc(y1,1:4)
% ans =
%     14    14    14    14 <- equalized to minimum observation number

count0=histc(y,classIDs);
[s,si]=sort(count0);  %Be careful that 
m=find(s>0,1,'first');
mn=s(m);
scv=classIDs(si); 
rc=scv(m+1:end); %sorted classes to downsample
crc=s(m+1:end);
ind=true(size(y));
if any(m)
    %Selectively downsample other categories to match minimum category
    for i=1:length(rc)
        c=(rc(i));
        n0=count0(classIDs==c); %Original count
%         n0=count0(si==c); %Original count
        d=n0-mn; %Difference
        ind0=find(y==c);
        removeInd=datasample(ind0,d,'Replace',false);%ind0(randperm(n0,d));
        ind(removeInd)=false;
    end
end