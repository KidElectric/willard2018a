function [h,pAll,mseO]=barMeanSig2(data,x,cols2Compare,colStr,excludeOutliers,...
    isPaired,anovaFirst,mc,wilcox,negError,plotMedians)
% [h,pAll,mseO]=barMeanSig2(data,x,cols2Compare,colStr,excludeOutliers,isPaired,anovaFirst,mc,wilcox,negError,plotMedians)
%Categorial scatter plot with category means, 95% confidence intervals and ttests.
%OUT:
% pAll{1}=ANOVA (even if anovaFirst=false)using anova1()
% pAll{2}=ttest value (paired or not)
% Uses plotSpread function.
%data - if cell, each cell is data to use for each category in x.
%       if matrix, Each column of data is a category corresponding to the labels found
%       in x. Rows are observations (NaN values are excluded from analysis).
%x - cell arra of strings labeling data categories in data variable columns. 
%   E.g. {'WhiskInAir','0.00mm','2.00mm','4.00mm','8.00mm','10.00mm'}
%cols2Compare - cell aray of column pairs to perform ttest2 comparisions
%          E.g. {[1,2],[3,2]}.  Same format used by sigstar() function
%colStr - for 3 bars: {'k','r','b'} plots bars with white fill, colored
%         outline. To color fill and outline: format: {{bar1Edge,bar1Fill}
%         ... etc}
%         i.e. {{'k','k'},{'r','r'},{'b','b'}} 
%excludeOutliers - leave empty []
%isPaired - are data compared paired data
%anovaFirst= default = true (do anova first before post-hoc tests
%mc = default true (control for multiple comparisons (length of
%   cols2compare) using sidak correction
%wilcox = if ==1 and ispaired == 1,   uses a paired signrank(), 
%          if ==1 and ispaired ==0,    uses ranksum().
%negError = by default 0 (only positive error bars plotted).
%            if == 1 positive and negative error bars plotted.
%plotMedians = plot median +/- boostrapped 95CI -- for use with wilcox
% data
hold on
errorBarMean=zeros(1,size(data,2));
errorBarSE=errorBarMean;
if ~exist('excludeOutliers','var') || isempty(excludeOutliers)
    excludeOutliers=0;
end
if ~exist('wilcox','var') || isempty(wilcox)
    wilcox=0;
end
if ~exist('anovaFirst','var') || isempty(anovaFirst)
    anovaFirst=1;
end
if ~exist('mc','var') || isempty(mc)
    mc=1;
end
if ~exist('isPaired','var') || isempty(isPaired)
    isPaired=false;
end

% d2=zeros(1,sum(~isnan(data(:))));
if iscell(data)
    a=max(cellfun(@length,data));
    y=nan(a,length(data));
    grpInd=zeros(size(y));
    for j=1:length(data)
       dTemp=data{j};
       ll=length(dTemp);
       y(1:ll,j)=dTemp;
       grpInd(:,j)=j;
    end
else
    y=nan(size(data));
    grpInd=zeros(size(data));
    for j=1:size(data,2)        
        y(:,j)=data(:,j);
        grpInd(:,j)=j;
    end
end

if ~exist('plotMedians','var') || isempty(plotMedians)
    plotMedians=false;
end

xn=1:size(y,2);
if plotMedians
    mseO=ci95median(y,1,1);
else
    mseO=ci95mean(y,1,1);
end
if ~exist('colStr','var') || isempty(colStr)
    colStr={'w'};
end
hold on
if size(colStr,2)>1    
    for i=1:size(y,2)
        if iscell(colStr{i})
            h(i)=bar(i,mseO(i,1),'EdgeColor',colStr{i}{1},'FaceColor',colStr{i}{2},'LineWidth',2);     
        else
            h(i)=bar(i,mseO(i,1),'EdgeColor',colStr{1,i},'FaceColor','none','LineWidth',2);     
        end
    end
else
    h=bar(xn,mseO(:,1),'EdgeColor',colStr{1,1},'FaceColor','w','LineWidth',1); %Default
end
if plotMedians
    mse=ci95median(y,1,0);
else
    mse=ci95mean(y,1,0);
end
% mse=ci95mean(y,1,0);
% low=zeros(size(mse(:,1)));
% high=zeros(size(mse(:,1)));
% low(mse(:,1)<0)=mse(mse(:,1)<0,3);%Lower error bars
% high(mse(:,1)>0)=mse(mse(:,1)>0,2);%Higher error bars
% herr=errorbar(xn,mse(:,1),zeros(size(mse(:,2))),mse(:,2),'LineStyle','none','Color','k');
for i=1:length(xn)
    herr(i,1)=plot([xn(i) xn(i)],[mse(i,1),mse(i,2)],'Color','k');
    if exist('negError','var') && negError==1
        herr(i,2)=plot([xn(i) xn(i)],[mse(i,1),mse(i,3)],'Color','k');
    end
end
% for k=1:length(herr)
    set(herr,'LineWidth',1,'UserData','Ex')
%     ch=get(herr,'Children');
%     for j=1:length(ch)
%         set(ch(j),'LineWidth',1.5,'UserData','Ex')
%     end
% end
%     
f=gcf
if ~isempty(cols2Compare) &&  strcmp(cols2Compare{1},'all')
    if wilcox==1
        [pAnova,~,stats] = kruskalwallis(y,x,'off');
    else
        [pAnova,~,stats]=anova1(y,x,'off');
    end
    f2=figure
    [c,m,h,nms] = multcompare(stats);
%     c
    figure(f)
else
    if wilcox==1
        [pAnova,~,stats] = kruskalwallis(y,x,'off');
    else
        [pAnova,~,stats]=anova1(y,[],'off');
    end
%     pAnova=anova1(y,[],'off');   
    p=nan(size(cols2Compare));
    if pAnova < 0.05 || (length(p)==1 && isPaired==true) || anovaFirst == 0%if doing 1 paired t-test, do not do anova first.
        for i=1:length(cols2Compare)
            data1=y(:,cols2Compare{i}(1));
            data2=y(:,cols2Compare{i}(2));
            if isPaired==true
                if wilcox==1
                    [p(i),~]=signrank(data1,data2);
                else
                    [~,p(i)]=ttest(data1,data2);
                end
            else
                if wilcox==1
                    [p(i),~]=ranksum(data1,data2);
                else
                    [~,p(i)]=ttest2(data1,data2);
                end
            end
        end
    else
        for i=1:length(cols2Compare)
            p(i)=pAnova;
        end
    end
   
    cc=vertcat(cols2Compare{:});
    if ~isempty(cols2Compare)
        for i= 1:length(cols2Compare)
            alpha=[1E-3 1E-2 0.05];
            if length(cols2Compare)>1 && mc==1
                %Sidak correction:
                m=round(sum(ismember(cc(:),cols2Compare{i}))/2);
                alpha=sidak(alpha,m);
            end
            try
                sigstar(cols2Compare{i},p(i),0,alpha)
            catch err
               throw(err) 
            end
        end
    end
    pAll{2}=p;
end

set(gca,'fontsize',16)
set(gca,'xtick',xn)
set(gca,'xticklabel',x)
pAll{1}=pAnova;
% hold off
% ylim([0 max(get(gca,'ylim'))])
% xlabel('Stimulus')
% ylabel('Mean Response (%, n=9-16 days)')
