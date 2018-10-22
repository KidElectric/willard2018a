function aw2018_figure_3(pn)
% function aw2018_figure_2(pn)
%   pn = data path. if empty, assumes data are in a folder 'data' such that:
%   current directory is 'code,' and 'data' and 'code' are in 'fig' folder
% Load physiology data
% Read AW snr data .xlsx
if ~exist('pn','var') || isempty(pn)
    pn=['.' filesep '..' filesep 'data' filesep];
end
%

%% Figure 3A data load (note takes a while to load):
t=dir([pn 'Raw data*.xlsx']);
fns={t(:).name};
disp('Loading movement (x,y) data')
for i=1:length(fns)
    fn=fns{i};
    pnfn=[pn fn];
    disp(pnfn)
    [~,~,d]=xlsread(pnfn);    
    keep(i).x=cell2mat(d(41:end,3));
    keep(i).y=cell2mat(d(41:end,4));
    keep(i).vel=cell2mat(d(41:end,14));
end
display('Finished')

%% Figure 3A: B&W 5-minutes:
val={'G5','G60','Ctl'};
r=3;c=1;
reorder=[3 2 1];
markerSize=5;
dur=[0 60*5];
useTimes=[dur; dur+60*3; dur + 60*13]; %seconds
    figure;
    for i=1:length(fns)
        subplot(r,c,i)
        useTime=useTimes(i,:);
        time=linspace(0,length(keep(i).x)/(29.97),length(keep(i).x));
        ind=time > useTime(1) & time < useTime(end);
        vel=keep(reorder(i)).vel(ind);
        markerSize=(1./((vel)+ 0.001)) ;
        markerSize=((markerSize-min(markerSize))/max(markerSize))*30+5;
        scatter(keep(reorder(i)).x(ind)+35,keep(reorder(i)).y(ind)-22,...
            3,'ok','filled',...
            'MarkerEdgeColor','none','MarkerFaceAlpha',1)
        colormap(colorblind_cm)
        set(gca,'clim',[0 1.5],'xtick',-20:10:20,'ytick',-20:10:20,...
            'xlim',[-25 25],'ylim',[-25 25])
        axis off
    end
    bi_Plot_Corrections
    if c>r
        set(gcf,'pos',[507   420   963   198])
    else
        set(gcf,'pos',[   507    69   171   592])
    end
% end
%% Load data for remaining figures:
fn2='4_log_behavs_grad_ba+raw_no_ints_3';
[~,~,d]=xlsread([pn sprintf('%s.xlsx',fn2)]);
load([pn sprintf('%s.mat',fn2)])
fr=cell2mat(d(2:end,10));
sb=cell2mat(d(2:end,9));
cv=cell2mat(d(2:end,11));
syn=cell2mat(d(2:end,12));
phys=cell2mat(d(2:end,4:6));
zphys1=nanzscore(phys(:,1),[],1);
zphys2=nanzscore(phys(:,2),[],1);
zphys3=nanzscore(phys(:,3),[],1);
behav=cell2mat(d(2:end,7:8));
zbehave1=nanzscore(behav(:,1),[],1);
type=d(2:end,3);
px=(phys(:,1));
th=cell2mat(d(2:end,2));
th(th<=0)=1e-10;
anID=d(2:end,1);
for i=1:length(anID)
   anID{i}=anID{i}(1:8);    
end
asyn=ismember(type,{'A30','A60','A85'});
grad=ismember(type,{'G30','G60','G85'});
ctl=ismember(type,{'Ctrl'});
ba=ismember(type,{'BA'});
isEndStage=ismember(type,{'BA','G05'});
allgrad=ismember(type,{'G05','G30','G60','G85'});

phys(:,1)=-phys(:,1);
physCoeff(:,1)=-physCoeff(:,1);
display('Finished loading')
%% Figure 3B: Plot all behavior grad, ctl and BA 

r=2; c=2;
figure
useLog=[0 0 0 0];
useCols=13:16;
x=0:100;
ylabs=d(1,useCols);
fits={'poly2','poly2','poly2','poly3'};
fitInd= allgrad | ctl ;
lims={[0 1.25],[0 2.25],[0 2.5],[0 3]};
for i=1:length(useCols)    
    dat=cell2mat(d(2:end,useCols(i)));
    ylab=ylabs{i};
    if useLog(i)==1
        dat=log10(dat+1);
        dat(isinf(dat))=0;
        ylab=['log' ylab];
    end    
    [xData, yData] = prepareCurveData( th(fitInd),dat(fitInd) );
    ft = fittype(fits{i});
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'off';
    opts.Normalize = 'on';
    [fitresult2, gof] = fit( xData, yData, ft, opts );

    subplot(r,c,i)
    plot(x,fitresult2(x),'Color',[151 36 123]./255)
    hold on
    t=[mean(dat(ctl)),stdErr(dat(ctl),1)];
    scatter(mean(th(ctl)),t(1),'filled','d','MarkerFaceColor',[188 190 192]./255);
    
    plot([1 1].*mean(th(ctl)),[t(2) -t(2)]+t(1),'Color',[188 190 192]./255)    
    plot([0 100],t(1).*[1 1],'--','Color',[188 190 192]./255)    
    h=scatter(th(allgrad),dat(allgrad,1),10,'o','filled',...
        'MarkerFaceColor',[151 36 123]./255 );
    t=[mean(dat(ba)),stdErr(dat(ba),1)];
    scatter(mean(th(ba)),t(1),'filled','^','MarkerFaceColor',[80 45 140]./255);
    plot([1 1].*mean(th(ba)),[t(2) -t(2)]+t(1),'Color',[80 45 140]./255)
    title(sprintf('Adj R-Sqr: %1.3f',gof.adjrsquare))
    xlabel('%TH Remaining')
    ylabel(ylab)
    set(gca,'Xdir','reverse','xtick',0:25:100,'ytick',0:1:3.5)
    ylim(lims{i})
end
bi_Plot_Corrections
set(gcf,'pos',[680   572   414   406])

%% Figure 3C:
x=0:100;
figure
r=2;c=1;
subinds=subplot_indexer(r,c);
ind=allgrad | ctl ;
fittypes={'poly3','poly1','poly1'};
Robust={'off','Off','off'};

ylab={'V','R','TP','WH'};
reOrder=[1:4]; 
useInts=0;
d0={};
flipSign=0; %of all
usePCs=1:c;
out=behavCoeff(:,usePCs);
ind2=asyn | allgrad | ba | ctl;
endStage=ind2 & th<=5;
naive=ind2 & th>=90;
endStageScore=phys(endStage,:);
naiveScore=phys(naive,:);
zscale=1;
lab={'5%','30','60','85%','Ctl'};
bins=[0 15 48 73 100];
useBins=[1 3 5];
lab=lab(useBins);


for usePC=1:c     
    subaxis(r,c,subinds(1,usePC),'SV',0.1,'PL',0.05,'MR',0.1)
    barh((behavCoeff(reOrder,usePC)),'k')
    set(gca,'ytick',1:4,...
        'yticklabel',ylab(reOrder),'ydir','reverse')
    ylim([0.5 4.5])
    xlim([-1 1])
    if flipSign==1
        out(:,usePC)=-out(:,usePC);
    end  
    
    subaxis(r,c,subinds(2,usePC),'PL',0)
    [~,ind1]=histc(th,bins);
    temp=[];
    for ii=1:length(useBins)
        d0=behav(ind1==useBins(ii) & ind,usePC);
        temp=[temp;mean(d0*out(:,usePC)',1)];

    end
    imagesc(fliplr(temp(:,reOrder)'))
    set(gca,'ytick',1:size(out,1),...
        'yticklabel',ylab(reOrder),'xtick',1:length(bins),...
        'xticklabel',fliplr(lab),'clim',[-1 1].*zscale)
    colormap colorblind_cm
    
end
bi_Plot_Corrections
set(gcf,'pos',[    680   522   116   295])
%% Figure 3D

x=0:100;
figure
ind= ctl|  allgrad ;
i=1;
fittypes={'poly2','poly1'};
%Fit and plot behavior PC1 scores vs th
[xData, yData] = prepareCurveData( th(ind),behav(ind,i) );
ft = fittype( fittypes{i} );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'off';
opts.Normalize = 'off';
[fitresult2, gof] = fit( xData, yData, ft, opts );
p2 = predint(fitresult2,x,0.95,'functional','off');
hold on
plot(x,fitresult2(x),'k')
plot(x,p2,'--k')

t=[mean(behav(ctl,i),1),stdErr(behav(ctl,i),1)];
h=scatter(mean(th(ctl)),t(1),'filled','d','MarkerFaceColor',[188 190 192]./255);
plot([1 1].*mean(th(ctl)),[t(2) -t(2)]+t(1),'Color',[188 190 192]./255)
    

plot([0 100],mean(behav(ctl,1)).*[1 1],'--k')
scatter(th(allgrad),behav(allgrad,1),'o','filled','MarkerFaceColor',[151 36 123]./255)

t=[mean(behav(ba,i),1),stdErr(behav(ba,i),1)];
h=scatter(mean(th(ba)),t(1),'filled','^','MarkerFaceColor',[80 45 140]./255);
plot([1 1].*mean(th(ba)),[t(2) -t(2)]+t(1),'Color',[80 45 140]./255)

title(sprintf('Adj R-Sqr: %1.3f',gof.adjrsquare))
xlabel('%TH Remaining')
ylabel(sprintf('Behavior PC%d Score',i))
set(gca,'Xdir','reverse','xtick',0:25:100)
    
bins=[0 15 48 73 100];
lab={'Dep','30','60','85','Ctl'};
bi_Plot_Corrections

xlim([0 105])
set(gcf,'pos',[680   650   337   328])

%% Figure 3E: PC1 and PC2 fits compared to behavior 
zphys=bsxfun(@rdivide,phys-mean(phys(:)),std(phys(:)));
abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
    'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
r=2;c=2;
x=0:100;
figure
subinds=subplot_indexer(r,c);

ind= ctl |  allgrad;
zd={zphys(ind,1), zphys(ind,2), zphys(ind,3)};
fit_types={'poly3','poly1','poly1'};
Robust={'off','off','off'};
clear cc
for i=1:c
    % NORMALIZED scores and fits:
    [xData, yData] = prepareCurveData( th(ind), zd{i} );
    ft = fittype(fit_types{i} );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
     
    opts.Normalize = 'off';
    opts.Robust=Robust{i};
    [fitresult1, gof1] = fit( xData, yData, ft, opts );
    p1 = predint(fitresult1,x,0.95,'functional','off');
    [xData, yData] = prepareCurveData( th(ind), zbehave1(ind) );
    ft = fittype( 'poly2' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'LAR';
    opts.Normalize = 'on';
    [fitresult2, gof2] = fit( xData, yData, ft, opts );
    p2 = predint(fitresult2,x,0.95,'functional','off');

    y1=fitresult1(x)-fitresult1(100);
    h=plot(x,y1,'k');
    hold on
    h.DisplayName= 'Physiology PC1';
    y2=fitresult2(x) - fitresult2(100);
    h=plot(x,y2,'b');
    h.DisplayName= sprintf('Behavior PC%d',i);
    plot([0 100],[0 0],'--b')
    ylabel(sprintf('Norm Score',i))
    xlabel('%TH Remaining')
     set(gca,'Xdir','reverse','xtick',0:25:100)
    ylim([-3 1.5])
end
% bi_linkaxes('x')

 set(gca,'Xdir','reverse','xtick',0:25:100)
set(gcf,'pos',[917   431   264   276])
bi_Plot_Corrections
% 
