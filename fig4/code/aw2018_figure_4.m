function aw2018_figure_4(pn)
% function aw2018_figure_4(pn)
%   pn = data path. if empty, assumes data are in a folder 'data' such that:
%   current directory is 'code,' and 'data' and 'code' are in 'fig' folder
% Load physiology data
% Read AW snr data .xlsx
if ~exist('pn','var') || isempty(pn)
    pn=['.' filesep '..' filesep 'data' filesep];
end
%
fn2='4_behavs_asyn_only+raw_no_ints_3';
%% Load these data
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
bFlipSign=1;
if bFlipSign==1
    behav=-behav; %flip sign of scores
    behavCoeff=-behavCoeff; %flip sign of coeffs
end
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

disp('Finished loading')

%% Figure 4C (manually added fits from 3B):
r=2; c=2;
figure
useLog=[0 0 0 0];
useCols=13:16;
x={0:100,25:100};
ylabs=d(1,useCols);
fits={'poly2','poly2','poly2','poly3';...
      'poly2','poly2','poly2','poly1'};
fitInd= {allgrad | ctl, asyn | ctl};
lims={[0 1.25],[0 2.25],[0 2.5],[0 3]};
for j=1:2
    for i=1:length(useCols)        
        dat=cell2mat(d(2:end,useCols(i)));
        ylab=ylabs{i};
        if useLog(i)==1
            dat=log10(dat+1);
            dat(isinf(dat))=0;
            ylab=['log' ylab];
        end        
        [xData, yData] = prepareCurveData( th(fitInd{j}),dat(fitInd{j}) );
        ft = fittype(fits{j,i});
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = 'off';
        opts.Normalize = 'on';
        [fitresult2, gof] = fit( xData, yData, ft, opts );
        
        subplot(r,c,i)
        if j==2
            plot(x{j},fitresult2(x{j}),'Color',[90 187 234]./255)
            
            hold on
            t=[mean(dat(ctl)),stdErr(dat(ctl),1)];
            scatter(mean(th(ctl)),t(1),'filled','d','MarkerFaceColor',[188 190 192]./255);
            plot([1 1].*mean(th(ctl)),[t(2) -t(2)]+t(1),'Color',[188 190 192]./255)
            plot([0 100],t(1).*[1 1],'--','Color',[188 190 192]./255)
            h=scatter(th(asyn),dat(asyn,1),10,'o','filled',...
                'MarkerFaceColor',[90 187 234]./255 );
            t=[mean(dat(ba)),stdErr(dat(ba),1)];
            if strcmp(fits{j,i},'poly1')
                title(sprintf('R^2: %1.3f',gof.rsquare))
            else
                title(sprintf('Adj R^2: %1.3f',gof.adjrsquare))
            end            
            xlabel('TH Remaining')
            ylabel(ylab)
            set(gca,'Xdir','reverse','xtick',0:25:100,'ytick',0:1:3.5)
            ylim(lims{i})
        else
            plot(x{j},fitresult2(x{j}),'--','Color',[151 36 123]./255)
            hold on
        end
    end
end
bi_Plot_Corrections
set(gcf,'pos',[680   572   414   406])
%% Figure 4D
x=0:100;
figure
r=2;c=1;
subinds=subplot_indexer(r,c);
ind=asyn| ctl ;
fittypes={'poly3','poly1','poly1'};
Robust={'off','Off','off'};

ylab={'V','R','TP','WH'};
reOrder=[1:4]; 
useInts=0;
d0={};
flipSign=1; %of all
usePCs=1:c;
out=behavCoeff(:,usePCs);
ind2 = asyn | allgrad | ba | ctl;
endStage = ind2 & th<= 5;
naive = ind2 & th>= 90;
endStageScore=phys(endStage,:);
naiveScore=phys(naive,:);
zscale=1;
lab={'5%','30','60','85%','Ctl'};
bins=[0 15 48 73 100];
useBins=[1 2 3 4 5];
lab=lab(useBins);
for usePC=1:c     
    subaxis(r,c,subinds(1,usePC),'SV',0.1,'PL',0.05,'MR',0.1)
%     if flipSign==1
%         out(:,usePC)=-out(:,usePC);
%     end  
    barh((out(reOrder,usePC)),'k')
    set(gca,'ytick',1:4,...
        'yticklabel',ylab(reOrder),'ydir','reverse')
    ylim([0.5 4.5])
    xlim([-1 1])    
    
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

%% Figure 4E

x=25:100;
figure
inds ={ctl | asyn};
fittypes={'poly2','poly2'};
for i=1:1
    %Fit and plot behavior PC1 scores vs th
    [xData, yData] = prepareCurveData( th(inds{i}),behav(inds{i},1) );
    ft = fittype( fittypes{i} );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'off';
    opts.Normalize = 'on';
    [fitresult2, gof] = fit( xData, yData, ft, opts );
    p2 = predint(fitresult2,x,0.95,'functional','off');
    hold on
    plot(x,fitresult2(x),'k')
    if i==1
        title(sprintf('Adj R-Sqr: %1.3f',gof.adjrsquare))
        t=[mean(behav(ctl,1),1),stdErr(behav(ctl,1),1)];
        h=scatter(mean(th(ctl)),t(1),'filled','d','MarkerFaceColor',[188 190 192]./255);
        plot([1 1].*mean(th(ctl)),[t(2) -t(2)]+t(1),'Color',[188 190 192]./255)
        plot(x,p2,'--k')
        plot([0 100],mean(behav(ctl,1)).*[1 1],'--k')
          scatter(th(asyn),behav(asyn,1),'o','filled','MarkerFaceColor',[90 187 234]./255)
    end
    xlabel('%TH Remaining')
    ylabel(sprintf('Behavior PC%d Score',1))
    set(gca,'Xdir','reverse','xtick',0:25:100)
end
bins=[0 15 48 73 100];
lab={'Dep','30','60','85','Ctl'};
bi_Plot_Corrections
% ylim([-0.5 2])
xlim([0 105])
set(gcf,'pos',[680   650   337   328])