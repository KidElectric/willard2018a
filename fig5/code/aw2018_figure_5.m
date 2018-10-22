function aw2018_figure_5(pn)
% function aw2018_figure_5(pn)
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
load([pn sprintf('%s.mat',fn2)],'behavCoeff','physCoeff')
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

%% Figure 5 E-H

x=25:100;
figure
r=2;c=4;
subinds=subplot_indexer(r,c);
ind=asyn | ctl ;
physFits={'poly3','poly1','poly1'};
Robust={'off','Off','off'};

ylab={'B','FR','IR','SYN'};
reOrder=[2 3 1 4]; 
useInts=0;
d0={};
% flipSign=0; %of all
usePCs=1:r;
if useInts==1
    reOrder=[reOrder length(reOrder)+ (1:size(ints,1))];
    for i=1:size(ints,1)
        newLab{i}=sprintf('%s x %s',ylab{ints(i,1)},ylab{ints(i,2)});
    end
    allLab=[ylab newLab];
else
    allLab=ylab;
end
out=physCoeff(:,usePCs);
ind2=asyn | ctl;
endStage=ind2 & th<=5;
naive=ind2 & th>=90;
endStageScore=phys(endStage,:);
naiveScore=phys(naive,:);
zscale=0.5;
lab={'5','30','60','85','Ctl'};
bins=[0 15 48 73 100];
useBins=[1 3 5];
lab=lab(useBins);

for usePC=1:r        
    subaxis(r,c,subinds(usePC,1),'SV',0.1,'PL',0.05,'MR',0.1)
    barh((physCoeff(reOrder,usePC)),'k')
    set(gca,'ytick',1:size(out,1),...
        'yticklabel',allLab(reOrder),'ydir','reverse')
    ylim([0.5 size(physCoeff,1)+0.5])
    xlim([-1 1])
    
    subaxis(r,c,subinds(usePC,2),'PL',0)
    [~,ind1]=histc(th,bins);
    temp=[];
    for ii=1:length(useBins)
        d0=phys(ind1==useBins(ii) & ind,usePC);
        temp=[temp;mean(d0*out(:,usePC)',1)];

    end
    imagesc(fliplr(temp(:,reOrder)'))
    set(gca,'ytick',1:size(out,1),...
        'yticklabel',{},'xtick',1:length(bins),...
        'xticklabel',fliplr(lab),'clim',[-1 1].*zscale)
    colormap colorblind_cm
    
    %Fit and plot behavior PC1 scores vs th
    subaxis(r,c,subinds(usePC,3:4),'PL',0.1)
    [xData, yData] = prepareCurveData( th(ind),phys(ind,usePC) );
    
    ft = fittype( physFits{usePC} );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = Robust{usePC};
    opts.Normalize = 'off';
    [fitresult2, gof] = fit( xData, yData, ft, opts );
    p2 = predint(fitresult2,x,0.95,'functional','off');
    hold on
    plot(x,fitresult2(x),'k')
    plot(x,p2,'--k')
    
    t=[mean(phys(ctl,usePC),1),stdErr(phys(ctl,usePC),1)];
    h=scatter(mean(th(ctl)),t(1),'filled','d','MarkerFaceColor',[188 190 192]./255);
    plot([1 1].*mean(th(ctl)),[t(2) -t(2)]+t(1),'Color',[188 190 192]./255)
    

    plot([0 100],mean(phys(ctl,usePC)).*[1 1],'--k')
    scatter(th(asyn),phys(asyn,usePC),'o','filled','MarkerFaceColor',[90 187 234]./255)
  
    if usePC==1
        title(sprintf('Adj R^2: %1.3f',gof.adjrsquare))
    elseif usePC==2 
        title(sprintf('R^2: %1.3f',gof.rsquare))
        xlabel('%TH Remaining')
    end
    ylabel(sprintf('Physiology PC%d',usePC))
    set(gca,'Xdir','reverse','xtick',0:25:100)
    xlim([0 105])
end


bi_Plot_Corrections

set(gcf,'pos',[ 680   591   474   387])

%% Load these data
fn2='4_behavs_all_states+raw_no_ints_3'; %For figure 5I-J
[~,~,d]=xlsread([pn sprintf('%s.xlsx',fn2)]);
load([pn sprintf('%s.mat',fn2)],'behavCoeff','physCoeff')
fr=cell2mat(d(2:end,10));
sb=cell2mat(d(2:end,9));
cv=cell2mat(d(2:end,11));
syn=cell2mat(d(2:end,12));
phys=cell2mat(d(2:end,4:6));
zphys1=nanzscore(phys(:,1),[],1);
zphys2=nanzscore(phys(:,2),[],1);
zphys3=nanzscore(phys(:,3),[],1);
behav=cell2mat(d(2:end,7:8));
bFlipSign=0;
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
%% Figure 5I-J
% Need to load in a file with grad and asyn run simultaneously in PC or
% z-score all 
x0=25:100;
x2=0:100;
asyn_inds=asyn | ctl;
grad_inds=allgrad | ctl;
inds={asyn_inds,grad_inds};
dep = (allgrad | ba ) & th < 15;
cols={'c','m';'c','m'};
physFits={'poly3','poly1';...
    'poly3','poly1'};
behaveFit={'poly3','poly3'};
Robust={'off','off','off'};

for usePC=1:2
    figure   
    clear h
    plot([25 100],[0 0],'--k')
    hold on
    base=mean(phys(ctl,:),1);    
    for i=1:2 % 1 == alphasyn, 2 == 6ohda
        if i==1
            x=x0;
        else
            x=x2;
        end
        ind=inds{i};
        [xData, yData] = prepareCurveData( th(ind),phys(ind,usePC) );
        ft = fittype( physFits{i,usePC} );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = Robust{usePC};
        opts.Normalize = 'off';
        [fitresult2, gof] = fit( xData, yData, ft, opts );
        %         p2 = predint(fitresult2,x,0.95,'functional','off');
        y=fitresult2(x);
        h(i,1)=plot(x,y-y(end),cols{usePC,i});
        
        %    Plot behavior
%         if usePC==1
            [xData, yData] = prepareCurveData( th(ind),behav(ind,1) );
            ft = fittype( behaveFit{i} );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
            opts.Robust = Robust{1};
            opts.Normalize = 'on';
            [fitresult2, gof] = fit( xData, yData, ft, opts );
            %             p2 = predint(fitresult2,x,0.95,'functional','off');
            y=fitresult2(x);
            h(i,2)=plot(x,y-y(end),['--' cols{1,i}]);
%         end
    end
    
    t=[mean(phys(dep,1))-base, stdErr(phys(dep,1),1)];
    scatter(mean(th(dep)),t(1),'filled','^','MarkerFaceColor',[80 45 140]./255);
    plot([1 1].*mean(th(dep)),[t(2) -t(2)]+t(1),'Color',[80 45 140]./255)
    
    % legend([h(1,:) h(2,:)],{'Asyn Phys PC1', 'Asyn Behav PC1','Grad Phys PC1','Grad Behave PC1'})
    xlim([0 105])
    
    set(gca,'Xdir','reverse','xtick',0:25:100)
    % ylim([-3 1.5])
    xlabel('%TH Remaining')
    ylabel('Norm Score')
    bi_Plot_Corrections
    set(gcf,'pos',[  680   344   532   473])
    
end
