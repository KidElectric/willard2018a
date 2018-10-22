function aw2018_figure_2(pn)
% function aw2018_figure_2(pn)
%   pn = data path. if empty, assumes data are in a folder 'data' with 'code' (one level up)
%% Load physiology data
% Read AW snr data .xlsx
if ~exist('pn','var') || isempty(pn)
    pn=['.' filesep '..' filesep 'data' filesep];
end
%
fn='SNrInVivoData_wRest_FINAL_v7.xlsx';
fn2='SynchronyData_6-6-18_useable_updated.xlsx';
fn3='GradualProtocolDuration.xlsx';
fn4='4_behavs_grad_ba+raw_no_ints_2';

%Protocol Duration:
[~,~,gradpdur]=xlsread([pn fn3],'Gradual','A2:D30');
[~,~,asynpdur]=xlsread([pn fn3],'Alpha-Syn','A2:D11');
dur_ID=[gradpdur(:,2) ; asynpdur(:,2)];
dur_val=cell2mat([gradpdur(:,4) ; asynpdur(:,4)]);

%New synch:
[~,~,psync]=xlsread([pn fn2],'Sheet1','A2:D81');
ps_ID=psync(:,2);
ps_val=cell2mat(psync(:,3));
contra=cell2mat(psync(:,4)); %2 == contra

sheets={'Grad85','Grad60','Grad30','Gradual5'...
    'Acute',...
    'ASyn30','ASyn60','ASyn85',...
    'UniDepl',...
    'Ipsi(Depl)Asym','ContraAsym','UniCtl','Ctl'};

rows={'M3:Z274','M3:Z320','M3:Z320','M3:Z295',...
    'L3:Y247',...
    'M3:Z152','M3:Z157','M3:Z100',...
    'M3:Z138',...
    'M3:Z360','M3:Z267','M3:Z96','I3:V264',};

ids={'AA3:AA1000','AA3:AA1000','AA3:AA1000','AA3:AA1000',...
    'Z3:Z1000',...
    'AA3:AA1000','AA3:AA1000','AA3:AA1000',...
    'AA3:AA1000',...
    'AA3:AA1000','AA3:AA1000','AA3:AA1000','W3:W1000'};

type=1:length(sheets);

cols={'%SpikesInBurst','MeanFreq','CVisi','Synch'};
data=[];
types=[];
pnfn=[pn fn];
mouseID={};
mcount=0;
for i=1:length(sheets)
    sheet=sheets{i};
    dat=xlsread(pnfn,sheet,rows{i});
    data=[data;dat];
    types=[types; ones(size(dat,1),1)*type(i)];
    [~,mnum]=xlsread(pnfn,sheet,ids{i});
    um=unique(mnum);
    mouseID=[mouseID; mnum];
    if length(mnum) ~= size(dat,1)
        fprintf('Warning: Mismatch in %s.\n',sheet)
    end
    mcount=mcount+length(um);
end

%add percent synchronous pairs as column 15
data=[data nan(size(data,1),1)];
for i=1:length(ps_ID)
    mID=mouseID;
    if contra(i)==2
        useID={sprintf('%s_%d',ps_ID{i},11), sprintf('%s_%d',ps_ID{i},12)};
    else
        useID=ps_ID{i};
    end
    ind=ismember(mID,useID);
    tt(i)=sum(ind);
    data(ind,15)=ps_val(i);
end

%add protocol duration as column 16:
data=[data nan(size(data,1),1)];
for i=1:length(dur_ID)
    ind=ismember(mouseID,dur_ID{i});
    data(ind,16)=dur_val(i);
end

%%%%
useSeparateCond_As_Mice=1;
if useSeparateCond_As_Mice==1
    for i=1:length(mouseID)
        mouseID{i}=[mouseID{i} sprintf('_%d',types(i))];
    end
end
useResponse=data(:,9); % Velocity
xlab={'%SB','FR','CV','Synch'};

%include vel as variable for PCA:
xlab=[xlab 'Vel'];

% Load in raw %TH values

%All Categories:
% sheets={'Grad85','Grad60','Grad30','Gradual5'...
%     'Acute',...
%     'ASyn30','ASyn60','ASyn85',...
%     'UniDepl','1injIntactAsym','2injIntactAsym','5injIntactAsym'...
%     '1injDeplAsym','2injDeplAsym','5injDeplAsym','UniCtl','Ctl'};
rows={'G3:I274','G3:I320','G3:I320','G3:I295',...
    'F3:H247',...
    'G3:I152','G3:I157','G3:I100',...
    'G3:I138',...
    'G3:I360','G3:I267','G3:I96','G3:I264','none'};
data2=[];
side=['L','R'];
skeep=[];
for i=1:length(sheets)
    sheet=sheets{i};
    [~,~,dat]=xlsread(pnfn,sheet,rows{i});
    if i==length(sheets) %control
        data2=[data2;ones(size(dat,1),1)*100];
        temp=cell(size(dat,1),1);
        temp(:)={'LR'};
        skeep=[skeep; temp];
    else
        sides=dat(:,1);
        lval=cell2mat(dat(:,2));
        rval=cell2mat(dat(:,3));
        dat=zeros(size(dat,1),1);
        dat(ismember(sides,side(1)))=lval(ismember(sides,side(1)));
        dat(ismember(sides,side(2)))=rval(ismember(sides,side(2)));
        data2=[data2;dat];
        skeep=[skeep;sides];
    end
end

useSeparateHems_As_Mice=0;
if useSeparateHems_As_Mice==1
    for i=1:length(mouseID)
        mouseID{i}=[mouseID{i} skeep{i}];
    end
end

th_values=data2;
display('Finished loading physiology')

%% Figure 2D Confirmation of statistics
%types variable are numbers 1-13 corresponding to the following groups:
abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
    'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
useOrder=[13, 1:4];
grad=ismember(types,useOrder);
params={'%SB','FR','CV','Synch'};
paramCols=[1,4,6,15];
cvISI=data(grad,paramCols(3));
an=mouseID(grad);
groups=types(grad);

for i=1:length(useOrder)
    datCell{i}=cvISI(groups==useOrder(i));
end

[p,tbl,stats] = kruskalwallis(cvISI,groups);
figure;
combs=nchoosek(1:5,2);
for i=1:length(combs)
    ind1=groups==useOrder(combs(i,1));
    ind2=groups==useOrder(combs(i,2));
    [p(i),h]=ranksum(cvISI(ind1),cvISI(ind2));
end
sigs=p < sidak(0.05,size(combs,1)+1);
sigs=p < 0.05/(size(combs,1)+1);
combs(sigs,:)
[c,m,h,gnames]=multcompare(stats,'Display','on','Ctype','bonferroni')
% sigs=c(c(:,end)<sidak(0.05,comb_n_choose_k(5,2)),[1 2])

%% Load behavior data in:
% pn='C:\Users\Gittis\Dropbox\Gittis Lab Data\AWSNR\';
% pn='C:\Dropbox\Gittis Lab Data\AWSNR\';
% fn2='4_behavs_with_all_states+raw_all_ints_2';
% fn2='4_behavs_grad_only+raw_no_ints_2';
[~,~,d]=xlsread([pn sprintf('%s.xlsx',fn4)]);
load([pn sprintf('%s.mat',fn4)])

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

% Only use paired unilateral mice:
uniContra=ismember(type,'UCtrl');
uniIpsi=ismember(type,'UDep');
anTemp=unique([anID(uniIpsi);anID(uniContra)]);
coan=anTemp(ismember(anTemp,anID(uniContra)) & ismember(anTemp,anID(uniIpsi)));
isPaired= ismember(anID,coan);
uniIpsi= uniIpsi & isPaired;
uniContra= uniContra & isPaired;

asymContra=ismember(type,{'1IIA','2IIA','5IIA','AsymCtrl'});
asymIpsi=ismember(type,{'1IDA','2IDA','5IDA','AsymDep'});
anTemp=unique([anID(asymIpsi);anID(asymContra)]);
coan=anTemp(ismember(anTemp,anID(asymContra)) & ismember(anTemp,anID(asymIpsi)));
isPaired= ismember(anID,coan);
asymIpsi= asymIpsi & isPaired;
asymContra= asymContra & isPaired;
phys(:,1)=-phys(:,1);
physCoeff(:,1)=-physCoeff(:,1);
display('Finished')
flipSign=1; %Assign custom meaning to relative PC scores, negative = more like depleted
%% Figure 2G-J: Gradual PC1 and PC2 vs TH

x=0:100;
figure
r=2;c=4;
subinds=subplot_indexer(r,c);
ind=allgrad | ctl ;
fittypes={'poly3','poly1','poly1'};
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
ind2=asyn | allgrad | ba | ctl;
endStage=ind2 & th<=5;
naive=ind2 & th>=90;
endStageScore=phys(endStage,:);
naiveScore=phys(naive,:);
zscale=0.5;
lab={'Dep','30','60','85','Ctl'};
bins=[0 15 48 73 100];
useBins=[1 3 5];
lab=lab(useBins);


for usePC=1:r
    
    subaxis(r,c,subinds(usePC,1),'SV',0.1,'PL',0.05,'MR',0.1)
    %     subplot(r,c,subinds(usePC,1))
    barh((physCoeff(reOrder,usePC)),'k')
    set(gca,'ytick',1:size(out,1),...
        'yticklabel',allLab(reOrder),'ydir','reverse')
    ylim([0.5 size(physCoeff,1)+0.5])
    xlim([-1 1])
    if flipSign==1
        out(:,usePC)=-out(:,usePC);
    end
    
    
    subaxis(r,c,subinds(usePC,2),'PL',0)
    [~,ind1]=histc(th,bins);
    temp=[];
    for ii=1:length(useBins)
        d0=phys(ind1==useBins(ii) & ind,usePC);
        %         t0=d0*out(:,usePC)';
        %         t1=mean(t0,1)./std(t0,[],1);
        %         temp=[temp;t1];
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
    
    ft = fittype( fittypes{usePC} );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = Robust{usePC};
    opts.Normalize = 'off';
    [fitresult2, gof] = fit( xData, yData, ft, opts );
    p2 = predint(fitresult2,x,0.95,'functional','off');
    hold on
    plot(x,fitresult2(x),'k')
    plot(x,p2,'--k')
    
    t=[mean(phys(ctl,usePC),1),stdErr(phys(ctl,usePC),1)];
    scatter(mean(th(ctl)),t(1),'dk','filled')
    plot([1 1].*mean(th(ctl)),[t(2) -t(2)]+t(1),'k')
    plot([0 100],mean(phys(ctl,usePC)).*[1 1],'--k')
    scatter(th(allgrad),phys(allgrad,usePC),'om','filled')
    t=[mean(phys(ba,usePC),1),stdErr(phys(ba,usePC),1)];
    scatter(mean(th(ba)),t(1),'^b','filled')
    plot([1 1].*mean(th(ba)),[t(2) -t(2)]+t(1),'b')
    
    
    if usePC==1
        title(sprintf('Adj R-Sqr: %1.3f',gof.adjrsquare))
    elseif usePC == 2
        title(sprintf('R-Sqr: %1.3f',gof.rsquare))
    end
    if usePC==2
        xlabel('%TH Remaining')
    end
    ylabel(sprintf('Physiology PC%d',usePC))
    set(gca,'Xdir','reverse','xtick',0:25:100)
    xlim([0 105])
end


bi_Plot_Corrections

set(gcf,'pos',[ 680   591   474   387])
