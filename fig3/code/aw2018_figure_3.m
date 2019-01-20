function aw2018_figure_3(pn_local,pn_common,reProc)
% function aw2018_figure_3(pn_local,pn_common,reProc)
%
%   pn_local = data path for figure. if empty, assumes data are in a folder 'data' such that:
%   current directory is 'code,' and 'data' and 'code' are in 'fig' folder
%
%   pn_local = common data folder.  if empty, assumes data are in a folder 'common_data' such that:
%   current directory is 'code,' and 'common_data' is 2 levels above 'code'
%
%   reProc = rerun processing from common_data. default is 0.
%
% Load physiology data
% Read AW snr data .xlsx
if ~exist('pn_local','var') || isempty(pn_local)
    pn_local=['.' filesep '..' filesep 'data' filesep];
end
if ~exist('pn_common','var') || isempty(pn_common)
    pn_common=['.' filesep '..' filesep '..' filesep 'common_data' filesep];
end
if ~exist('reProc','var') || isempty(pn_common)
    reProc=0;
end

%% Figure 3A data load (note takes a while to load):
t=dir([pn_local 'Raw data*.xlsx']);
fns={t(:).name};
disp('Loading movement (x,y) data')
for i=1:length(fns)
    fn=fns{i};
    pnfn=[pn_local fn];
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
if reProc==1
    %% Load raw data for gradual behavior PCA
    disp('Loading data...')
    fn='SNrInVivoData_wRest_FINAL_v8_NewSynchronyIndex.xlsx';
    fn2='SynchronyData_12-04-18_useable.xlsx';
    
    %New synch:
    [~,~,psync]=xlsread([pn_common fn2],'Sheet1','A2:D87');
    ps_ID=psync(:,2);
    isNans= cellfun(@ischar,psync(:,3));
    temp=psync(:,3);
    temp(isNans)={nan};
    ps_val=cell2mat(temp);
    unilat=cell2mat(psync(:,4)); %1= bilat, 3== ipsi, 2 == contra
    % abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
    %     'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
    %All Categories:
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
    pnfn=[pn_common fn];
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
    mID=mouseID;
    params.useSeparateCond_As_Mice=1;
    if params.useSeparateCond_As_Mice==1
        for i=1:length(mouseID)
            mouseID{i}=[mouseID{i} sprintf('_%d',types(i))];
        end
    end
    
    %
    %add percent synchronous pairs as column 15
    data=[data nan(size(data,1),1)];
    
    % abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
    %     9 'UDep',10 'AsymDep',11 'AsymCtrl',12 'UCtrl',13 'Ctrl'};
    for i=1:length(ps_ID)
        if unilat(i)==2
            useID={sprintf('%s_%d',ps_ID{i},11), sprintf('%s_%d',ps_ID{i},12)};
            ind=ismember(mouseID,useID);
        elseif unilat(i)==3
            useID={sprintf('%s_%d',ps_ID{i},9), sprintf('%s_%d',ps_ID{i},10)};
            ind=ismember(mouseID,useID);
        else
            useID=ps_ID{i};
            ind=ismember(mID,useID);
        end
        
        tt(i)=sum(ind);
        data(ind,15)=ps_val(i);
    end
    
    if any(tt==0)
        error('Mouse synchrony data dropped erroneously.')
    end
    
    %%%%
    xlab={'%SB','FR','CV','Synch'};
    
    %All Categories:
    
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
    
    params.useSeparateHems_As_Mice=0;
    if params.useSeparateHems_As_Mice==1
        for i=1:length(mouseID)
            mouseID{i}=[mouseID{i} skeep{i}];
        end
    end
    
    th_values=data2;
    
    disp('Finished loading physiology')
    
    %% Generate physiology PC data for Figure 3H-I
    clear params
    disp('Performing PCA')
    abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
        'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
    usePredictors=[1,4,6,15];
    phys_labs={'%SB','FR','CV','%Sync'};
    params.phys_predict_cols=usePredictors;
    params.phys_predict_labs=phys_labs;
    
    % grad=ismember(types,[1:4,13]);
    grad_ba=ismember(types,[1:5,13]);
    % asyn_grad30=ismember(types,[3,6:8,13]);
    ind=grad_ba;
    
    params.use_conditions_num=unique(types(ind));
    params.use_conditions_label=abrev(params.use_conditions_num);
    
    xx=any(isinf(data(:,usePredictors)) | isnan(data(:,usePredictors)),2); %Exclude bad physiology data only
    
    % Which mice having units exlcuded & why?
    uex=unique(mouseID(xx));
    clear unitsEx colsBd typeM
    for i=1:length(uex)
        mi=ismember(mouseID,uex{i});
        unitsEx(i)=sum(xx & mi)./sum(mi);
        [r,c]=find(isinf(data(mi,usePredictors)) | isnan(data(mi,usePredictors)));
        colsBd{i}=unique(c);
        typeM{i}=abrev{types(find(mi,1,'first'))};
    end
    params.excluded_mice=uex;
    params.reason_mice_ex=phys_labs(c);
    
    dat=data(~xx & ind,usePredictors);
    tempTH=th_values;
    useTH=th_values(~xx & ind);
    isSkew=[1 0 1 1];
    params.phys_take_log=isSkew;
    
    datRaw=dat; %Before Log()
    for i=1:length(isSkew)
        if isSkew(i)
            dat(:,i)=log10(dat(:,i));
            dat(isinf(dat(:,i)),i)=0;
        end
    end
    phys_labs={'%SB','FR','CV','%SyncP'};
    t=types(~xx & ind);
    params.types_after_indexing=unique(t);
    useData=dat;
    params.phys_interactions=0;
    if params.phys_interactions==1
        ints=combntns(1:(size(useData,2)),2);
        for i=1:size(ints,1)
            useData=[useData useData(:,ints(i,1)).*useData(:,ints(i,2))];
        end
    else
        ints=[];
    end
    
    w = 1./var(useData);
    [wcoeff,physscore,latent,tsquared,explained_phys] = pca(useData,...
        'VariableWeights',w,'centered','on');
    physCoeff = inv(diag(std(useData)))*wcoeff;
    physUseData=useData;
    behave=[9:10, 13, 14]; %Velocity, rearing, total pole, wire hang
    params.behav_cols=behave;
    params.behave_labs={'Vel','Rear','Total Pole','Wire Hang'};
    behave_raw=data(~xx & ind,behave);
    
    useMouse=mouseID(~xx & ind);
    % durs=data(~xx & ind,16);
    t=types(~xx & ind);
    useLog=[0 1 1 1];
    params.behave_take_log=useLog;
    clear behave_log
    for i=1:length(useLog)
        if useLog(i)==1
            behave_log(:,i)=log10(behave_raw(:,i));
            behave_log(isinf(behave_log(:,i)),i)=0;
        else
            behave_log(:,i)=behave_raw(:,i);
        end
    end
    
    a=unique(useMouse);
    clear data_mouse score_mouse type_mouse keep_score keep_mouse ...
        data_mouse_raw p_dur th_mouse  mouse_behave_raw  mouse_phys_raw mouse_behave_norm
    
    for j=1:length(a)
        mouse=ismember(useMouse,a{j});
        type_mouse(j)=mode(t(mouse));
        mouse_behave_norm(j,:)=nanmean(behave_log(mouse,:),1); %Log taken
        for i=1:length(usePredictors)
            mouse_behave_raw(j,i)=mode(behave_raw(mouse,i));
            if isSkew(i)
                mouse_phys_raw(j,i)=nanmedian(datRaw(mouse,i),1);
            else
                mouse_phys_raw(j,i)=nanmean(datRaw(mouse,i),1);
            end
        end
        keep_mouse{j}=a{j};
        keep_score(j,:)=nanmean(physscore(mouse,:),1);
        th_mouse(j,1)=nanmean(useTH(mouse));
    end
    
    sign_invert=0; %Depends on whether sign convention assigned to chosen
    %data subset matches sign convention of the paper
    if sign_invert==1
        keep_score=-keep_score;
        physCoeff=-physCoeff;
    end
    
    params.bhv_interactions=0;
    % useInts=0; %Decide whether interactions used for behavior data (typically not)
    if params.bhv_interactions==1
        ints=combntns(1:(size(mouse_behave_norm,2)),2);
        for i=1:size(ints,1)
            mouse_behave_norm=[mouse_behave_norm mouse_behave_norm(:,ints(i,1)).*mouse_behave_norm(:,ints(i,2))];
        end
    end
    
    
    % We must deal with nan-values in behavior
    % Must identify, remove for PCA, then re-insert after pca:
    nan_ind=find(any(~isnan(mouse_behave_norm),2));
    temp=mouse_behave_norm(nan_ind,:);
    score_temp=nan(size(mouse_behave_norm));
    w = 1./var(temp);
    [wcoeff,score,latent,tsquared,explained_behav] = pca(temp,...
        'VariableWeights',w);
    behavCoeff=inv(diag(std(temp)))*wcoeff;
    
    score_temp(nan_ind,:)=score;
    
    allDat=[keep_mouse',num2cell(th_mouse),abrev(type_mouse)',...
        num2cell(keep_score(:,1)), num2cell(keep_score(:,2)),...
        num2cell(keep_score(:,3)),...
        num2cell(score_temp(:,1)), num2cell(score_temp(:,2)),...
        num2cell(mouse_phys_raw),num2cell(mouse_behave_raw),...
        ];  %num2cell(p_dur)
    cols=[cellstr('Mouse'),cellstr('%TH'), cellstr('Type'), ...
        cellstr('Phys PC1 Score'), ...
        cellstr('Phys PC2 Score'),cellstr('Phys PC3 Score'), ...
        cellstr('Behav PC1 Score'),cellstr('Behav PC2 Score'),...
        phys_labs,'Vel','Rear','TotePole','WireHang']; %'Protocol Duration'
    allDat=[cols; allDat];
    % Save data to output folder
    fnSave='2019_01_08_V2';
    
    mkdir([pn_local 'output'])
    xlswrite([pn_local 'output' filesep fnSave '.xlsx'],allDat)
    behaveCoeff=wcoeff;
    mkdir([pn_local 'output'])
    save([pn_local 'output' filesep fnSave '.mat'],'allDat','physCoeff',...
        'behavCoeff','usePredictors','pn_local','fn','fn2','mouse_phys_raw',...
        'mouse_behave_raw','physUseData','params',...
        'useTH','explained_behav','explained_phys','ints')
    display('Finished')
end
%% Load data for remaining figures:

fn4='2019_01_08_V2';
[~,~,d]=xlsread([pn_local 'output' filesep sprintf('%s.xlsx',fn4)]);
load([pn_local 'output' filesep sprintf('%s.mat',fn4)])


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
use=ismember(type,{'G30','G60','G85'});
ctl=ismember(type,{'Ctrl'});
ba=ismember(type,{'BA'});
isEndStage=ismember(type,{'BA','G05'});
allgrad=ismember(type,{'G05','G30','G60','G85'});

phys(:,1)=-phys(:,1);
physCoeff(:,1)=-physCoeff(:,1);
display('Finished loading')

%% Plot 3B-E (Binned): Plot raw behavior grad, ctl and BA
groups={'Ctrl','G85','G60','G30','G05','BA'}; %bilateral acute at end
xlabs={'Ctl','85%','60%','30%','5%','A5%'};
%Generate color scale:
col=zeros(length(groups),3);
for j=1:length(groups)
        if j==1
            t=[188 190 192]./255; %Gray
        elseif j==2
            t=[243 207 227]./255; %
        elseif j==3            
            t=[230 164 201]./255; %
        elseif j==4
            t=[208 90 161]./255; %
        elseif j==5
            t=[177 65 152]./255; %            
        else
            t=[80 45 140]./255; %Purple
        end
        col(j,:)=t;
end
r=2; c=2;
figure
% useLog=[1 1 1 1];
useLog=[0 0 0 0];
useCols=13:16;
x=0:100;
ylabs=d(1,useCols);
fits={'poly2','poly2','poly2','poly3'};
fitInd= allgrad | ctl ;
lims={[0 1.25],[0 2.25],[0 2.5],[0 3]};
out={};
cc=1;
for i=1:length(useCols)    
    dat=cell2mat(d(2:end,useCols(i)));
    ylab=ylabs{i};
    if useLog(i)==1
        dat=log10(dat+1);
        dat(isinf(dat))=0;
        ylab=['log' ylab];
    end    

    subplot(r,c,i)
    hold on
    for j=1:length(groups)
        ind=ismember(type,groups{j});
        if useLog(i)==1
            m(j)=nanmean(dat(ind));
        else
            m(j)=nanmedian(dat(ind));
        end
        h=bar(j,m(j),'k');
        set(h,'FaceColor',col(j,:))
        h=plot(j,dat(ind),'ks');
        set(h,'MarkerFaceColor','w')
        cc=cc+1;
    end   
    xlabel('%TH Remaining')
    ylabel(ylab)
    set(gca,'xtick',1:length(groups),'xticklabel',xlabs)
%      set(gca,'xtick',1:length(groups),'xticklabel',xlabs,'ytick',0:1:3.5)
%     ylim(lims{i})
end
bi_Plot_Corrections
set(gcf,'pos',[   680   433   760   545])
%% Output data

out={};
cc=1;
for i=1:length(useCols)    
    dat=cell2mat(d(2:end,useCols(i)));
    mid=d(2:end,1);
    ylab=ylabs{i};
    if useLog(i)==1
        dat=log10(dat+1);
        dat(isinf(dat))=0;
        ylab=['log' ylab];
    end    
    for j=1:length(groups)
        ind=find(ismember(type,groups{j}));
        for k=1:length(ind)
            out{cc,1}=mid{ind(k)};
            out{cc,2}=ylab;
            out{cc,3}=groups{j};
            out{cc,4}=dat(ind(k));
            cc=cc+1;
        end
    end   
end
disp('Finished')
out=[{'Mouse','Behavior', 'Depletion','Value'};out];
fnSave='grad_behave_by_mouse_v2_raw';
xlswrite([pn_local 'output' filesep fnSave '.xlsx'],out)
%% Figure 3F:
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
%% Figure 3G

x=0:100;
figure
ind= ctl|  allgrad  ;
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

t=[nanmean(behav(ba,i),1),stdErr(behav(ba,i),1)];
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

%% Figure 3H: Behavior vs. TH vs. Physiology. Note: physiology fit as X, dopamine fit as Y and transposed for consistent visualization.
use={'G85','G60','G30','G05','Ctrl','BA'};

ind=ismember(type,use);

physn=-phys(ind,1);
bhvn=behav(ind,1);
zeroBehave=median(bhvn(ismember(type(ind),'Ctrl')));
% zeroBehave=median(bhvn(ismember(type(ind),'Ctrl')));
thn=th(ind);
[xData, yData, zData] = prepareSurfaceData( thn, physn, bhvn);

% Set up fittype and options.
ft = fittype( 'poly23' );

% Fit model to data.
[fitresult1, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Create a figure for the plots.
figure( 'Name', 'untitled fit 1' );
n1=100;
n2=100;
dat1=linspace(-5,105,n1);
dat2=linspace(-1.5,1.5,n2);
tempZ=zeros(n1,n2);
for i=1:n1
    for j=1:n2
        tempZ(i,j)=fitresult1([dat1(j),dat2(i)]);
    end
end

% h = plot( fitresult1, [xData, yData], zData, 'Style', 'Contour' );
% h=contourf(dat1,dat2,logmodulus(tempZ'-zeroBehave),[-1.3:0.2:1.3]);
% temp
ll=-3.5;
aa=tempZ-zeroBehave;
aa(aa<ll)=ll;
h=contourf(dat1,dat2,aa,[-3.5:.75:2]);
hold on
scatter(xData,yData,60,[0 0 0],'filled')
colormap(colorblind_cm)
h=colorbar;
% set(gca,'clim',[-2.5 2.5])
h.Label.String='Behavior PC1';
xlabel('TH Remaining')
ylabel('Phys PC1')
grid off
set(gca,'xdir','reverse')
hold on
title(sprintf('Adj-R^2 = %1.3f',gof.adjrsquare))
[xData, yData] = prepareCurveData( thn, physn );
ft = fittype( 'poly3' );
[fitresult, gof] = fit( xData, yData, ft );
hold on
plot(0:100,fitresult(0:100),'--k')
set(gca,'xtick',0:25:100)
bi_Plot_Corrections


