function aw2018_figure_5(pn_local,pn_common,reProc)
% function aw2018_figure_5(pn_local,pn_common,reProc)
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

if reProc==1
    %%
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
        %     mnum=xlsread(pnfn,sheet,'B3:B1000');
        [~,mnum]=xlsread(pnfn,sheet,ids{i});
        %     if ismember(shees{i},{'UniDepl','1injIntactAsym','2injIntactAsym','5injIntactAsym'...
        %        '1injDeplAsym','2injDeplAsym','5injDeplAsym','UniCtl'})
        %         mnum=[
        %     else
        um=unique(mnum);
        mouseID=[mouseID; mnum];
        %     end
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
    
    % %add protocol duration as column 16:
    % data=[data nan(size(data,1),1)];
    % for i=1:length(dur_ID)
    %     ind=ismember(mouseID,dur_ID{i});
    %     data(ind,16)=dur_val(i);
    % end
    % mouseID0=mouseID;
    
    
    %%%%
    
    
    
    useResponse=data(:,9); % Velocity
    xlab={'%SB','FR','CV','Synch'};
    % [1,4,6,15];
    %include vel as variable for PCA:
    xlab=[xlab 'Vel'];
    % t=types(~xx & ind);
    
    % Load in raw %TH values
    
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
    
    display('Finished loading physiology')
    
    %% Generate physiology PC data for Figure 5
    clear params
    disp('Performing PCA')
    abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
        'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
    usePredictors=[1,4,6,15];
    xlab={'%SB','FR','CV','%Sync'};
    params.phys_predict_cols=usePredictors;
    params.phys_predict_labs=xlab;
    
    % grad=ismember(types,[1:4,13]);
    % grad_ba=ismember(types,[1:5,13]);
    asyn=ismember(types,[6:8,13]);
    ind=asyn;
    
    params.use_conditions_num=unique(types(ind));
    params.use_conditions_label=abrev(params.use_conditions_num);
    
    % xx=any(isinf(data(:,[usePredictors,9])) | isnan(data(:,[usePredictors,9])),2); %Exclude bad data and others
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
    params.reason_mice_ex=xlab(c);
    
    dat=data(~xx & ind,usePredictors);
    tempTH=th_values;
    useTH=th_values(~xx & ind);
    isSkew=[1 0 1 1];
    params.phys_take_log=isSkew;
    
    datRaw=dat;
    for i=1:length(isSkew)
        if isSkew(i)
            dat(:,i)=log10(dat(:,i));
            dat(isinf(dat(:,i)),i)=0;
        end
    end
    xlab={'%SB','FR','CV','%SyncP'};
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
    useData=data(~xx & ind,behave);
    
    useMouse=mouseID(~xx & ind);
    % durs=data(~xx & ind,16);
    t=types(~xx & ind);
    useLog=[0 1 1 1];
    params.behave_take_log=useLog;
    for i=1:length(useLog)
        if useLog(i)==1
            useData(:,i)=log10(useData(:,i));
            useData(isinf(useData(:,i)),i)=0;
        end
    end
    
    behaveUseData=useData;
    a=unique(useMouse);
    clear data_mouse score_mouse type_mouse keep_score keep_mouse ...
        data_mouse_raw p_dur th_mouse
    
    for j=1:length(a)
        mouse=ismember(useMouse,a{j});
        type_mouse(j)=mode(t(mouse));
        data_mouse(j,:)=nanmean(behaveUseData(mouse,:),1);
        for i=1:length(usePredictors)
            if isSkew(i)
                data_mouse_raw(j,i)=nanmedian(datRaw(mouse,i),1);
            else
                data_mouse_raw(j,i)=nanmean(datRaw(mouse,i),1);
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
        ints=combntns(1:(size(data_mouse,2)),2);
        for i=1:size(ints,1)
            data_mouse=[data_mouse data_mouse(:,ints(i,1)).*data_mouse(:,ints(i,2))];
        end
    end
    
    
    % We must deal with nan-values in behavior
    % Must identify, remove for PCA, then re-insert after pca:
    nan_ind=find(any(~isnan(data_mouse),2));
    temp=data_mouse(nan_ind,:);
    score_temp=nan(size(data_mouse));
    w = 1./var(temp);
    [wcoeff,score,latent,tsquared,explained_behav] = pca(temp,...
        'VariableWeights',w);
    behavCoeff=inv(diag(std(temp)))*wcoeff;
    
    score_temp(nan_ind,:)=score;
    
    allDat=[keep_mouse',num2cell(th_mouse),abrev(type_mouse)',...
        num2cell(keep_score(:,1)), num2cell(keep_score(:,2)),...
        num2cell(keep_score(:,3)),...
        num2cell(score_temp(:,1)), num2cell(score_temp(:,2)),...
        num2cell(data_mouse_raw),num2cell(data_mouse),...
        ];  %num2cell(p_dur)
    cols=[cellstr('Mouse'),cellstr('%TH'), cellstr('Type'), ...
        cellstr('Phys PC1 Score'), ...
        cellstr('Phys PC2 Score'),cellstr('Phys PC3 Score'), ...
        cellstr('Behav PC1 Score'),cellstr('Behav PC2 Score'),...
        xlab,'Vel','Rear','TotePole','WireHang']; %'Protocol Duration'
    allDat=[cols; allDat];
    
    % fnSave='4_behavs_all_states+raw_no_ints_2';
    fnSave='2019_01_14_V1';
    
    mkdir([pn_local 'output'])
    xlswrite([pn_local 'output' filesep fnSave '.xlsx'],allDat)
    behaveCoeff=wcoeff;
    mkdir([pn_local 'output'])
    save([pn_local 'output' filesep fnSave '.mat'],'allDat','physCoeff',...
        'behavCoeff','usePredictors','pn_local','fn','fn2','data_mouse',...
        'physUseData','params',...
        'useTH','explained_behav','explained_phys','ints')
    display('Finished')
end
%% Load these data
fn2='2019_01_14_V1';
[~,~,d]=xlsread([pn_local 'output' filesep sprintf('%s.xlsx',fn2)]);
load([pn_local 'output' filesep sprintf('%s.mat',fn2)])
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
phys(:,1)=phys(:,1);
physCoeff(:,1)=physCoeff(:,1);
disp('Finished loading')

%% Figure 5F-G 
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