function aw2018_figure_7(pn_local,pn_common,reProc)
% function aw2018_figure_7(pn_local,pn_common,reProc)
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
    display('Finished')
    
    
    %% Compile all PC scores for summary
    disp('Performing PCA')
    abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
        'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
    usePredictors=[1,4,6,8,15]; %What if include neuron-specific synchrony index
    xlab={'%SB','FR','CV','NSync','%Sync'};
    params.phys_predict_cols=usePredictors;
    params.phys_predict_labs=xlab;
    
    grad=ismember(types,[1:4,13]);
    grad_ba=ismember(types,[1:5,13]);
    asyn=ismember(types,[6:8,13]);
    ind=ismember(types,[1:13]); %All conditions or use subset
    
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
    isSkew=[1 0 1 0 1];
    params.phys_take_log=isSkew;
    
    datRaw=dat;
    for i=1:length(isSkew)
        if isSkew(i)
            dat(:,i)=log10(dat(:,i));
            dat(isinf(dat(:,i)),i)=0;
        end
    end
    % xlab={'%SB','FR','CV','Synch'};
    xlab={'%SB','FR','CV','NSync','%SyncP'};
    t=types(~xx & ind);
    params.types_after_indexing=unique(t);
    useData=dat;
    params.phys_interactions=1;
    useInts=1; %Use interactions for physiology data?
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
        %     p_dur(j,:)=mode(durs(mouse,1));
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
    fnSave='2018_12_10_V1';
    
    mkdir([pn_local 'output'])
    xlswrite([pn_local 'output' filesep fnSave '.xlsx'],allDat)
    behaveCoeff=wcoeff;
    mkdir([pn_local 'output'])
    save([pn_local 'output' filesep fnSave '.mat'],'allDat','physCoeff',...
        'behavCoeff','usePredictors','pn_local','fn','fn2','fn3','data_mouse',...
        'physUseData','params',...
        'useTH','explained_behav','explained_phys','ints','useInts')
    display('Finished')
end
%% Load data for Figure 7A, B, E, F and Supp Fig 1A
fn='2018_12_10_V1.mat'; %See
ddd=load([pn_local fn]);
d=ddd.allDat;
display('Finished loading.')

%% Supp Figure 1A
th=cell2mat(d(2:end,2));
phys1=cell2mat(d(2:end,4));
behv1=cell2mat(d(2:end,7));
bilat_ind=ismember(d(2:end,3),{'G85','G60','G30','G05','BA','A30','A60','A85','Ctrl'});
bilat_or_uni={'G85','G60','G30','G05','BA','A30','A60','A85','Ctrl','UDep','AsymDep'}; %'AsymDep'
figure
ind=ismember(d(2:end,3),bilat_or_uni);
physn=-phys1(ind);
thn=th(ind);
labs=d(2:end,3);
labs=labs(ind);

preds={'Int','B','FR','IR','NSYN','SYN'};
if ddd.useInts==1
    ylab=preds(2:end);
    for i=1:length(ddd.ints)
        newLab{i}=sprintf('%s x %s',ylab{ddd.ints(i,1)},ylab{ddd.ints(i,2)});
    end
    allLab=[ylab newLab];
else
    allLab=predictors;
end

usePC=2;
r=1;c=usePC;
reord=[2 3 1 5 4];
allLab(1:length(reord))=allLab(reord);
for i=1:usePC
    if i==1
        out=-ddd.physCoeff(:,i);
    end
    out(1:length(reord))=out(reord);
    subplot(r,c,i)
    barh(out,'k')
    set(gca,'ytick',1:size(out,1),...
        'yticklabel',allLab,'ydir','reverse')
    xlim([-1 1])
    hold on
    title(sprintf('%2.1f%% explained, PC%d',ddd.explained_phys(i),i))
end
%    set(gcf,'pos',[ 309        1104         634         510])
bi_Plot_Corrections

%% Figure 7A
th=cell2mat(d(2:end,2));
phys1=cell2mat(d(2:end,4));
bilat_or_uni={'G85','G60','G30','G05','BA','A30','A60','A85','Ctrl','UDep','AsymDep'}; %'AsymDep'
figure
ind=ismember(d(2:end,3),bilat_or_uni);
physn=-phys1(ind);
thn=th(ind);
labs=d(2:end,3);
labs=labs(ind);

col=zeros(length(thn),3);
groups={{'G85','G60','G30','G05'},{'BA'},{'Ctrl'},{'A30','A60','A85'},...
    {'UDep'},{'AsymDep'}};
for i=1:length(groups)
    f=find(ismember(labs,groups{i}));
    if i==1 %Grad 6OHDA
        for ii=1:length(f)
            col(f(ii),:)=[177 65 152]./255; %purple
        end
    elseif i==2 % 6OHDA Acutge ES
        for ii=1:length(f)
            col(f(ii),:)=[102 45 145]./255; %dark purple
        end
    elseif i==3
        for ii=1:length(f)
            col(f(ii),:)=[240 240 240]'./255; %grey
        end
    elseif i==4
        for ii=1:length(f)
            col(f(ii),:)=[163 222 249]'./255; %blue
        end
    elseif i==5 %Uni ipsi
        for ii=1:length(f)
            col(f(ii),:)=[0 104 56]'./255; %dark green
        end
    elseif i==6 %Asym Ipsi
        for ii=1:length(f)
            col(f(ii),:)=[0 104 56]'./255; %also dark green
            %             col(f(ii),:)=[0 104 56]'.*2./255; %light green
        end
    end
end

scatter(thn,physn,80,col,'filled')
ylabel('Phys PC1'); xlabel('%TH')
hold on
[xData, yData] = prepareCurveData( thn, physn );

% Set up fittype and options.
ft = fittype( 'poly5' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
plot(0:100,fitresult(0:100),'k')
p2 = predint(fitresult,0:100,0.95,'functional','off');
plot(0:100,p2,'k')
title(sprintf('Adj R^2 = %1.2f',gof.adjrsquare))
set(gca,'xdir','reverse','xtick',0:25:100)


bi_Plot_Corrections
set(gcf,'pos',[680   357   440   621])

tind=thn < 75 & thn > 35;
m=mean(thn(tind));
localMin=55;
pm=mean(physn(tind));
scatter(m,pm,70,'k','filled')
plot([localMin-20 localMin+20],[pm pm],'k')

%% Figure 7B
% Plot PC1 vs. TH by High Medium Low
th=cell2mat(d(2:end,2));
th=abs(th);
phys1=cell2mat(d(2:end,4));
use={'G85','G60','G30','G05','BA','A30','A60','A85','Ctrl','UDep','AsymDep'};
figure
ind=ismember(d(2:end,3),use);
physn=-phys1(ind);
thn=th(ind);

col=zeros(size(th));
col(th> 75)=1; %Ctrl
col(th <= 75 & th > 35)=2; %Intermediate capt
col(th <= 35) = 3; % High
scatter(thn,physn,80,col(ind),'filled')
ylabel('Phys PC1'); xlabel('%TH')
set(gca,'xdir','reverse','xtick',0:25:100)
bi_Plot_Corrections
set(gcf,'pos',[   680   708   602   270])
if reProc==1
    %% Fit model for 7C
    abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
        'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
    
    typeTemp=zeros(size(types));
    
    %If not using TH by hems, average across hemispheres first:
    th_mouse_val=th_values;
    um=unique(mouseID);
    nn=zeros(size(th_mouse_val));
    for i=1:length(um)
        ind=ismember(mouseID,um{i});
        temp=mean(unique(abs(th_values(ind))));
        th_mouse_val(ind)=temp;
        nn(ind)=sum(ind);
    end
    
    type_include=ismember(types,[1:10,13]);% Asym and uni IPSI included
    
    useG=[1 2 3];
    
    % +/- 20% TH from local min (TH=55)
    typeTemp(th_mouse_val > 75)=1; %Ctrl
    typeTemp(th_mouse_val <= 75 & th_mouse_val >  35 )=3; %Intermediate capt
    typeTemp(th_mouse_val <=  35 ) = 2; % High
    
    % 4 group version (Supp Figure 1C).
    % typeTemp(th_mouse_val > 75)=1; %Ctrl
    % typeTemp(th_mouse_val <= 75 & th_mouse_val > 35)=2; %Intermediate capt
    % typeTemp(th_mouse_val <= 35 & th_mouse_val > 10)=4; %gets confused with ES
    % typeTemp(th_mouse_val <= 10) = 3; % High
    
    ut=unique(typeTemp);
    useResponse=ut(ut~=0);
    
    typeTemp(~type_include)=0; %Exclude unilaterals for now
    
    % typeTemp(nn < 10)=0; %Exclude low-neuron mice
    um=unique(mouseID);
    for i=1:length(um)
        tt(i)=length(unique(typeTemp(ismember(mouseID,um{i}))));
    end
    
    usePredictors=[1, 4, 6, 8, 15]; %Synchrony index = ; Percent synchronous pairs = 15
    useCol=usePredictors;
    useData=data(:,useCol); %Put Ctrl in right-most column to be "reference category"
    useInts=1; %Use interactions? 1 = yes, 0 = no
    
    %Include Behavior
    useBehavior=0;
    if useBehavior == 1
        behave=[9, 10, 13, 14]; %Velocity, rearing, total pole, wire hang
        usePredictors=[usePredictors, behave];
        tempData=data(:,behave);
        useData=[useData, tempData];
    end
    
    %Include TH as a predictor
    useTH=0; %Since we are predicting TH, don't include TH as predictor
    if useTH==1
        useData=[useData, th_values];
    end
    
    originalTypes=types;
    useTypes=typeTemp;
    typeInd=ismember(useTypes,useResponse);
    
    useMID=mouseID;
    
    % Only keep neurons of correct type and with all required data fields:
    originalTypes(any(isnan(useData),2) | ~typeInd)=[];
    useTypes(any(isnan(useData),2) | ~typeInd)=[];
    useMID(any(isnan(useData),2) | ~typeInd,:)=[];
    useData(any(isnan(useData),2) | ~typeInd,:)=[];
    uInputs=unique(useMID);
    collectVals=zeros(size(useMID,1),1);
    
    useTypes=useTypes+1000;
    for i=1:length(useResponse)
        useTypes(useTypes== (useResponse(i)+1000 ))=i;
    end
    
    clear neuronsPerMouse typePerMouse collectVals origPerMouse th_mouse
    for i=1:length(uInputs)
        ind=ismember(useMID,uInputs{i});
        neuronsPerMouse(i)=sum(ind);
        collectVals(ind)=neuronsPerMouse(i);
        typePerMouse(i)=mean(useTypes(ind)); %all have same value
        origPerMouse(i)=mean(originalTypes(ind));
        th_mouse(i)=mean(th_values(ind));
    end
    
    if useTH==1
        isSkew=[1 0 1 1 1 0];
        zScore=[1 1 1 1 1 1];
        isCategorical=[0 0 0 0 0];
    elseif useBehavior==1
        isSkew=[1 0 1 1 1 0 1 1];
        zScore=[1 1 1 1 1 1 1 1];
        isCategorical=[0 0 0 0 0 0 0 0];
    else
        if length(usePredictors)==4
            isSkew=[1 0 1 1  ];
            zScore=[1 1 1 1 ];
        elseif length(usePredictors)==5
            isSkew=[1 0 1 1 1];
            zScore=[1 1 1 1 1];
            %         isSkew=[1 1 0 1 1];
            %         zScore=[1 1 1 1 1];
        elseif length(usePredictors)==6
            isSkew=[1 1 0 1 1 1];
            zScore=[1 1 1 1 1 1];
        end
        %     isCategorical=[0 0 0 0 ];
    end
    
    for i=1:size(useData,2)
        if isSkew(i)
            useData(:,i)=log10(useData(:,i));
            useData(isinf(useData(:,i)),i)=0;
        end
        if zScore(i)
            useData(:,i)=nanzscore(useData(:,i),[],1);
        end
    end
    if useBehavior==1
        useCol=1:4;
    else
        useCol=1:size(useData,2);
    end
    
    if useInts==1
        ints=combntns(useCol,2);
        for i=1:size(ints,1)
            useData=[useData useData(:,ints(i,1)).*useData(:,ints(i,2))];
        end
        useCol=1:size(useData,2);
    end
    
    perms=500; %Note: set to 500 for final analyses but use 10 or so to test
    clear pc conf_mat_counts conf_mat_norm bKeep conf_mat_norm_keep ...
        pk bk correct actual predicted predicted_null bkn pkn pc_null
    type=useResponse;
    nstim=length(type);
    pKeep=[]; pKeep_null=[]; bKeepNull=[];
    fitType='mrfit';
    options.alpha=1;
    % options.parallel=1;
    n=length(uInputs);
    pk=[];
    type_names=abrev(useResponse);
    targetN=50; %round(median(neuronsPerMouse));
    runNull=1;
    local_pred=zeros(length(uInputs),perms);
    local_correct=local_pred; local_pred_null=local_pred;
    local_correct_null=local_pred;
    disp('Starting models...')
    for ii=1:n % For each mouse, drop it out, fit model on remaining mice etc.
        ex=uInputs{ii};    %Cell to exclude in jack-knife
        useMouse=find(~ismember(uInputs,ex));
        use=~ismember(useMID,ex);
        clear pred_tot
        act=typePerMouse(ismember(uInputs,ex));
        bKeep=[];
        for i=1:perms
            tic
            % 1) Randomly downsample remaining cells to have equal examples of
            % each class
            ec=equalize_class_quant(typePerMouse(useMouse),1:nstim);
            useMouseID=uInputs(useMouse(ec)); %Bingo.
            
            % 2) Randomly resample cells of each mouse to targetN
            cj=[];
            for j=1:length(useMouseID)
                ind=ismember(useMID,useMouseID{j});
                cc=find(ind);
                if sum(ind) < targetN
                    cj=[cj; datasample(cc,targetN,'Replace',true)]; %Up-sample
                else
                    cj=[cj; datasample(cc,targetN,'Replace',false)]; %Down-sample
                end
            end
            permPredictors=useData(cj,:);
            permResponse=useTypes(cj);
            iind=find(ismember(useMID,ex));
            if length(iind) < targetN
                ti=datasample(find(ismember(useMID,ex)),targetN,'Replace',true);
            else
                ti=datasample(find(ismember(useMID,ex)),targetN,'Replace',false);
            end
            testDat=useData(ti,:);
            if strcmp(fitType,'glmnet')
                %No longer supported
                error('glmnet no longer supported')
            else
                [B,dev,stats] = mnrfit(permPredictors,permResponse,'Interactions','on',...
                    'model','nominal');
                prob=mnrval(B,testDat,stats,'model','nominal','Interactions','on');%'model','ordinal','ineractions','on');
                pl=nansum(prob,1);
                
                %Winner take all:
                if any(pl)
                    pred=find(pl==max(pl));
                else
                    warning('Nothing predicted!')
                    pred=randperm(3,1);
                end
                
                pKeep= stats.p;
                
                % Null model:
                if runNull==1
                    permClassNull=permResponse(randperm(length(permResponse),length(permResponse)));
                    [B_null,dev_null,stats_null] = mnrfit(permPredictors,permClassNull,'Interactions','on',...
                        'model','nominal');
                    prob=mnrval(B_null,testDat,stats_null);%'model','ordinal','ineractions','on');
                    pl=nansum(prob,1);
                    
                    %Winner take all:
                    pred_null=find(pl==max(pl));
                    pKeep_null=cat(3,pKeep_null, stats_null.p);
                    bKeepNull=cat(3,bKeepNull,B_null);
                    
                    local_correct_null(ii,i)=pred_null==act;
                    local_pred_null(ii,i)=pred_null;
                end
                
            end
            
            local_correct(ii,i)=pred==act;
            local_pred(ii,i)=pred;
            bKeep=cat(3,bKeep, B);
            
            %Confusion matrix contruction
            %Rows = Actual, Columns = Predicted
            conf=hist3([act pred],'edges',{1:nstim,1:nstim});
            conf_mat_counts(:,:,i)=conf;
            
            %Normalize
            normC=histc(act,1:nstim);%Normalization constants for each stimulus type (number of times each stimulus was used as a test trial)
            conf_norm=bsxfun(@rdivide,conf,normC);
            conf_norm(isinf(conf_norm))=nan;
            conf_mat_norm(:,:,i)=conf_norm;
            %         fprintf('Jack-knife mouse %d, Perm %d: ',ii,i)
            %         fprintf('%s Mouse predicted as %s.\t',type_names{act},type_names{pred})
            %         toc
        end
        pc(ii)=sum(local_correct(ii,:))/perms * 100;
        fprintf('\n%2.1f%% correct, for mouse %d of %d.\n',pc(ii),ii,length(uInputs))
        
        conf_mat_norm_keep(:,:,ii)=mean(conf_mat_norm,3);
        bk(:,:,ii)=mean(bKeep,3);
        
        %Evaluate holdout
        actual(ii) = act;
        predicted(ii)= mode(local_pred(ii,:));
        if runNull==1
            pc_null(ii)=sum(local_correct_null(ii,:))/perms * 100;
            predicted_null(ii)=mode(local_pred_null(ii,:));
            bkn(:,:,ii)=mean(bKeepNull,3);
            fprintf('%2.1f%% null correct, for mouse %d of %d.\n',pc_null(ii),ii,length(uInputs))
        end
        name=ex;
        if strcmp(fitType,'glmnet')
            %no longer supported
        else
            pk(:,:,ii)=pKeep; %ALL
            if runNull==1
                pkn(:,:,ii)=median(pKeep_null,3);
            end
        end
        fprintf('%d mouse predicted as %d.\n\n',act,predicted(ii))
    end
    
    
    conf=hist3([actual' predicted'],'edges',{1:nstim,1:nstim});
    
    %Normalize
    normC=histc(actual,1:nstim)';%Normalization constants for each stimulus type (number of times each stimulus was used as a test trial)
    conf_norm=bsxfun(@rdivide,conf,normC);
    conf_norm(isinf(conf_norm))=nan;
    conf_mat_norm_total=conf_norm; %Actual x Predicted
    
    if runNull==1
        conf_null=hist3([actual' predicted_null'],'edges',{1:nstim,1:nstim});
        normC=histc(actual,1:nstim)';%Normalization constants for each stimulus type (number of times each stimulus was used as a test trial)
        conf_norm=bsxfun(@rdivide,conf_null,normC);
        conf_norm(isinf(conf_norm))=nan;
        conf_null_mat_norm_total=conf_norm;
    end
    
    beta=bk;
    beta_p=pk;
    
    responses=type_names;
    predictors=['int' cols];
    disp('Finished')
    %
    if runNull==1 && perms > 100
        beta_p_null=pkn;
        beta_null=bkn;
        save([pn_local , sprintf('mnrfit_%d_perm_jackknife_100-75-35_ints_pred_phys_med_rel_rerun.mat',perms)],...
            'pc','pc_null','conf_mat_norm_keep','beta','beta_p','beta_p_null','beta_null','conf_mat_norm_total',...
            'responses','useResponse','origPerMouse','local_pred','local_pred_null','nstim',...
            'predictors','sheets','cols','pnfn','fitType','conf_null_mat_norm_total',...
            'predicted','actual','ints','useInts','predicted_null','params')
        disp('Saved')
    end
    
    figure
    imagesc(conf_mat_norm_total)
end
%% Load model from 7C
load([pn_local 'mnrfit_500_perm_jackknife_100-75-35_ints_pred_phys_med_rel.mat'])

%% Figure 7C: Confusion Matrix
%Plot confusion matrix from all iteraions? Polot conusion matrix mode
local_act=bsxfun(@times,ones(size(local_pred)),actual');
conf=hist3([local_act(:),local_pred(:)],...
    'edges',{1:nstim,1:nstim}); %Confusion matrix performed on all iterations
conf_null=hist3([local_act(:),local_pred_null(:)],...
    'edges',{1:nstim,1:nstim});
normC=histc(local_act(:), 1:nstim);
conf_norm=bsxfun(@rdivide,conf,normC);

%Reorder (model fit relative to medium, put medium in middle for confusion matrix):
rr=[1 3 2]; %Low medium high
conf_temp=zeros(nstim,nstim);
for i=1:nstim
    for j=1:nstim
        conf_temp(i,j)=conf_norm(rr(i),rr(j));
    end
end
conf_norm=conf_temp;
figure
imagesc(conf_norm*100)
set(gca,'xtick',1:3,'ytick',1:3,'yticklabel',{'High','Medium','Low'},...
    'xticklabel',{'High','Medium','Low'})
xlabel('Predicted')
ylabel('Actual')
colormap('bone')

conf_final=hist3([actual', predicted'],'edges',{1:nstim,1:nstim});
null_conf_all=zeros(nstim,nstim,size(local_pred,2));
for i=1:size(local_pred,2)
    null_conf_all(:,:,i)=hist3([actual', local_pred_null(:,i)],'edges',{1:nstim,1:nstim});
end

for i=1:nstim
    for j=1:nstim
        crit=conf_final(j,i);
        null_dist=null_conf_all(j,i,:);
        null_dist(1)=crit;
        p(j,i)=sum(null_dist >=crit)/length(null_dist);
        text(i,j,sprintf('%2.1f%%',conf_norm(j,i)*100),'HorizontalAlign','Center','Color',[1 1 1])
    end
end
hold on

colormap(colorblind_cm)
colorbar
bi_Plot_Corrections

% Permutation test on significance of a confusion matrix
local_act=bsxfun(@times,ones(size(local_pred)),actual');
conf_total=hist3([local_act(:),local_pred(:)],...
    'edges',{1:nstim,1:nstim});
crit=multiclass_mcc(conf_total); %observed value
mcc_null(1)=crit;
set(gca,'clim',[10 62])
for i=2:size(local_pred,2)
    p=local_pred_null(:,i);
    act=actual';
    conf=hist3([act,p],'edges',{1:nstim,1:nstim});
    mcc_null(i)=multiclass_mcc(conf);
end
p_value=sum(mcc_null>=crit)/size(local_pred,2);

title(sprintf('%2.1f%% Correct, MCC=%1.2f, p=%1.3f',...
    sum(actual==predicted)/length(predicted)* 100,...
    crit,p_value))

%% Figure 7D: Plot Betas
preds={'Int','B','FR','IR','NSYN','%SYN'};
figure
if useInts==1
    ylab=preds(2:end);
    for i=1:length(ints)
        newLab{i}=sprintf('%s x %s',ylab{ints(i,1)},ylab{ints(i,2)});
    end
    allLab=[ {'int'} ylab newLab];
else
    allLab=predictors;
end


if nstim==3
    cats={'High','Low'};
    rel='Medium';
else
    cats={'High','Medium','Low'};
    rel='Med-High';
end

r=1;c=nstim-1;
for ii=1:size(beta,2)
    pk=median(beta_p(:,ii,:),3); %Very skewed- use median?
    b=ci95mean(squeeze(beta(:,ii,:)),2,0);
    if ii==1
        b=-b;
    end
    subplot(r,c,ii)
    hold on
    for i=1:size(b,1)
        plot(b(i,1),i,'ok','UserData','Ex')
        if pk(i) < 0.001
            h=text(b(i,1),i - (.2),'***');
        elseif pk(i) < 0.01
            h=text(b(i,1),i - (.2),'**');
        elseif pk(i) < 0.05
            h=text(b(i,1),i - (.2),'*');
        else
            h=text(b(i,1),i - (.2),'');
        end
        set(h,'HorizontalAlignment','Center')
        plot([0 b(i,1)],[i i],'k')
    end
    set(gca,'ydir','reverse','ytick',1:size(b,1),...
        'yticklabel',allLab)
    ylim([0,size(b,1) + 0.5])
    plot([0 0],[1, size(beta,1)+0.5],'--k')
    xlabel('Coefficients')
    if ii==1
        title(sprintf('%s coefficients relative to %s',rel,cats{ii}))
    else
        title(sprintf('%s coefficients relative to %s',cats{ii},rel))
    end
    xlim([-0.7 0.7])
end
set(gcf,'pos',[    680   557   835   421])
bi_Plot_Corrections

%% Supp figure 1D: Evaluate the fraction of mice in each depletion model where mean boot-strapped performance exceeded mean null performance
local_act=bsxfun(@times,ones(size(local_pred)),actual');
overall=sum(local_pred==local_act,2)./size(local_pred,2) * 100; % percent of iterations with correct predictions
overall_null=sum(local_pred_null==local_act,2)./size(local_pred_null,2) * 100;

asyn=ismember(origPerMouse,[6 7 8])';
grad=ismember(origPerMouse,[1 2 3 4])';
ctrl=ismember(origPerMouse,13)';
ba=ismember(origPerMouse,5)';
ipsi=ismember(origPerMouse,[9 10])';
labs={'Ctrl','Asyn','Grad','BA','Ipsi'};

inds={ctrl asyn grad ba ipsi};
cor=actual' == predicted';
clear obs exp tot dd mm prop_cor prop_chance
figure
hold on
for i=1:length(inds)
    ind=inds{i};
    tot(i)=sum(ind);
    dd{i}=overall(inds{i}); %Mean performance over iterations approach
    mm(i)=mean(dd{i});
    exp(i)= sum((overall_null(inds{i}) >  1/length(useResponse)*100)); %Number of mice where boot-strapped performance exceeded chance by chance
    prop_cor(i)=sum(dd{i}> 1/length(useResponse)*100); %Number of mice where boot-strapped performance exceeded chance
    text(i,90,sprintf('%d/%d',prop_cor(i),tot(i)),'horizontalalign','center')  % %of iterations that the correct outcome was chosen
end
[p,c]=chi2_test_oe(prop_cor,exp,[]);
bar(1:length(inds),mm,'k')

h=plotSpread(dd,[],[],labs);
set(h{1},'Marker','s','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7])
plot([0 length(inds)+1],1/length(useResponse)*100.*[1 1],'--r')
title(sprintf('Chi^2 = %2.1f, p = %1.4f',c,p))
set(gca,'ytick',0:25:100,'ylim',[0 100])
set(gca,'xtick',1:length(inds),'xticklabel',labs)
ylabel('% Correct')

%% Load model for Supp. 1C (4 Group confusion matrix)
load([pn_local 'mnrfit_500_perm_jackknife_100-75-35-10_ints_pred_phys_2ndmed_rel.mat'])

%Plot confusion matrix from all iteraions? Polot conusion matrix mode
local_act=bsxfun(@times,ones(size(local_pred)),actual');
conf=hist3([local_act(:),local_pred(:)],...
    'edges',{1:nstim,1:nstim}); %Confusion matrix performed on all iterations
conf_null=hist3([local_act(:),local_pred_null(:)],...
    'edges',{1:nstim,1:nstim});
normC=histc(local_act(:), 1:nstim);
conf_norm=bsxfun(@rdivide,conf,normC);

%Reorder (model fit relative to medium, put medium in middle for confusion matrix):
rr=[1 2 4 3]; %Low medium high
conf_temp=zeros(nstim,nstim);
for i=1:nstim
    for j=1:nstim
        conf_temp(i,j)=conf_norm(rr(i),rr(j));
    end
end
conf_norm=conf_temp;
figure
imagesc(conf_norm*100)
set(gca,'xtick',1:nstim,'ytick',1:nstim,'yticklabel',{'High','Medium','Med-Low','Low'},...
    'xticklabel',{'High','Medium','Med-Low','Low'})
xlabel('Predicted')
ylabel('Actual')
% colormap('bone')

conf_final=hist3([actual', predicted'],'edges',{1:nstim,1:nstim});
null_conf_all=zeros(nstim,nstim,size(local_pred,2));
for i=1:size(local_pred,2)
    null_conf_all(:,:,i)=hist3([actual', local_pred_null(:,i)],'edges',{1:nstim,1:nstim});
end

for i=1:nstim
    for j=1:nstim
        crit=conf_final(j,i);
        null_dist=null_conf_all(j,i,:);
        null_dist(1)=crit;
        p(j,i)=sum(null_dist >=crit)/length(null_dist);
        text(i,j,sprintf('%2.1f%%',conf_norm(j,i)*100),'HorizontalAlign','Center','Color',[1 1 1])
    end
end
hold on

colormap(colorblind_cm)
set(gca,'clim',[0 60])
colorbar
bi_Plot_Corrections

% Permutation test on significance of a confusion matrix
local_act=bsxfun(@times,ones(size(local_pred)),actual');
conf_total=hist3([local_act(:),local_pred(:)],...
    'edges',{1:nstim,1:nstim});
crit=multiclass_mcc(conf_total); %observed value
mcc_null(1)=crit;
for i=2:size(local_pred,2)
    p=local_pred_null(:,i);
    act=actual';
    conf=hist3([act,p],'edges',{1:nstim,1:nstim});
    mcc_null(i)=multiclass_mcc(conf);
end
p_value=sum(mcc_null>=crit)/size(local_pred,2);

title(sprintf('%2.1f%% Correct, MCC=%1.2f, p=%1.3f',...
    sum(actual==predicted)/length(predicted)* 100,...
    crit,p_value))

%% Figure 7E:

th=cell2mat(d(2:end,2));
phys1=cell2mat(d(2:end,4));
bilat_or_uni={'G85','G60','G30','G05','BA','A30','A60','A85','Ctrl','UDep','AsymDep'}; %'AsymDep'
figure
ind=ismember(d(2:end,3),bilat_or_uni);
physn=-phys1(ind);
thn=th(ind);
labs=d(2:end,3);
labs=labs(ind);

r=3; c=2;

fr=cell2mat(d(2:end,10));
sb=cell2mat(d(2:end,9));
cv=cell2mat(d(2:end,11));
syn=cell2mat(d(2:end,13)); %updated
nsync=cell2mat(d(2:end,12));
%
% raw={fr syn nsync cv sb};
% xlab={'FR','SYN','NSYNC','IR','B'};
% useLog=[0 1 0 1 1];


raw={fr syn cv sb};
xlab={'FR','SYN','IR','B'};
useLog=[0 1 1 1];

x=0:100;

saveFits=[];
useZ=1;
usePC=1:length(xlab); %i
figure
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.001; % smooth
for i=1:length(raw)
    subplot(r,c,i)
    if useLog(i)==1
        y=log10(raw{i}(ind));
        y(isinf(y))=0;
    else
        y=raw{i}(ind);
    end
    if useZ==1
        y=nanzscore(y,[],1);
    end
    if i==1
        fr_i=y;
    elseif i==2
        cv_i=y;
    elseif i==3
        sb_i=y;
    elseif i==4
        syn_i=y;
    end
    [xData, yData] = prepareCurveData( th(ind), y );
    [fitresult, gof1] = fit( xData, yData, ft, opts );
    plot(xData,yData,'ok','UserData','Ex')
    hold on
    plot(x,fitresult(x),'r')
    saveFits(:,i)=fitresult(x);
    if useZ==1
        ylabel(['z' xlab{i}])
    else
        ylabel(xlab{i})
    end
    title(sprintf('Adj r-sqr: %1.2f',gof1.adjrsquare))
    set(gca,'xdir','reverse','ylim',[-2.5 2.5],'xtick',0:25:100)
end
xlabel('%TH remaining')

figure
r=5; c=2;
% Summary
dat=saveFits;
endStage=dat(1:10,usePC);
earlyStage=dat(end-5:end,usePC);
cc2=zeros(size(dat));
pad=ones(30,length(usePC));
padDat=[bsxfun(@times,pad,dat(1,:));...
    dat(:,usePC);...
    bsxfun(@times,pad,dat(end-1,:))];
tt=size(pad,1)- 2;%round(size(endStage,1)/2)
clear temp2
for i=1:length(usePC)
    temp=xcorr(padDat(:,i),endStage(:,i),'none');
    temp2(:,i)=temp(size(padDat,1):end);
    cc2(:,i)=temp2(tt+1:(tt+size(dat,1)),i);
end

cc2=bsxfun(@minus,cc2,min(cc2,[],1));
cc2=bsxfun(@rdivide,cc2,max((cc2),[],1));
cc2=2*(cc2-0.5);

if length(raw)==5
    ccE=xcorr2([bsxfun(@times,pad(:,2:end),dat(1,2:end));...
        dat(:,2:end);...
        bsxfun(@times,pad(:,2:end),dat(end,2:end))],earlyStage(:,2:end));
else
    ccE=xcorr2([bsxfun(@times,pad(:,1:end),dat(1,1:end));...
        dat(:,1:end);...
        bsxfun(@times,pad(:,1:end),dat(end,1:end))],endStage);
end
ccE=ccE((tt+1:(tt+size(dat,1))),:);
ccE=-ccE;

subplot(r,c,[5,6])
imagesc(x,usePC,cc2')

set(gca,'xdir','reverse','xtick',0:25:100,'ytick',1:5,'yticklabel',xlab)
title('Similarity to End Stage')
set(gca,'clim',[-1 1].*1)
colorbar
subplot(r,c,[7,8])
imagesc(x,usePC,cc2')
set(gca,'xdir','reverse','xtick',0:25:100,'ytick',1:5,'yticklabel',xlab)
title('Similarity to End Stage')
colormap(colorblind_cm)

set(gca,'clim',[-1 1].*1)
xlabel('%TH Remaining')

ylabel('z')
set(gcf,'pos',[   680   249   353   729])
bi_Plot_Corrections

if reProc==1
    %% Fig 7F: Generate PCA identically to figure 3 with expanded data
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
    display('Finished')
    
    
    %% Compile all PC scores for summary
    disp('Performing PCA')
    abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
        'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
    usePredictors=[1,4,6,15]; %What if include neuron-specific synchrony index
    xlab={'%SB','FR','CV','%Sync'};
    params.phys_predict_cols=usePredictors;
    params.phys_predict_labs=xlab;
    
    grad=ismember(types,[1:4,13]);
    grad_ba=ismember(types,[1:5,13]);
    asyn=ismember(types,[6:8,13]);
    ind=ismember(types,[1:13]); %All conditions or use subset
    
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
    xlab={'%SB','FR','CV','Synch'};
    t=types(~xx & ind);
    params.types_after_indexing=unique(t);
    useData=dat;
    params.phys_interactions=0;
    useInts=1; %Use interactions for physiology data?
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
        %     p_dur(j,:)=mode(durs(mouse,1));
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
    
    fnSave='2019_01_18_V1';
    
    mkdir([pn_local 'output'])
    xlswrite([pn_local 'output' filesep fnSave '.xlsx'],allDat)
    behaveCoeff=wcoeff;
    mkdir([pn_local 'output'])
    save([pn_local 'output' filesep fnSave '.mat'],'allDat','physCoeff',...
        'behavCoeff','usePredictors','pn_local','fn','fn2','data_mouse',...
        'physUseData','params',...
        'useTH','explained_behav','explained_phys','ints','useInts')
    display('Finished')
end
%% Load data for Figure 7F
fn='2019_01_18_V1.mat'; %See
ddd=load([pn_local 'output' filesep fn]);
d=ddd.allDat;
display('Finished loading.')

%% Figure 7F  Plot Phys. vs. TH. vs Behavior (All)
use={'G85','G60','G30','G05','BA','A30','A60','A85','Ctrl'};
th=cell2mat(d(2:end,2));
type=d(2:end,3);
phys1=cell2mat(d(2:end,4));
behav=cell2mat(d(2:end,7));
figure
ind=ismember(type,use);
physn=phys1(ind,1);
bhvn=behav(ind,1);
zeroBehave=median(bhvn(ismember(type(ind),'Ctrl')));
thn=th(ind);
[xData, yData, zData] = prepareSurfaceData( thn,physn, bhvn);
ft = fittype( 'poly23' );

% Fit model to data.
[fitresult1, gof1] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );%'Robust','Bisquare'

% Create a figure for the plots.
n1=100;
n2=100;
dat1=linspace(-5,105,n1);
dat2=linspace(-1.4,1.6,n2);
tempZ=zeros(n1,n2);
for i=1:n1
    for j=1:n2
        tempZ(i,j)=fitresult1([dat1(j),dat2(i)]);
    end
end
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
title(sprintf('Adj-R^2 = %1.3f',gof1.adjrsquare))
[xData, yData] = prepareCurveData( thn, physn );
ft = fittype( 'poly3' );
[fitresult, gof] = fit( xData, yData, ft );
hold on
plot(0:100,fitresult(0:100),'--k')
set(gca,'xtick',0:25:100)
bi_Plot_Corrections

