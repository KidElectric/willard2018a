function aw2018_figure_1(pn_local,pn_common,reProc)
% function aw2018_figure_1(pn)
%
%   pn_local = data path for figure. if empty, assumes data are in a folder 'data' such that:
%   current directory is 'code,' and 'data' and 'code' are in 'fig' folder
%
%   pn_local = common data folder.  if empty, assumes data are in a folder 'common_data' such that:
%   current directory is 'code,' and 'common_data' is 2 levels above 'code'
%
%   reProc = rerun processing from common_data. code needs a few manual changes
%   in model fitting section to reproduce all models (described below).
%   default = 0 (may take over an hour)
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
%
if reProc==1
    %% Read AW snr data .xlsx
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
    useResponse=data(:,9); % Velocity
    xlab={'%SB','FR','CV','Synch'};
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
    
    params.useSeparateHems_As_Mice=0;
    if params.useSeparateHems_As_Mice==1
        for i=1:length(mouseID)
            mouseID{i}=[mouseID{i} skeep{i}];
        end
    end
    
    th_values=data2;
    display('Finished')
    
    %% Categorical multinomial for 2 groups: control vs graduals, endstage
    abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
        'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
    % abrev2={'Int','%SB','FR','CV','SYN'};
    useResponse=[1 2]; %For example, discriminate between Gradual 5% and bilateral accute
    
    typeTemp=zeros(size(types));
    origResp=useResponse;
    
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
    
    %Uncomment this code to discriminate gradual 5% from acute 5%: mnrfit_500_perm_jackknife_wta_grad5_vs_ba
    % typeTemp(types==4)=1; %Gradual 5% 
    % typeTemp(types==5)=2; %Bilateral Acute
    %Uncomment this code to discriminate gradual 5% + acute 5% from ctl: mnrfit_500_perm_jackknife_wta_depgt20_vs_ctrl
    typeTemp(types==13)=2; %Ctrl
    typeTemp(ismember(types,[4,5]))=1; %Gradual 5% OR bilateral acute
    
    %Note: change save filename below to avoid over-writing
    
    typeTemp(~type_include)=0; %Exclude unilaterals for now
    
    um=unique(mouseID);
    for i=1:length(um)
        tt(i)=length(unique(typeTemp(ismember(mouseID,um{i}))));
    end
    
    usePredictors=[1, 4, 6, 15]; %Synchrony index = ; Percent synchronous pairs = 15
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
    
    isSkew=[1 0 1 1];
    zScore=[1 1 1 1];
    
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
    fitType='mnrfit';
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
        save([pn_local , sprintf('mnrfit_%d_perm_jackknife_wta_depgt20_vs_ctrl_rerun.mat',perms)],...
            'pc','pc_null','conf_mat_norm_keep','beta','beta_p','beta_p_null','beta_null','conf_mat_norm_total',...
            'responses','useResponse','origResp','origPerMouse','local_pred','local_pred_null','nstim',...
            'predictors','sheets','cols','pnfn','fitType','conf_null_mat_norm_total',...
            'predicted','actual','ints','useInts','predicted_null','params')
        disp('Saved')
    end
end
%% Figure 1J & 1K -- compare classifications
%'mnrfit_500_perm_jackknife_wta_grad5_vs_ba.mat'
fns={'mnrfit_500_perm_jackknife_wta_grad5_vs_ba.mat',...
     'mnrfit_500_perm_jackknife_wta_depgt20_vs_ctrl.mat'};
labs={'G5% vs A5%','5% vs Ctrl'};
for ii=1:length(fns)
    load([pn_local fns{ii}],'pc','pc_null','conf_mat_norm_keep','beta','beta_p',...
        'beta_p_null','beta_null','conf_mat_norm_total',...
        'responses','useResponse','origResp','origRespLabels',...
        'combine','origPerMouse','predictors','sheets','cols',...
        'fitType','conf_null_mat_norm_total','predicted','actual',...
        'ints','useInts','predicted_null','local_pred','local_pred_null');
    order = [ 1 2];

    % Single box plot with all correct mice  
    figure    
    hold on
    exOut=[]; isPaired=1; anovaFirst=0; mc=0; useWilcox=0; negError=1; plotMedians=0;
    [h,pAll,mseO]=barMeanSig2({pc,pc_null},{'Actual','Chance'},{[1 2]},{'k','b'},...
        exOut,isPaired,anovaFirst,mc,useWilcox,negError,plotMedians);
    chance(ii,:)=stdErrMean(pc_null,2,1);
    perf(ii,:)=stdErrMean(pc,2,1);
   
    h=plotSpread(pc',[],[],labs(ii));
    set(h{1},'Marker','s','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7])
    set(gca,'ytick',0:25:100)
    xlim([0.5 2.5])
    set(gcf,'pos',[   680   482   297   496])
    ylabel('Accuracy (%)')
    bi_Plot_Corrections
end
%% For text: Evaluate the fraction of mice in each depletion model where mean boot-strapped performance exceeded mean null performance
local_act=bsxfun(@times,ones(size(local_pred)),actual');
overall=sum(local_pred==local_act,2)./size(local_pred,2) * 100; % percent of iterations with correct predictions
overall_null=sum(local_pred_null==local_act,2)./size(local_pred_null,2) * 100;
useResponse=[1 2];
grad=ismember(origPerMouse,[1 2 3 4])';
ba=ismember(origPerMouse,5)';
ctrl=ismember(origPerMouse,13)';
labs={'G5%', 'BA', 'Ctrl'};
inds={grad ba ctrl};
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
bar(1:length(inds),mm,'k')
h=plotSpread(dd,[],[],labs);
set(h{1},'Marker','s','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7])
plot([0 length(inds)+1],1/length(useResponse)*100.*[1 1],'--r')
set(gca,'ytick',0:25:100,'ylim',[0 100])
set(gca,'xtick',1:length(inds),'xticklabel',labs)
ylabel('% Correct')
%% Figure 1L - Plot betas
load([pn_local fns{2}],'beta','beta_p');
preds={'Int','B','FR','IR','SYN'};
reOrder=[3 4 2 5]-1;
if useInts==1
    ylab=preds(2:end);
    for i=1:length(ints)
        newLab{i}=sprintf('%s x %s',ylab{ints(i,1)},ylab{ints(i,2)});
    end
    allLab=[ {'int'} ylab(reOrder) newLab];
%     allLab=[ {'int'} ylab newLab];
else
    allLab=predictors;
end
reOrder=[1 3 4 2 5];
% reOrder=1:5;
beta_temp=beta;
% beta=beta*2; %because logistic, *2
beta_temp(1:5,:,:)=-beta_temp(reOrder,:,:); %if relative to ES flip sign
beta_temp=logmodulus(beta_temp);
beta_p(1:5,:,:)=beta_p(reOrder,:,:);
for ii=1:size(beta_temp,2)
    pk=median(beta_p(:,ii,:),3); %Very skewed- use median?
    b=ci95mean(squeeze(beta_temp(:,ii,:)),2,0);
    figure
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
    plot([0 0],[1, size(beta_temp,1)+0.5],'--k')
    xlabel('logmod(beta)')
    set(gcf,'pos',[ 680   616   267   362])
    bi_Plot_Corrections

end
   set(gcf,'pos',[680   616   238   362])