function  aw2018_figure_6(pn_local,pn_common,reProc)
% function aw2018_figure_6(pn_local,pn_common,reProc)
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
    
    %% Categorical multinomial for 2 groups: control different types of depleted
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
    
    % type_include=ismember(types,[4,5,9,10]);% Asym and uni IPSI included
    % type_include=ismember(types,[4,5,13]);% Asym and uni IPSI included
    type_include=ismember(types,[4,5,9,10,13]);% Asym and uni IPSI included
    %Predict quantile behavior
    
    % Figure 6N:
    % typeTemp(types==4)=1;
    % typeTemp(types==5)=2;
    % typeTemp(types==9)=3;
    % typeTemp(types==10)=4;
    
    % Figure 6O
    % typeTemp(ismember(types,[4,5]))=1;
    % typeTemp(types==13)=2;
    
    % Figure 6P:
    typeTemp(ismember(types,[4,5,9,10]))=1;
    typeTemp(types==13)=2;
    
    typeTemp(~type_include)=0;
    
    ut=unique(typeTemp);
    useResponse=ut(ut~=0);
    
    
    % typeTemp(nn < 10)=0; %Exclude low-neuron mice
    um=unique(mouseID);
    for i=1:length(um)
        tt(i)=length(unique(typeTemp(ismember(mouseID,um{i}))));
    end
    
    % origRespLabels=abrev(origResp);
    
    usePredictors=[1, 4, 6, 15]; %Synchrony index = 8; Percent synchronous pairs = 15
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
            isSkew=[1 0 1 1];
            zScore=[1 1 1 1];
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
    n=length(uInputs);
    pk=[];
    type_names=abrev(useResponse);
    targetN=50;
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
        save([pn_local , sprintf('mnrfit_%d_perm_dep_ES_Vs_ctrl_pred_each.mat',perms)],...
            'pc','pc_null','conf_mat_norm_keep','beta','beta_p','beta_p_null','beta_null','conf_mat_norm_total',...
            'responses','useResponse','origPerMouse','local_pred','local_pred_null','nstim',...
            'predictors','sheets','cols','pnfn','fitType','conf_null_mat_norm_total',...
            'predicted','actual','ints','useInts','predicted_null','params')
        disp(sprintf('mnrfit_%d_perm_dep_ES_Vs_ctrl_pred_each.mat -Saved!',perms))
    end
    
    figure
    imagesc(conf_mat_norm_total)
    
end
%% Figure 6 Unilateral model comprison plots K,L,M

fns={'mnrfit_500_perm_dep_4way_ints_pred_each.mat',...
    'mnrfit_500_perm_dep_ipsi_Vs_ctrl_pred_each.mat',...
    'mnrfit_500_perm_dep_ES_Vs_ctrl_pred_each.mat'};
clear dd n
cc=[0.25 0.5 0.5].*100;
for ii=1:length(fns)
    load([pn_local fns{ii}]);
    % Plot each mouse in each class vs chance
    ct=unique(actual);
    ci_null=ci95mean(pc_null,2,0);
    alpha=0.05;
    n(ii)=length(actual);
    [h,p(ii)]=ttest(pc ,cc(ii),'tail','right');
    chance(ii,:)=stdErrMean(pc_null,2,1);
    perf(ii,:)=stdErrMean(pc,2,0);
    dd{ii}=pc;
    cor(ii)=sum(pc>cc(ii));
    tot(ii)=length(actual);
end
% Single box plot with all correct mice
figure
hold on
bar(1:3,perf(:,1),'b')
h=plotSpread(dd,[],[],{'4-way','Ipsi vs Ctl','ES vs Ctl'});
for i=1:length(fns)
    plot([i i],perf(i,2:3),'k')
end

set(h{1},'Marker','s','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7])
set(gca,'ytick',0:25:100)
set(gcf,'pos',[ 557   158   431   364])
ylabel('Accuracy (%)')
bi_Plot_Corrections()
%% Figure 6 ES depletion model performance comparison

load([pn_local 'mnrfit_500_perm_dep_ES_Vs_ctrl_pred_each.mat']);
act_local=bsxfun(@times,ones(size(local_pred)),actual');
overall=sum(local_pred==act_local,2)./size(local_pred,2) * 100;
overall_null=sum(local_pred_null==act_local,2)./size(local_pred_null,2) * 100;

uni=ismember(origPerMouse,[9 10])';
grad=ismember(origPerMouse,[1 2 3 4])';
ctrl=ismember(origPerMouse,13)';
ba=ismember(origPerMouse,5)';
labs={'Ctrl','Ipsi','Grad','BA'};

inds={ctrl uni grad ba};
cor=actual' == predicted';
clear obs exp tot dd mm prop_cor prop_chance pt
figure
hold on
for i=1:length(inds)
    ind=inds{i};
    tot(i)=sum(ind);
    dd{i}=overall(inds{i}); %Mean performance over iterations approach
    mm(i)=mean(dd{i});
    exp(i)= sum((overall_null(inds{i}) >  1/length(useResponse)*100));
    prop_cor(i)=sum(dd{i}> 1/length(useResponse)*100); %Number of mice where boot-strapped performance exceeded chance
    h=text(i,90,sprintf('%d/%d',prop_cor(i),tot(i)),'horizontalalign','center');  % %of iterations that the correct outcome was chosen
    if i==4
        set(h,'Color','w')
    end
end
bar(1:length(inds),mm,'k')

h=plotSpread(dd,[],[],labs);
set(h{1},'Marker','s','MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7])
plot([0 length(inds)+1],(1/length(useResponse))*100.*[1 1],'--r')
set(gca,'ytick',0:25:100,'ylim',[0 100])
set(gca,'xtick',1:length(inds),'xticklabel',labs)
ylabel('% Correct')
bi_Plot_Corrections()