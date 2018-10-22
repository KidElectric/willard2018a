function cv_model(pn)
% function cv_model(pn)
%   pn = data path. if empty, assumes data are in a folder 'data' such that:
%   current directory is 'code,' and 'data' and 'code' are in 'fig' folder
% Load physiology data
% Read AW snr data .xlsx
if ~exist('pn','var') || isempty(pn)
    pn=['.' filesep '..' filesep 'data' filesep];
end
%
%% Read AW snr data .xlsx
fn='SNrInVivoData_wRest_FINAL_v7.xlsx';
fn2='SynchronyData_6-6-18_useable_updated.xlsx';
fn3='GradualProtocolDuration.xlsx';
disp('Loading data...')
%Protocol Duration:
[~,~,gradpdur]=xlsread([pn fn3],'Gradual','A2:D30');
[~,~,asynpdur]=xlsread([pn fn3],'Alpha-Syn','A2:D11');
dur_ID=[gradpdur(:,2) ; asynpdur(:,2)];
dur_val=cell2mat([gradpdur(:,4) ; asynpdur(:,4)]);

%New synch:
[~,~,psync]=xlsread([pn fn2],'Sheet1','A2:C65');
ps_ID=psync(:,2);
ps_val=cell2mat(psync(:,3));

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


%add percent synchronous pairs col
data=[data nan(size(data,1),1)];
for i=1:length(ps_ID)
    ind=ismember(mouseID,ps_ID{i});
    data(ind,15)=ps_val(i);
end

%add protocol duration column:
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

useSeparateHems_As_Mice=0;
if useSeparateHems_As_Mice==1
   for i=1:length(mouseID)
      mouseID{i}=[mouseID{i} skeep{i}]; 
   end
end

th_values=data2;
display('Finished loading data.')


%% Categorical multinomial: 1) Normalize neurons per mouse, 2) Jack-knife normalized mice
abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
    'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
abrev2={'Int','%SB','FR','CV','SYN'};
useResponse=[4 5]; %For example, discriminate between Gradual 5% and bilateral accute
typeTemp=types;
origResp=useResponse;

%Example for how to combine groups:
% combine=[9 10];
% typeTemp(ismember(types,combine))=9; %Combine Unilateral Dep and Asym Dep
% combine2=[4 5];
% typeTemp(ismember(types,combine2))=5; %Combine G5%, BA, 
% origResp=[combine combine2];
% useResponse=[5 9]; %New groups to compare
combine=[];
origRespLabels=abrev(origResp);
useTH=0;
usePredictors=[1,4,6,15]; %Percent synchronous pairs = 15
useCol=usePredictors;
useData=data(:,useCol); %Put Ctrl in right-most column to be "reference category"

%Include TH
if useTH==1
    useData=[useData, th_values];
end
originalTypes=types;
useTypes=typeTemp;
typeInd=ismember(useTypes,useResponse);
% useResponse=[5 9]; %For combining
useMID=mouseID;
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

clear neuronsPerMouse typePerMouse collectVals origPerMouse
for i=1:length(uInputs)
    ind=ismember(useMID,uInputs{i});    
    neuronsPerMouse(i)=sum(ind);
    collectVals(ind)=neuronsPerMouse(i);
    typePerMouse(i)=mean(useTypes(ind)); %all have same value   
    origPerMouse(i)=mean(originalTypes(ind));
end


if useTH==1
    isSkew=[1 0 1 1 0];
    zScore=[1 1 1 1 1];
    isCategorical=[0 0 0 0 0];
else
    isSkew=[1 0 1 1];
    zScore=[1 1 1 1];
    isCategorical=[0 0 0 0];
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
useCol=1:size(useData,2);


useInts=1; %Use interactions? 1 = yes, 0 = no
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
% fitType='glmnet';  %Harder to get running these days
fitType='mrfit';
options.alpha=1;
options.parallel=1;
n=length(uInputs);
pk=[];
type_names=abrev(useResponse);
targetN=round(median(neuronsPerMouse));
runNull=1;
for ii=1:n % For each mouse, drop it out, fit model on remaining mice etc.
    ex=uInputs{ii};    %Cell to exclude in jack-knife
    useMouse=find(~ismember(uInputs,ex));
    use=~ismember(useMID,ex);    
    clear pred_tot local_correct local_pred local_correct_null local_pred_null
    act=typePerMouse(ismember(uInputs,ex));
    bKeep=[]; 
    for i=1:perms
        tic
        % 1) Randomly downsample remaining cells to have equal examples of
        % each class        
        ec=equalize_class_quant(typePerMouse(useMouse),1:nstim);
        useMouseID=uInputs(useMouse(ec)); %Bingo.
%         figure;hist(typePerMouse(useMouse(ec))) % check equal

        % 2) Randomly resample puncta of each cell to targetPunct   
        cj=[];
        for j=1:length(useMouseID)
           ind=ismember(useMID,useMouseID{j});
           cc=find(ind); 
           if sum(ind) < targetN
               cj=[cj;datasample(cc,targetN,'Replace',true)]; %Up-sample
           else
               cj=[cj;datasample(cc,targetN,'Replace',false)]; %Down-sample
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
%             display('1')
            cvfit=cvglmnet(permPredictors,permResponse,'multinomial',options,[],3);
%             cvfit=glmnet(permDat,permClass,'multinomial',options);
            B0= cvglmnetCoef(cvfit,cvfit.lambda_1se);
%             B0=glmnetCoef(cvfit);
            B=[B0{:}];
            pred=cvglmnetPredict(cvfit,testDat,cvfit.lambda_1se,'response');
%             pred=glmnetPredict(cvfit,testDat,[],'response');
            pl=sum(pred,1);
        
            %Winner take all:
            pred=find(pl==max(pl));
%             pred_tot(i,:)=sum(cvglmnetPredict(cvfit,useData(~use,useCol),cvfit.lambda_1se,'response'),1);
            % Null model:
            if runNull==1
%                 display('2')
                permClassNull=datasample(permResponse,length(permResponse),'Replace',false);
                cvfit=cvglmnet(permPredictors,permClassNull,'multinomial',options,[],3);
                B0= cvglmnetCoef(cvfit,cvfit.lambda_1se);
                %             B0=glmnetPredict(lfit,testDat,[],'coefficient'); %cvfit.lambda_1se
                B_null=[B0{:}];
                pred_null=cvglmnetPredict(cvfit,testDat,cvfit.lambda_1se,'response');
                pl=sum(pred_null,1);
        
                %Winner take all:
                pred_null=find(pl==max(pl));
                if length(pred_null) > 1
                    pred_null=datasample(pred_null,1); %Pick randomly
                end
                %Winner take all:
                bKeepNull=cat(3,bKeepNull,B_null);
                
                local_correct_null(i)=pred_null==act;
                local_pred_null(i)=pred_null;
            end  
        else
            [B,dev,stats] = mnrfit(permPredictors,permResponse,'Interactions','on',...
                'model','nominal');
            pred=mnrval(B,testDat,stats);%'model','ordinal','ineractions','on');
            pl=sum(pred,1);
        
            %Winner take all:
            pred=find(pl==max(pl));
%             pKeep=cat(3,pKeep, stats.p);
            pKeep= stats.p;
            
            % Null model:
            if runNull==1
                permClassNull=permResponse(randperm(length(permResponse),length(permResponse)));                
                [B_null,dev_null,stats_null] = mnrfit(permPredictors,permClassNull,'Interactions','on',...
                    'model','nominal');
                pred_null=mnrval(B_null,testDat,stats_null);%'model','ordinal','ineractions','on');
                pl=sum(pred_null,1);

                %Winner take all:
                pred_null=find(pl==max(pl));
                pKeep_null=cat(3,pKeep_null, stats_null.p);
                bKeepNull=cat(3,bKeepNull,B_null);
                
                local_correct_null(i)=pred_null==act;
                local_pred_null(i)=pred_null;
            end           

        end
        
        local_correct(i)=pred==act;
        local_pred(i)=pred;
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
        fprintf('Jack-knife %d, Perm %d: ',ii,i)
        fprintf('%s Mouse predicted as %s.\t',type_names{act},type_names{pred})
        toc
    end
    pc(ii)=sum(local_correct)/perms * 100;
    fprintf('\n%2.1f%% correct, for mouse %d of %d.\n',pc(ii),ii,length(uInputs))

    conf_mat_norm_keep(:,:,ii)=mean(conf_mat_norm,3);
    bk(:,:,ii)=mean(bKeep,3);
    
    %Evaluate holdout
    actual(ii) = act;
    predicted(ii)= mode(local_pred);
    if runNull==1
        pc_null(ii)=sum(local_correct_null)/perms * 100;
        predicted_null(ii)=mode(local_pred_null);
        bkn(:,:,ii)=mean(bKeepNull,3);
        fprintf('%2.1f%% null correct, for cell %d of %d.\n',pc_null(ii),ii,length(uInputs))
    end
    name=ex;
    if strcmp(fitType,'glmnet')

    else
        pk(:,:,ii)=pKeep; %ALL
        if runNull==1
           pkn(:,:,ii)=median(pKeep_null,3); 
        end
    end
    fprintf('%s Cell predicted as %s.\n\n',type_names{act},type_names{predicted(ii)})
end

conf=hist3([actual' predicted'],'edges',{1:nstim,1:nstim});

%Normalize
normC=histc(actual,1:nstim)';%Normalization constants for each stimulus type (number of times each stimulus was used as a test trial)
conf_norm=bsxfun(@rdivide,conf,normC);
conf_norm(isinf(conf_norm))=nan;
conf_mat_norm_total=conf_norm;

if runNull==1
   conf_null=hist3([actual' predicted_null'],'edges',{1:nstim,1:nstim}); 
   normC=histc(actual,1:nstim)';%Normalization constants for each stimulus type (number of times each stimulus was used as a test trial)
   conf_norm=bsxfun(@rdivide,conf_null,normC);
   conf_norm(isinf(conf_norm))=nan;
   conf_null_mat_norm_total=conf_norm;
end

beta=bk;
beta_p=pk;
beta_p_null=pkn;
beta_null=bkn;
responses=type_names;
predictors=['int' cols];
if runNull==1
    save([pn , sprintf('mnrfit_%d_perm_jackknife_mice_g5_v_ba.mat',perms)],...
    'pc','pc_null','conf_mat_norm_keep','beta','beta_p','beta_p_null','beta_null','conf_mat_norm_total',...
    'responses','useResponse','origResp','origRespLabels','combine','origPerMouse',...
    'predictors','sheets','cols','pnfn','fitType','conf_null_mat_norm_total',...
    'predicted','actual','ints','useInts','predicted_null')
    display('Saved')
end

