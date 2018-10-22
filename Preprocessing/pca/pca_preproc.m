function pca_preproc(pn)
% function pca_preproc(pn)
%   pn = data path. if empty, assumes data are in a folder 'data' such that:
%   current directory is 'code,' and 'data' and 'code' are in same folder
% Load physiology data
% Read AW snr data .xlsx
if ~exist('pn','var') || isempty(pn)
    pn=['.' filesep '..' filesep 'data' filesep];
end
disp('Loading data...')
fn='SNrInVivoData_wRest_FINAL_v7.xlsx';
fn2='SynchronyData_6-6-18_useable_updated.xlsx';
fn3='GradualProtocolDuration.xlsx';

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



%
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
% mouseID0=mouseID;

%%%%
useSeparateCond_As_Mice=1;
if useSeparateCond_As_Mice==1
   for i=1:length(mouseID)
      mouseID{i}=[mouseID{i} sprintf('_%d',types(i))]; 
   end
end

useResponse=data(:,9); % Velocity
xlab={'%SB','FR','CV','Synch'};
% [1,4,6,15];
%include vel as variable for PCA:
xlab=[xlab 'Vel'];
% t=types(~xx & ind);

% Load in raw %TH values

%All Categories:
% sheets={'Grad85','Grad60','Grad30','Gradual5'...
%     'Acute',...
%     'ASyn30','ASyn60','ASyn85',...
%     'UniDepl','1injIntactAsym','2injIntactAsym','5injIntactAsym'...
%     '1injDeplAsym','2injDeplAsym','5injDeplAsym','UniCtl','Ctl'};
% 
% rows={'G3:I274','G3:I320','G3:I320','G3:I295',...
%     'F3:H247',...
%     'G3:I152','G3:I157','G3:I100',...
%     'G3:I138','G3:I112','G3:I136','G3:I86'...
%     'G3:I175','G3:I177','G3:I114','G3:I96','G3:I264','none'};
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
display('Finished')


%% Compile all PC scores for summary
% pn='C:\Users\Gittis\Dropbox\Gittis Lab Data\AWSNR\';
disp('Performing PCA')
abrev={'G85','G60','G30','G05','BA','A30','A60','A85',...
    'UDep','AsymDep','AsymCtrl','UCtrl','Ctrl'};
usePredictors=[1,4,6,15];

grad=ismember(types,[1:4,13]);
grad_ba=ismember(types,[1:5,13]);
asyn=ismember(types,[6:8,13]);
ind=ismember(types,[1:13]); %All conditions or use subset
xx=any(isinf(data(:,[usePredictors,9])) | isnan(data(:,[usePredictors,9])),2); %Exclude bad data and others
dat=data(~xx & ind,usePredictors);
tempTH=th_values;
useTH=th_values(~xx & ind);
isSkew=[1 0 1 1];
datRaw=dat;
for i=1:length(isSkew)
    if isSkew(i)
        dat(:,i)=log10(dat(:,i));
        dat(isinf(dat(:,i)),i)=0;
    end
end
xlab={'%SB','FR','CV','Synch'};
t=types(~xx & ind);
useData=dat;
useInts=1; %Use interactions for physiology data?
if useInts==1
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
useData=data(~xx & ind,behave);
useMouse=mouseID(~xx & ind);
durs=data(~xx & ind,16);
t=types(~xx & ind);
useLog=[1 0 1 1];
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
    for i=1:4
        if isSkew(i)
            data_mouse_raw(j,i)=nanmedian(datRaw(mouse,i),1);
        else
            data_mouse_raw(j,i)=nanmean(datRaw(mouse,i),1);
        end
    end
    keep_mouse{j}=a{j};
    keep_score(j,:)=nanmean(physscore(mouse,:),1);
    p_dur(j,:)=mode(durs(mouse,1));
    th_mouse(j,1)=nanmean(useTH(mouse));
end

sign_invert=0; %Depends on whether sign convention assigned to chosen
%data subset matches sign convention of the paper
if sign_invert==1
    keep_score=-keep_score;
    physCoeff=-physCoeff;
end

useInts=0; %Decide whether interactions used for behavior data (typically not)
if useInts==1
    ints=combntns(1:(size(data_mouse,2)),2);
    for i=1:size(ints,1)
        data_mouse=[data_mouse data_mouse(:,ints(i,1)).*data_mouse(:,ints(i,2))];
    end
end

w = 1./var(data_mouse);
[wcoeff,score,latent,tsquared,explained_behav] = pca(data_mouse,...
    'VariableWeights',w);
behavCoeff=inv(diag(std(data_mouse)))*wcoeff;

allDat=[keep_mouse',num2cell(th_mouse),abrev(type_mouse)',...
    num2cell(keep_score(:,1)), num2cell(keep_score(:,2)),...
    num2cell(keep_score(:,3)),...
    num2cell(score(:,1)), num2cell(score(:,2)),...
    num2cell(data_mouse_raw),num2cell(data_mouse),...
    num2cell(p_dur)]; 
cols=[cellstr('Mouse'),cellstr('%TH'), cellstr('Type'), ...
    cellstr('Phys PC1 Score'), ...
    cellstr('Phys PC2 Score'),cellstr('Phys PC3 Score'), ...
    cellstr('Behav PC1 Score'),cellstr('Behav PC2 Score'),...
    xlab,'Vel','Rear','TotePole','WireHang','Protocol Duration'];
allDat=[cols; allDat];
% fnSave='4_behavs_all_states+raw_no_ints_2';
fnSave='4_behavs_all_states+all_ints_5';

mkdir([pn 'output'])
xlswrite([pn 'output' filesep fnSave '.xlsx'],allDat)
behaveCoeff=wcoeff;
mkdir([pn 'output'])
save([pn 'output' filesep fnSave '.mat'],'allDat','physCoeff',...
    'behavCoeff','usePredictors','pn','fn','fn2','fn3','data_mouse',...
    'physUseData',...
    'useTH','explained_behav','explained_phys','ints','useInts')
display('Finished')


