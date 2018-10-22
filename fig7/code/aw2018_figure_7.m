function aw2018_figure_7(pn)
% function aw2018_figure_7(pn)
%   pn = data path. if empty, assumes data are in a folder 'data' such that:
%   current directory is 'code,' and 'data' and 'code' are in 'fig' folder
% Load physiology data
% Read AW snr data .xlsx
if ~exist('pn','var') || isempty(pn)
    pn=['.' filesep '..' filesep 'data' filesep];
end
%
%%
fn2='4_behavs_all_states+all_ints_4'; % used with 2 sample t-test on all hemispheres (

[~,~,d]=xlsread([pn fn2 '.xlsx']);
load([pn sprintf('%s.mat',fn2)],'physCoeff','explained_phys','ints')
% [~,~,d]=xlsread([pn '\MouseLevel\4_behavs_with_all_states+raw_no_ints_nlog_synch_1.xlsx']);
fr=cell2mat(d(2:end,10));
sb=cell2mat(d(2:end,9));
cv=cell2mat(d(2:end,11));
syn=cell2mat(d(2:end,12));

% [~,~,d]=xlsread([pn '\MouseLevel\mouse_3x2score_summary_nobehaveorphys_ints.xlsx']);
% phys=logmodulus(cell2mat(d(2:end,4:6)));
phys=(cell2mat(d(2:end,4:6)));

if strcmp(fn2,'4_behavs_all_states+raw_no_ints_2') ...
        || strcmp(fn2,'4_behavs_grad_ba_asyn+raw_no_ints_2')
    phys(:,1)=-phys(:,1);
    physCoeff(:,1)=-physCoeff(:,1);
end
% physCoeff=physCoeff;
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
% zbehave2=nanzscore(behav(:,2),[],1);


asyn=ismember(type,{'A30','A60','A85'});
grad=ismember(type,{'G30','G60','G85'});
ctl=ismember(type,{'Ctrl'});
ba=ismember(type,{'BA'});
isEndStage=ismember(type,{'BA','G05'});
allgrad=ismember(type,{'G05','G30','G60','G85'});

uniContra=ismember(type,'UCtrl');
uniIpsi=ismember(type,'UDep');
asymContra=ismember(type,{'1IIA','2IIA','5IIA','AsymCtrl'});
asymIpsi=ismember(type,{'1IDA','2IDA','5IDA','AsymDep'});

% Only use paired unilateral mice:
% 
uniContra=ismember(type,'UCtrl');
uniIpsi=ismember(type,'UDep');
anTemp=unique([anID(uniIpsi);anID(uniContra)]);
coan=anTemp(ismember(anTemp,anID(uniContra)) & ismember(anTemp,anID(uniIpsi)));
isPaired= ismember(anID,coan);
uniIpsi= uniIpsi & isPaired;
uniContra= uniContra & isPaired;

anTemp=unique([anID(asymIpsi);anID(asymContra)]);
coan=anTemp(ismember(anTemp,anID(asymContra)) & ismember(anTemp,anID(asymIpsi)));
isPaired= ismember(anID,coan);


display('Finished loading')
%% Figure 7A & 7C
ylab={'B','FR','IR','SYN'};
reOrder=[2 3 1 4]; 
useInts=1;
r=4;c=1;
d0={};
subinds=subplot_indexer(r,c);
flipSign=0; %of all
usePCs=1:3;
% ints=combntns(1:(size(useData,2)),2);
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
ind=asyn | allgrad | ba | ctl;
zscale=0.5;
lab={'Dep','30','60','85','Ctl'};
bins=[0 15 48 73 100];
useBins=[1 3 5];
lab=lab(useBins);
for i=usePCs       
    if flipSign==1
        out(:,i)=-out(:,i);
    end
    subplot(r,c,subinds(i,1))
    [~,ind1]=histc(th,bins);
    temp=[];
    for ii=1:length(useBins)
        d0=phys(ind1==useBins(ii) & ind,i);
        temp=[temp;mean(d0*out(:,i)',1)];
    end
    imagesc(fliplr(temp(:,reOrder)'))
    set(gca,'ytick',1:size(out,1),...
        'yticklabel',allLab(reOrder),'xtick',1:length(bins),...
        'xticklabel',fliplr(lab),'clim',[-1 1].*zscale)
end
 subplot(r,c,subinds(i+1,:))
     imagesc(temp(reOrder)')
    set(gca,'ytick',1:size(out,1),'xtick',usePCs,...
    'clim',[-1 1].*zscale)
colorbar
colormap colorblind_cm
bi_Plot_Corrections

set(gcf,'pos',[680   215   184   763])
 figure; 
 subplot(1,2,1)
 barh((physCoeff(reOrder,1)),'k')
 title('PC1')
set(gca,'ytick',1:size(out,1),'yticklabel',allLab(reOrder),'ydir','reverse')

 subplot(1,2,2)
 title('PC2')
 hold on
 barh((physCoeff(reOrder,2)),'k')
 set(gca,'ytick',1:size(out,1),'yticklabel',allLab(reOrder),'ydir','reverse')
%  bi_linkaxes('x')
 xlim([-1 1])
 bi_Plot_Corrections

%% Figure 7B & 7D: Plot combined ipsi and contra mice onto the physiology trends 
r=1;c=2;
x=0:100;
figure

ind= allgrad | asyn| ctl | ba;
fit_types={'poly3','poly1',};
Robust={'off','off',};
zd={phys(:,1), phys(:,2)};
clear p
for i=1:c   
    % NORMALIZED scores and fits:
    [xData, yData] = prepareCurveData( th(ind), zd{i}(ind) );
    ft = fittype(fit_types{i} );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Normalize = 'off';
    opts.Robust=Robust{i};
    [fitresult1, gof1] = fit( xData, yData, ft, opts );
    p1 = predint(fitresult1,x,0.95,'functional','off');
    [xData, yData] = prepareCurveData( th(ind), zbehave1(ind) );
    
    ft = fittype( 'poly2' );    
    [fitresult2, gof2] = fit( xData, yData, ft, opts );
    p2 = predint(fitresult2,x,0.95,'functional','off');
    set(gca,'Xdir','reverse','xtick',0:25:100)
    
    subplot(r,c,i)
    pcShift=fitresult1(100); %Value at 100% dopamine    
    behavShift=fitresult2(100);
    h=plot(x,fitresult1(x)-pcShift,'r');   
    
    hold on
    plot(x,p1-pcShift,'--r')
    h.DisplayName= 'Physiology PC1';
    ylabel(sprintf('Normalized Scores',i))
    
    %Add on the mean unilateral + asymmetric
    ucth=th(uniContra);
    udth=th(uniIpsi);
    ucz=zd{i}(uniContra)-pcShift;
    udz=zd{i}(uniIpsi)-pcShift;   
    uczb=zbehave1(uniContra)-behavShift;
    ascth=th(asymContra);
    asith=th(asymIpsi);    
    ascz=zd{i}(asymContra)-pcShift;
    asiz=zd{i}(asymIpsi)-pcShift;   
    asczb=zbehave1(asymContra)-behavShift;
    
    ub=mean([uczb; asczb]); 
    ube=stdErr([uczb; asczb],1);
    [h,p(i)]=ttest2([ascz; ucz],[asiz; udz]); %unpaired hems 

    xdist=120;
    tt=ci95mean([ascz; ucz],1,0);
    tt2=stdErr([ascz; ucz],1);
    plot(xdist,tt(1),'ok','UserData','Ex')

    plot([0 xdist xdist 0 0],...
        [ [1 1].*tt(1)+tt2 [1 1].*tt(1)-tt2 tt(1)+tt2] ,'--k')
    
    plot([1 xdist-1],[1 1].*mean([ascz; ucz]),'--k')

    tt=ci95mean([asiz; udz],1,0);
    tt2=stdErr([asiz; udz],1);

    tt=ci95mean([asiz; udz],1,0);
    plot(0,tt(1),'og','UserData','Ex')
    plot([1 xdist-1],[1 1].*mean([asiz; udz]),'--g')
    plot([0 xdist xdist 0 0],...
            [ [1 1].*tt(1)+tt2 [1 1].*tt(1)-tt2 tt(1)+tt2] ,'--g')
    plot([0 xdist], [ 0 0],'--k')
    xlabel('%TH Remaining')
    set(gca,'Xdir','reverse','xtick',0:25:100)
    ylim([-2.5 1.5])
end
% ylim([-2 2])
% bi_linkaxes('x')
bi_Plot_Corrections
xlim([0 130])
set(gcf,'pos',[   680   790   512   188])

%% Figure 7E,7F,7G

fns={'mnrfit_100_perm_jackknife_endstage_mice_acuteVg05VudepVasymdep.mat',...
    'mnrfit_100_perm_jackknife_endstage_mice_acute+g05+udep+asymdep_vs_ctl.mat',...
    'mnrfit_100_perm_jackknife_endstage_mice_udep+asymdep_v_ctl.mat'};
for ii=1:length(fns)
    load([pn fns{ii}]);
    % Plot each mouse in each class vs chance
    ct=unique(actual);
    ci_null=ci95mean(pc_null,2,0);
    clear d p dn
    alpha=0.05;
    nr=length(responses);
    figure;imagesc(1:nr,1:nr,conf_mat_norm_total);set(gca,'xtick',1:nr,'xticklabel',responses,'ytick',1:nr,'yticklabel',responses,'ydir','normal')
    ttt=nanmean(conf_mat_norm_keep,3);
    sss=stdErr(conf_mat_norm_keep,3);
    ipsiContraConfusion=( sum(ismember(actual,3) & predicted==4 )+ sum(ismember(actual,4) & predicted==3)) / sum(ismember(actual,[3,4]));
    for i=1:length(ct)
        ind=actual==ct(i);
        d{i}=pc(ind);
        dn{i}=pc_null(ind);
        [h,pt(ii,i)]=ttest(d{i},dn{i},alpha,'right');
        
    end
    [h,p_tot(ii)]=ttest(pc,pc_null,[],'right');
    figure
    isGrad=0;
    if isGrad==1
        order=[4 1 2 3];
        colCel={{[0 0 0],[243 243 243]./255},...
            {[0 0 0], [243 207 227]./255},...
            {[0 0 0], [230 164 201]./255},...
            {[0 0 0], [208 90 161]./255},...
            };
        xlab={'85%','60%','30%','Ctl'};
    elseif isGrad==2
        order=[4 3 2 1];
        colCel={{[0 0 0],[243 243 243]./255},... %Ctrl
            {[0 0 0], [208 238 252]./255},... %85
            {[0 0 0], [163 222 249]./255},... %60
            {[0 0 0], [90 187 234]./255},... %30
            };
        xlab={'30%','60%','85%','Ctl'};
    else
        order=1:length(responses);
        colCel={'k'};
        xlab=responses;
    end
    
    exOut=[]; isPaired=0; anovaFirst=0; mc=0; useWilcox=0; negError=1; plotMedians=0;
    [h,pAll,mseO]=barMeanSig2(d(order),xlab(order),{[1 2]},colCel,...
        exOut,isPaired,anovaFirst,mc,useWilcox,negError,plotMedians);
    set(h,'LineWidth',0.5,'UserData','Ex')
    hold on
    tot_mean=ci95mean([d{:}],2,1);
    plot([1 1].*length(xlab)+1,[tot_mean(2),tot_mean(3)],'g')
    hold on
    plot(length(xlab)+1,tot_mean(1),'.g')
    x=[];y=[];
    for i=1:length(order)
        y=[y d{order(i)}];
        xtemp=i+(rand(1,length(d{order(i)}))-0.5)/3;
        x=[x xtemp];
    end

    scatter(x,y,4,'sk')

    set(gca,'xtick',1:(length(xlab)+1),'xticklabel',[xlab(order) 'All'])
    plot([0 length(xlab)+1],[ci_null(1) ci_null(1)],'--k','LineWidth',0.5,'UserData','Ex')

    ylabel('% Accuracy')
    ylim([0 100])
    for i=1:(length(d))
        obs(i)=sum(d{i} > dn{i});
        exp(i)= round(length(d{i})*(mean(dn{i})/100));
        tot(i)=length(d{i});
    end
    title(sprintf('%2.1f%% Correct (%d of %d mice)',...
        (sum(actual==predicted)/length(actual))*100,sum(actual==predicted),...
        length(actual)))
    bi_Plot_Corrections
    set(gcf,'pos',[688   590   346   175])
    set(gca,'ytick',0:25:100)
    xlabel('TH Remaining')
    % Chi2
    ct=unique(actual);
    ci_null=ci95mean(pc_null,2,1);
    clear d dn
    alpha=0.05;
    for i=1:length(ct)
        ind=actual==ct(i);
        d{i}=pc(ind);
        dn{i}=pc_null(ind);
    end
    clear obs exp tot
    for i=1:(length(d))
        obs(i)=sum(d{i} > dn{i});
        exp(i)= round(length(d{i})*(mean(dn{i})/100));
        tot(i)=length(d{i});
    end
    
    [pchi(ii),c]=chi2_test_oe(obs,exp,[])

    % Single box plot with all correct mice  
    figure    
    hold on
    exOut=[]; isPaired=1; anovaFirst=0; mc=0; useWilcox=0; negError=1; plotMedians=0;
    [h,pAll,mseO]=barMeanSig2({pc,pc_null},{'Actual','Chance'},{[1 2]},{'k','b'},...
        exOut,isPaired,anovaFirst,mc,useWilcox,negError,plotMedians);
    chance(ii,:)=stdErrMean(pc_null,2,1);
    perf(ii,:)=stdErrMean(pc,2,1);
   
  

    colstr={'or','om','oc','ob','ok'};
    if ii==4
        ut=unique(origPerMouse);
        cols=othercolor('StepSeq_25',length(ut));
        clear h
        for i=1:length(ut)            
            ind=origPerMouse==ut(i);
            h(i)=scatter(1+rand(sum(ind),1)/2,pc(ind),'s','filled','MarkerFaceColor',...
                cols(i,:),'MarkerEdgeColor',cols(i,:));
        end
        legend(h,origRespLabels)
    else
        ut=unique(actual);
        
        for i=1:length(ut)
            ind=actual==ut(i);
            h(i)=scatter(1+rand(sum(ind),1)/2,pc(ind),'s','filled','MarkerFaceColor',...
                [1 1 1],'MarkerEdgeColor',[0 0 0]);
        end
    end

    set(gca,'ytick',0:25:100)
    xlim([0.5 2.5])
    set(gcf,'pos',[   680   482   297   496])
    ylabel('Accuracy (%)')
    bi_Plot_Corrections
end

%% Load Figure 7H data:
fn2='4_behavs_all_states+raw_no_ints_3';
[~,~,d]=xlsread([pn sprintf('%s.xlsx',fn2)]);
% load([pn sprintf('%s.mat',fn2)])

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

display('Finished loading.')

%% Figure 7H
r=5; c=2;
ind= asyn | allgrad | ctl | ba;
raw={fr syn cv sb};
xlab={'FR','SYN','IR','B'};
th_i=th(ind);
zb_i=zbehave1(ind);
ss=subplot_indexer(r,c);
bins=[0 15 48 73 100];
lab={'A','30','60','85','C'};
x=0:100;
useWilcox=[0 1 1 1];
saveFits=[];
useZ=1;
modelType={'poly2','poly3','poly3','poly1'};
usePC=1:length(xlab); %if applicable, otherwise == columns
Robust={'off','off','off','off'};
figure
for i=1:length(raw)
    subplot(r,c,i)
    if length(raw)==4 || i~=1
        ft = fittype( 'smoothingspline' );
        opts = fitoptions( 'Method', 'SmoothingSpline' );
        opts.SmoothingParam = 0.001; % smooth
    else
        ft = fittype( modelType{i} );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust=Robust{i};
    end
    if useWilcox(i)==1
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

% Summary
dat=saveFits;
endStage=dat(1:10,usePC);
earlyStage=dat(end-5:end,usePC);
cc2=zeros(size(dat));
pad=ones(30,length(usePC));
padDat=[bsxfun(@times,pad,dat(1,:));...
    dat(:,usePC);...
    bsxfun(@times,pad,dat(end,:))];
tt=size(pad,1)+round(size(endStage,1)/2);
for i=1:length(usePC)
    temp=xcorr(padDat(:,i),endStage(:,i),'none');
    temp2=temp(size(padDat,1):end);
    cc2(:,i)=temp2(31:131);
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

figure
subplot(2,2,[1,2])
ft = fittype('poly2');
opts = fitoptions( 'Method', 'LinearLeastSquares' );
[xData, yData] = prepareCurveData( th(ind), behav(ind,1) );
[fitresult, gof1] = fit( xData, yData, ft, opts );


 %Normalize so that -1 = pathological state:
dat=fitresult(x);
dat=bsxfun(@minus,dat,min(dat,[],1));
dat=bsxfun(@rdivide,dat,max((dat),[],1));
dat=dat * 1.25;
dat=dat-1;

imagesc(x,sort(unique(ccE(:,4))),dat');
hold on
mT=ccE(:,4);
mT=bsxfun(@minus,mT,min(mT,[],1));
mT=bsxfun(@rdivide,mT,max(mT,[],1));
mT=2*(mT-0.5);

h=plot(x,mT,'k');
a=othercolor('BrBG9',100);
a=(a(10:45,:));
colormap(a)
set(gca,'ydir','normal')
hold on
legend(h,{'Similarity to Naive','Behavior Score'},'location','southwest')
set(gca,'xdir','reverse','xtick',0:25:100)
xlabel('TH Remaining')
ylim([-1 1])
ylabel('z')

subplot(2,2,[3,4])
hold on
title('Behavior PC1 trajectory colorscale')
imagesc(x,sort(unique(mT)),dat');
set(gca,'xdir','reverse','xtick',0:25:100)
colorbar

set(gcf,'pos',[  474        1455         476         353])
bi_Plot_Corrections