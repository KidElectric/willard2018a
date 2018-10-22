function aw2018_figure_1(pn)
% function aw2018_figure_1(pn)
%   pn = data path. if empty, assumes data are in a folder 'data' such that:
%   current directory is 'code,' and 'data' and 'code' are in 'fig' folder
% Load physiology data
% Read AW snr data .xlsx
if ~exist('pn','var') || isempty(pn)
    pn=['.' filesep '..' filesep 'data' filesep];
end
%
%% Figure 1I see:

%% Figure 1J & 1K -- compare classifications
fns={'mnrfit_500_perm_jackknife_wta_grad5_ba_psp.mat',...
    'mnrfit_500_perm_jackknife_endstage_mice_medneur_wta_depgt20_vs_ctrl.mat'};
for ii=1:length(fns)
    load([pn fns{ii}],'pc','pc_null','conf_mat_norm_keep','beta','beta_p',...
        'beta_p_null','beta_null','conf_mat_norm_total',...
        'responses','useResponse','origResp','origRespLabels',...
        'combine','origPerMouse','predictors','sheets','cols',...
        'fitType','conf_null_mat_norm_total','predicted','actual',...
        'ints','useInts','predicted_null');
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

%% Figure 1L - Plot betas
preds={'Int','B','FR','IR','SYN'};
reOrder=[3 4 2 5]-1;
if useInts==1
    ylab=preds(2:end);
    for i=1:length(ints)
        newLab{i}=sprintf('%s x %s',ylab{ints(i,1)},ylab{ints(i,2)});
    end
    allLab=[ {'int'} ylab(reOrder) newLab];
else
    allLab=predictors;
end
reOrder=[1 3 4 2 5];
beta(1:5,:,:)=beta(reOrder,:,:);
beta_p(1:5,:,:)=beta_p(reOrder,:,:);
for ii=1:size(beta,2)
    pk=median(beta_p(:,ii,:),3); %Very skewed- use median?
    b=ci95mean(squeeze(beta(:,ii,:)),2,0);
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
    plot([0 0],[1, size(beta,1)+0.5],'--k')
    xlabel('Coefficients')
    set(gcf,'pos',[    680   616   498   362])
    bi_Plot_Corrections
end