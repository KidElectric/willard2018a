function bi_Plot_Corrections(h,params)
%bi_Plot_Corrections(h) - apply basic plot corrections to improve
% quality of MATLAB plots
% I.e.- it finds all the objects (from all subplots,axes etc)
% and changes them to better defaults. BI 20015
%
% INPUT:
% h = figure handle (i.e. gcf) - Optional.
% params = work in progress structure of default axis features that can be fed
%          into to make specific changes
%       params.box = 'on' or 'off'
%       params.axisFontSize = set axis font size
%       params.labelFontSize =  axis label font size
%       params.TickLength = axis tick length
%       params.LineWidth = line widths
% Note: can also be called without any inputs, in which case it will
% simply grab the handle of the current matlab figure;
% See below for details.

if ~exist('h','var') || isempty(h)
    h=gcf;
end
hh=findobj(h,'Type','Line');
if exist('params','var') && isfield(params,'LineWidth')
    lw=params.LineWidth;
else
    lw=0.5;
end
if ~isempty(hh) 
    use=~strcmp(get(hh,'UserData'),'Ex');
    set(hh(use),'Markersize',10,'LineWidth',lw)
end
hM=findobj(h,'Marker','.');
if ~isempty(hM)
    use=~strcmp(get(hM,'UserData'),'Ex');
    if exist('params','var') && isfield(params,'MarkerSize')
        set(hM(use & isprop(hM,'MarkerSize')),'Markersize',params.MarkerSize)
    else
        set(hM(use & isprop(hM,'MarkerSize')),'Markersize',7)
    end
end
% set(hM,'MarkerFaceColor','none')

%Axis numbering
hA=findobj(h,'Type','axes');
if exist('params','var') && isfield(params,'box')
    set(hA,'box',params.box)
else
    set(hA,'box','off')
end
% if exist('params','var') && isfield(params,'axisLineWidth')
%     if ~isempty(params.axisLineWidth)
%         set(hA,'linewidth',params.axisLineWidth,'Color','none')
%     end
% else
%     set(hA,'linewidth',0.5,'Color','none')
% end
% hA=findobj(h,'Type','axes');
set(hA,'TickLength',[0.05, 0.01])
if exist('params','var') && isfield(params,'axisFontSize')
    set(hA,'FontSize',params.axisFontSize) %'FontWeight','bold'
else
    use=~strcmp(get(hA,'UserData'),'Ex');
    if any(use)
        if exist('params','var') && isfield(params,'FontSize')
            fs=params.FontSize;
        else
            fs=14;
        end
        set(hA(use),'FontSize',fs) %'FontWeight','bold'
    end
    
end

%All other font:
hT=findall(h,'Type','text');
if ~isempty(hT)
    use=~strcmp(get(hT,'UserData'),'Ex');
    if any(use)
        if exist('params','var') && isfield(params,'labelFontSize')
            set(hT(use),'FontSize',params.labelFontSize) %'FontWeight','bold'
        else
            set(hT(use),'FontSize',14) %'FontWeight','bold'
        end
    end
end
hL=findobj(h,'Type','axes','Tag','legend');
set(hL,'EdgeColor',[1 1 1])% 'box','off'
set(h,'Color',[1 1 1])
set(hA,'TickDir','out')
if exist('params','var') && isfield(params,'TickLength')
    set(hA,'TickLength',params.TickLength)
else 
    tl=get(hA,'TickLength');
    if iscell(tl)
        on=~sum(bsxfun(@eq,vertcat(tl{:}),[0 0]),2)==2;
    else
        on=~sum(bsxfun(@eq,tl,[0 0]),2)==2;
    end
    
    if any(on)
        set(hA,'TickLength',[0.02 0.02])
    end
end
if exist('params','var') && isfield(params,'Xoff')
    if params.Xoff==true
        set(hA,'xcolor',get(h,'color'))
    end
end

if exist('params','var') && isfield(params,'ylim')
    set(hA,'ylim',params.ylim)
end
% set(hA,'Xcolor',[1 1 1])