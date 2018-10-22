function  aw2018_figure_1h()
% Tim C Whalen, October 2018
% Plots three example normalized cross-correlograms - the second is the
% left panel in Fig 1H. See AWanalysis_fun for full details of CCG analysis

load('AWgradual_ex4sync','data'); % example processed data
data.sync = struct();
data.sync.toPlot = 1;
[data] = quantifySynchrony(data);
end