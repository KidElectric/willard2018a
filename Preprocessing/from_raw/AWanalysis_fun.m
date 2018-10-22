function data = AWanalysis_fun(mat,cond,append)
% Tim Whalen, Jan 2018. Last edited October 2018
% Master function for Willard et al - opens processed .mat and genrates
% xls files including firing rate, CV, surprise burst, synchrony, and LFP
% beta.
% Inputs: mat: the .mat data file to load
%         cond: name of condition for m, to use as xls sheet name
%         append: string (typically date) to append to xls filenames
% Outputs: data: input struct with burst, rest and sync sub-structs added
%          also, two xls files in your current directory, AWsync for
%          synchrony and AWvals for other values, each appended with
%          "_append.xls"

load(mat,'data');

%% Bursts, FR, CV, beta

disp('Detecting bursts...');
data.burst = struct('fac', 2,'sampling_rate', 1000, 'min_length_of_burst', 3, 'local_length', 0, 'surprise_cutoff', 2, 'max_ISI', Inf);
data = surpriseBurst(data); % Burst detection
data = restMoveVals(data); % Separate rest and movement

%% Synchrony

minSyncRest = 20; % seconds of rest (from sync windows) needed to use file for synchrony

data.sync = struct();
% data.sync.toPlot = 1;
disp('Calculating synchrony (10ms bin)...')
[data] = quantifySynchrony(data);
sync_indices10 = data.sync.rest.sync_indices;
sync_resttime = data.sync.rest.resttime; % rest is defined more stringently for sync than other analyses

disp('Writing synchrony to Excel...')
xlsfilename = ['AWsync_' append];
sync_list_long = data.sync.rest.sync_list(sync_resttime>=minSyncRest); % sync_list only with long enough periods of rest
lenlist = sum(cellfun(@(x) size(x,1),sync_list_long));
xl = cell(lenlist,9);
xl(1,:) = {'Filename','Unit A','Unit B' ...
    [int2str(data.sync.bin*1000*-2) ' ms'] ...
    [int2str(data.sync.bin*1000*-1) ' ms'] ...
    '0 ms' ...
    [int2str(data.sync.bin*1000) ' ms'] ...
    [int2str(data.sync.bin*1000*2) ' ms'] ...
    'p-value upper bound (0 lag)'};

count = 2;
for f = 1:data.nfiles
    if ~isempty(data.sync.rest.sync_list{f})
        for p = 1:size(data.sync.rest.sync_list_units{f},1)
            xl(count,:) = {data.files{f} data.sync.rest.sync_list_units{f}(p,1) data.sync.rest.sync_list_units{f}(p,2) ...
                data.sync.rest.sync_list{f}(p,1) ...
                data.sync.rest.sync_list{f}(p,2) ...
                data.sync.rest.sync_list{f}(p,3) ...
                data.sync.rest.sync_list{f}(p,4) ...
                data.sync.rest.sync_list{f}(p,5) ...
                data.sync.rest.pval0_list{f}(p)};
            count = count+1;
        end
    end
end
xlswrite(xlsfilename,xl,cond,'A1');

%% Write rest to excel

disp('Writing rest of data to Excel...')
% Make cell for excel export
nu = sum(cellfun(@(x)  size(x,1), data.ts)); % total @units
xl = cell(nu,30);
xl{1,1} = 'Filename';
xl{1,2} = 'Unit';
xl{1,3} = 'Channel';
xl{1,4} = 'Rate (Hz)';
xl{1,5} = 'Total Spikes';
xl{1,6} = 'Total Spikes in Bursts';
xl{1,7} = 'Mean Spikes/Burst';
xl{1,8} = 'Fraction Time in Bursts';
xl{1,9} = 'Mean Burst Duration (sec)';
xl{1,10} = 'Bursts/Sec';
xl{1,11} = 'Mean Intra-Burst ISI';
xl{1,12} = 'Mean Surprise';
xl{1,13} = 'Total Spikes (REST ONLY)';
xl{1,14} = 'Total Spikes in Bursts (REST ONLY)';
xl{1,15} = 'Mean Burst Duration (sec) (REST ONLY)';
xl{1,16} = 'Bursts/Sec (REST ONLY)';
xl{1,17} = 'Firing Rate';
xl{1,18} = 'Firing Rate (REST ONLY)';
xl{1,19} = 'Firing Rate (MOVE ONLY)';
xl{1,20} = 'CV ISI';
xl{1,21} = 'CV ISI (REST ONLY)';
xl{1,22} = 'CV ISI (MOVE ONLY)';
xl{1,23} = 'Fractional Beta Power';
xl{1,24} = 'Fractional Beta Power (REST ONLY)';
xl{1,25} = 'Fractional Beta Power (MOVE ONLY)';
xl{1,26} = 'Synchrony Index 10 ms (REST ONLY)';
xl{1,27} = 'Time at Rest Used for Index';

count = 2;
for f = 1:data.nfiles
    for u = 1:length(data.ts{f})
        xl{count,1} = data.files{f};
        xl{count,2} = u;
        xl{count,3} = data.spkchans{f}(u);
        xl{count,4} = data.rates{f}(u);
        xl{count,5} = length(data.ts{f}{u});
        xl{count,6} = data.burst.total_spikes_in_bursts{f}(u);
        xl{count,7} = data.burst.mean_spikes_per_burst{f}(u);
        xl{count,8} = data.burst.proportion_time_in_bursts{f}(u);
        xl{count,9} = mean(data.burst.burst_duration{f}{u});
        xl{count,10} = data.burst.num_bursts{f}(u)/data.T(f);
        xl{count,11} = mean(1./data.burst.rate{f}{u});
        xl{count,12} = mean(data.burst.surprise{f}{u});
        xl{count,13} = data.rest.nspikes{f}(u);
        xl{count,14} = data.rest.nbspikes{f}(u);
        xl{count,15} = data.rest.mbdur{f}(u);
        xl{count,16} = data.rest.bps{f}(u);
        xl{count,17} = length(data.ts{f}{u})/data.T(f);
        xl{count,18} = data.rest.nspikes{f}(u)/data.rest.resttime(f);
        xl{count,19} = data.rest.nmvspikes{f}(u)/data.rest.mvtime(f);
        
        ISI = diff(data.ts{f}{u});
        xl{count,20} = std(ISI)/mean(ISI);
        rISI = data.rest.restISI{f}{u};
        xl{count,21} = std(rISI)/mean(rISI);
        mISI = data.rest.mvISI{f}{u};
        xl{count,22} = std(mISI)/mean(mISI);
        xl{count,23} = data.rest.betafrac(f);
        xl{count,24} = data.rest.betafracrest(f);
        xl{count,25} = data.rest.betafracmv(f);
        if sync_indices10{f}(u) ~=0 && sync_resttime(f)>=minSyncRest
            xl{count,26} = sync_indices10{f}(u);
        end
        xl{count,27} = sync_resttime(f);
        count=count+1;
    end
end
xlsfilename = ['AWvals' append];
xlswrite(xlsfilename,xl,cond);
disp('Saving .mat...')
save(strcat(xlsfilename,'_',cond),'xl','cond','data','-v7.3');

end