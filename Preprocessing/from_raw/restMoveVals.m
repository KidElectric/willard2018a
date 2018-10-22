function [data] = restMoveVals( data )
% Tim C Whalen January 2018, last edited October 2018
% Calcute spiking (rate, CV, burst) and LFP (beta) vaues during rest and
% movement bouts (no burst values during movement)
%
% data must include: ts, movet_raw, movey_raw, lfp, lfpchans,
%                    outputs from surpriseBurst in burst sub-struct
%                    if do_lfpbeta==1, also include: lfp, lfpchans
%
% To specify parameters, add the sub-struct "rest" to your data struct
% containing any of the following you wish to define (all undefined will
% use defaults):
%   edge: time in sec to discard on either side of movement bout start and
%         stop (defualt = 0.5)
%   do_lfpbeta: 1 to calculate fractional LFP beta power, 0 to skip
%               (default = 1)
% 
% Output: Fields will be added to the rest sub-struct of data; see line 50.

%% Initializations

% Get parameters and defaults
names = {'edge','do_lfpbeta'};
defaults = {0.5, 1};
[data, edge, do_lfpbeta] = extractDataInputs(data,'rest',names,defaults);

% Get critical data; will throw error if not in struct
ts_all = data.ts;
starts_all = data.burst.begin;
blengths_all = data.burst.num_spikes;
movet_all = data.movet_raw;
movey_all = data.movey_raw;
if do_lfpbeta
    try
        lfp_all = data.lfp;
        fpchans_all = data.lfpchans;
    catch
        disp('LFP data or channels not found - skipping. To always skip, set do_lfpbeta to 0')
        do_lfpbeta = 0;
    end
end

if ~isfield(data,'rest')
    data.rest = struct();
end

nfiles = length(ts_all);

% Initialize outputs
data.rest.bps = cell(nfiles,1);
data.rest.nspikes = cell(nfiles,1); % number of spikes during rest
data.rest.nmvspikes = cell(nfiles,1); % number of spikes during movement
data.rest.nbspikes = cell(nfiles,1); % number of spikes in bursts
data.rest.nbursts = cell(nfiles,1);
data.rest.mbdur = cell(nfiles,1);
data.rest.restISI = cell(nfiles,1);
data.rest.mvISI = cell(nfiles,1);
data.rest.resttime = zeros(nfiles,1);
data.rest.mvtime = zeros(nfiles,1);
if do_lfpbeta
    data.rest.betafrac = zeros(nfiles,1);
    data.rest.betafracrest = zeros(nfiles,1);
    data.rest.betafracmv = zeros(nfiles,1);
end

%% Main loop over files
for f = 1:nfiles
    movet = movet_all{f};
    movey = movey_all{f};
    mv = movet(movey>=2);
    
    nu = size(ts_all{f},1);
    
    data.rest.bps{f} = zeros(nu,1);
    data.rest.nspikes{f} = zeros(nu,1);
    data.rest.nmvspikes{f} = zeros(nu,1);
    data.rest.nbspikes{f} = zeros(nu,1);
    data.rest.nbursts{f} = zeros(nu,1);
    data.rest.mbdur{f} = zeros(nu,1);
    data.rest.restISI{f} = cell(nu,1);
    data.rest.mvISI{f} = cell(nu,1);
    
    %% Mark movement bouts
    bouts = zeros(0,2);
    if ~isempty(mv)
        bouts(1,:) = [mv(1)-edge mv(1)+edge];
        for i = 2:length(mv)
            if mv(i)-edge>bouts(end,2) % if new mv bout is starting
                bouts(size(bouts,1)+1,:)=[mv(i)-edge mv(i)+edge];
            else
                bouts(end,2)=mv(i)+edge;
            end
        end
    end
    
    if isempty(bouts)
        b = [0 movet(end)];
    else
        bouts2 = bouts';
        boutslin = bouts2(1:end); % alternating start and end
        if iscolumn(boutslin)
            boutslin = boutslin';
        end
        if bouts2(1,1) <= 0
            b = boutslin(2:end);
        else
            b = [0 boutslin(1:end)];
        end
        if b(end)>movet(end)
            b = b(1:end-1);
        else
            b = [b movet(end)];
        end
    end
    rests = reshape(b,2,length(b)/2)'; % periods of rest
    Trest = sum(diff(rests,1,2)); % length of time in rest
    
    %% Loop over units in file
    for u = 1:nu
        ts = ts_all{f}{u};
        starts = starts_all{f}{u};
        blengths = blengths_all{f}{u};
        
        % find bursts
        if size(ts,1) ~=1
            ts = ts';
        end
        bs = zeros(1,sum(blengths)); % times of spikes in bursts
        count = 1; % one more than the current number of spikes in bursts, for indexing
        for i = 1:length(starts)
            bs(count:count+blengths(i)-1) = ts(starts(i):(starts(i)+blengths(i)-1));
            count = count+blengths(i);
        end
        
        nbursts = 0;
        nspikes = 0;
        nbspikes = 0;
        bdurs = [];
        restISI = [];
        resttime = 0;
        for i = 1:size(rests,1)
            wind = rests(i,:);
            nspikes = nspikes + sum(ts > wind(1) & ts < wind(2));
            nbspikes = nbspikes + sum(bs > wind(1) & bs < wind(2));
            somestarts = starts(ts(starts) > wind(1) & ts(starts) < wind(2));
            nbursts = nbursts + length(somestarts);
            ns = blengths(ts(starts) > wind(1) & ts(starts) < wind(2)); % nspikes in said bursts
            for b = 1:length(ns)
                bend = ts(somestarts(b)+ns(b)-1); % end of burst
                if bend > wind(2)+edge
                    break % burst leaks into movement
                else
                    bdurs = [bdurs bend-ts(somestarts(b))];
                end
            end
            % get CV, mean FR here
            restts = ts(ts > wind(1) & ts < wind(2));
            restISI = [restISI diff(restts)];
            resttime = resttime + abs(wind(2)-wind(1)); %quick fix for weird negative rests
            
        end
        mbdur = mean(bdurs); % mean burst duration
        
        mvISI = [];
        mvtime = 0;
        nmvspikes = 0;
        % each bout is at least 2*edge long, but real movement starts at
        % edge and ends at wind(2)-edge
        for i = 1:size(bouts,1)
            wind = bouts(i,:); % each bout by definition > edge (0.5 sec)
            wind = [wind(1)+edge wind(2)-edge];
            if diff(wind)<edge
                continue
            end
            mvtime = mvtime + diff(wind);
            mvts = ts(ts > wind(1) & ts < wind(2));
            mvISI = [mvISI diff(mvts)];
            nmvspikes = nmvspikes + size(mvts,2);
        end
        
        data.rest.bps{f}(u) = nbursts/Trest; % bursts/sec
        data.rest.nspikes{f}(u) = nspikes; % number of spikes in train during rest
        data.rest.nbspikes{f}(u) = nbspikes; % number of spikes in train that are in bursts
        data.rest.nbursts{f}(u) = nbursts;
        data.rest.mbdur{f}(u) = mean(bdurs); % mean burst duration
        data.rest.resttime(f) = resttime;
        data.rest.restISI{f}{u} = restISI;
        data.rest.nmvspikes{f}(u) = nmvspikes; % number of spikes in train during movement
        data.rest.mvtime(f) = mvtime;
        data.rest.mvISI{f}{u} = mvISI;
    end
    
    %% Beta power of LFP
    if do_lfpbeta
        disp(['Computing rest PSD of file ' int2str(f) ' of ' int2str(nfiles)])
        
        fpgood = lfp_all{f}(:,mod(fpchans_all{f},2)==1); % odd channels correctly referenced
        if isempty(fpgood)% no good FP channels
            data.rest.betafrac(f) = NaN;
            data.rest.betafracrest(f) = NaN;
            data.rest.betafracmv(f) = NaN;
            continue
        end
        fpm = mean(fpgood,2);
        [psd, freqs] = plomb(fpm, 1000, 100);
        bfrac = sum(psd(freqs>13 & freqs<30))/sum(psd(freqs>1 & freqs<100));
        
        % for rest and mvmt
        fpmrest = fpm;
        for i = 1:size(bouts,1)
            wind = bouts(i,:);
            if wind(1)< 0
                wind(1) = 0;
            end
            bouti = [floor(wind(1)*1000+1) ceil(wind(2)*1000+1)]; % indices to NaN
            if bouti(1) < size(fpmrest,1)
                if bouti(2) < size(fpmrest,1)
                    fpmrest(bouti(1):bouti(2)) = NaN;
                else
                    fpmrest(bouti(1):end) = NaN;
                    break % reached end of recording
                end
            else
                break
            end
        end
        try
            [psdrest, freqsrest] = plomb(fpmrest, 1000, 100);
            bfracrest = sum(psdrest(freqsrest>13 & freqsrest < 30))/sum(psdrest(freqsrest>1 & freqsrest < 100));
        catch
            bfracrest = NaN; % not enough data;
        end
        fpmmv = fpm;
        fpmmv(isnan(fpmrest)==0)=NaN; % everything not NaN in rest is NaN in mv
        try
            [psdmv, freqsmv] = plomb(fpmmv, 1000, 100);
            bfracmv = sum(psdmv(freqsmv>13 & freqsmv < 30))/sum(psdmv(freqsmv>1 & freqsmv < 100));
        catch
            bfracmv = NaN; % not enough data
        end
        data.rest.betafrac(f) = bfrac;
        data.rest.betafracrest(f) = bfracrest;
        data.rest.betafracmv(f) = bfracmv;
    end
end
end