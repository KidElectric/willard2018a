function data = surpriseBurst(data)
% Tim Whalen Aug 2017, modification of code by Thomas Wichmann July 2005
% Carries out the burst detection algorithm by Legendy and Salcman (J. Neurophysiology 53:926) on
% the input ISI data stream. As the initial criterion for burst detection, spikes
% have to occur at a frequency which is higher than the baseline frequency by a factor fac.  The
% user can modify the minimal length of a prospective burst (min_length_of_burst), and the surprise
% cutoff value.  In the original paper, this value was set to 10 to detect significant spikes.
% The burst discharge is compared to a segment of data immediately preceding the burst.  A local_length
% value of 0 indicates that the entire data stream should be used to yield the firing rate (can be used
% in case of homogeneous firing of units).  Other values indicate a local time interval that should be
% used (for instance, a value of 1 would mean that each spike/burst is compared to 1 second of data
% prior to the index spike).
%
% The function makes use of a subroutine called surprise_new3.  The function returns 100 for
% a very high surprise value and 0 for a very low one.  This means that the output of this rou-
% tine is not accurate for such very high or very low surprise values (although bursts are correctly
% detected).
%
% The function produces or adds to the "burst" struct which contains fields describing the bursts, including the
% onset and lengths of bursts (described as indices of the input ISI stream) the highest rate within
% the burst, the average discharge rate in the burst, the pre-burst baseline rate, as well as the
% surprise values for individual bursts.  In field 1 of the structure, summary parameters are defined,
% including, num_bursts (the total number of bursts detected), mean_spikes_per_burst, total_spikes_in_bursts,
% mean_intra_burst_frequency, proportion_time_in_bursts, and proportion_spikes_in_bursts.

% TCW: fixed many errors which poorly labeled starts and ends of bursts.
% Comments prefaced by TW indicate changes.
% Also added real timescale by specifying max_ISI, the largest ISI in msec
% that can be considered part of a burst. This keeps non-burst series of
% ISI's with high surprise from being considered bursts. However, with
% above fixes, this parameter is less useful and may not be necessary.
% If max_ISI is small enough, this generally overrides the fac parameter, 
% making fac mostly unnecessary
% Default = Inf, i.e. no maximum.

% data must include: ts
% To specify parameters, add the sub-struct "burst" to your data struct 
% containing any of the following you wish to define (all undefined will 
% use defaults):
%   fac = firing rate threshold (as multiple of baseline) to initially
%       consider burst (defualt = 1.5)
%   sampling_rate: rate (in Hz) to downsampole spike times (default = 10000)
%       min_length_of_burst: minimum number of ISI's a burst can have
%       (default = 2)
%   local_length: time (in sec) in past over which to define current
%       baseline firing rate. If 0, calculates over entire train (default = 0)
%   surprise_cutoff: minimum Poisson surprise for a burst to be kept
%       (default = 10)
%   max_ISI: maximum ISI allowed within burst, overriding fac (added by
%       TCW, default = Inf)

%% Initializations

% First get critical data; will throw error if not in struct
ts = data.ts;

% Get parameters and defaults
names = {'fac', 'sampling_rate', 'min_length_of_burst', 'local_length', 'surprise_cutoff', 'max_ISI'};
defaults = {1.5, 1000, 2, 0, 10, Inf}; % Wichmann
% defaults = {2, 1000, 4, 0, 2, Inf}; % Mastro et al. 2017
[data, fac, sampling_rate, min_length_of_burst, local_length, surprise_cutoff, max_ISI] = extractDataInputs(data,'burst',names,defaults);

nfiles = length(ts);

if ~isfield(data,'burst')
    data.burst = struct();
end

data.burst.begin = cell(nfiles,1);
data.burst.num_spikes = cell(nfiles,1);
data.burst.surprise = cell(nfiles,1);
data.burst.rate = cell(nfiles,1);
data.burst.max_rate = cell(nfiles,1);
data.burst.baseline_rate = cell(nfiles,1);
data.burst.num_bursts = cell(nfiles,1);
data.burst.mean_spikes_per_burst = cell(nfiles,1);
data.burst.median_spikes_per_burst = cell(nfiles,1);
data.burst.total_spikes_in_bursts = cell(nfiles,1);
data.burst.mean_intra_burst_frequency = cell(nfiles,1);
data.burst.median_intra_burst_frequency = cell(nfiles,1);
data.burst.proportion_time_in_bursts = cell(nfiles,1);
data.burst.proportion_spikes_in_bursts = cell(nfiles,1);
data.burst.burst_duration = cell(nfiles,1);

%% Loop over files
for f = 1:nfiles;
    nu = length(ts{f});
    
    data.burst.begin{f} = cell(nu,1);
    data.burst.num_spikes{f} = cell(nu,1);
    data.burst.surprise{f} = cell(nu,1);
    data.burst.rate{f} = cell(nu,1);
    data.burst.max_rate{f} = cell(nu,1);
    data.burst.baseline_rate{f} = cell(nu,1);
    data.burst.burst_duration{f} = cell(nu,1);
    
    data.burst.num_bursts{f} = zeros(nu,1);
    data.burst.mean_spikes_per_burst{f} = zeros(nu,1);
    data.burst.median_spikes_per_burst{f} = zeros(nu,1);
    data.burst.total_spikes_in_bursts{f} = zeros(nu,1);
    data.burst.mean_intra_burst_frequency{f} = zeros(nu,1);
    data.burst.median_intra_burst_frequency{f} = zeros(nu,1);
    data.burst.proportion_time_in_bursts{f} = zeros(nu,1);
    data.burst.proportion_spikes_in_bursts{f} = zeros(nu,1);
    for u = 1:nu
        ISI = diff(ts{f}{u});
        
        burst_num = 0;
        CA = cumsum(ISI);
        
        if local_length == 0,
            mean_FR = length(ISI)/(sum(ISI)/sampling_rate);
            fr_thr = min(sampling_rate/(fac*mean_FR), max_ISI/1000);   % calculation of frequency threshold
            beg_idx = 0;
        else
            beg_idx = find(CA < local_length*sampling_rate,1,'LAST');     % finds last index within the 'local length' - incremented by one, this will result in the first index that can be evaluate.
        end
        n = beg_idx;
        
        % ***** Main loop ****
        
        while n < length(ISI) - min_length_of_burst
            
            n = n+1;      
            if local_length > 0,
                I = ISI(find(CA > CA(n)-local_length*sampling_rate,1,'FIRST'):n-1);     % find the ISI data segment I that is fully contained within the local_length
                mean_FR = length(I)/(sum(I)/sampling_rate);
                fr_thr = min(sampling_rate/(fac*mean_FR),max_ISI/1000);                                   % calculation of frequency threshold
            end
            
            
            % ****** 1. selection step - find runs of short ISIs *******************
            if (ISI(n) < fr_thr) % finding areas of the spike train which fulfill the length_of_burst criterion
                
                % TW 08/2017 - simplified, inc_flag was unnecessary and
                % sometimes caused bursts = min_length_of_burst to go
                % undetected
                q = 0;           % running parameter that points to the number of spikes to be added
                while (n+q <= length(ISI)) & (ISI(n+q) < fr_thr)
                    q = q+1;
                end
                q = q-1;
                % TW - this is actually quite inefficient - could just stop
                % once q+1 < min_length - but didn't bother to change it...
                
                % at this point, the provisional burst starts at n and ends at n+q;
                % it has q+1 spikes in it
                
                
                % ******* 2. selection step - adjust length of burst to maximize surprise value ***************
                if q+1 >= min_length_of_burst,
                    
                    m = min_length_of_burst; % TW: try one spike more than necessary
                    % running parameter of the number of spikes to be added to n
                    while ((n+m <= length(ISI)) & ...
                            (ISI(n+m) < fr_thr) & ...
                            (surprise_new3(mean_FR,ISI(n:n+m),sampling_rate) >= surprise_new3(mean_FR,ISI(n:n+m-1),sampling_rate))),   % 2. burst criterion - surprise values are increasing when spikes are added
                        m = m+1; % TW: try one more spike
                    end
                    m = m-1; % TW: subtract last try, even w/o inc_flag
                    % here, and inc_flag was irrelevant.
                    
                    % at this point, the beginning of the burst is provisionally settled to be n, the end n+m
                    % the burst has m+1 spikes in it. % TW 08/2017: this
                    % condition causes many n:m's to change to n:m-1's
                    % (e.g., if m==1, there should be 2 spikes and 1 ISI,
                    % so n:n+m-1 = n:n is a 1-element vector of ISI's).
                    
                   
                    % ******* 3. selection step - test whether adding up to 10 more ISIs will enhance surprise value **********************
                    if n+m+10 <= length(ISI) % mmax is set to 10 unless one couldn't add 10 to n before reaching the end of FR
                        mmax = 10;
                    else
                        mmax = length(ISI)-(n+m);
                    end
                    
                    ind_long_ISI = find(ISI(n+m:n+m+mmax) > fr_thr,1,'FIRST'); % TW 08/2017: m+1 to m     % looking for 'slow spikes' within the next 10 spikes after the current burst end
                    if ~isempty(ind_long_ISI)                                           % pmax is set to be the index of the slow spike
                        pmax = ind_long_ISI-1;
                    else
                        pmax = mmax;
                    end
                    
                    S = surprise_new3(mean_FR,ISI(n:n+m-1),sampling_rate); % TW: m -> m-1
                    for p = 1:pmax                                  % forms array of surprise values for this burst, starting from the end of the burst to pmax values later
                        S(p+1) = surprise_new3(mean_FR,ISI(n:n+m-1+p),sampling_rate); % TW: m -> m-1
                    end

                    if n+m < length(ISI)
                        [max_S,ind_max_S] = max(S);
                        if (max_S > surprise_new3(mean_FR,ISI(n:n+m-1),sampling_rate)) % TW: m -> m-1    % check whether the maximal surprise value in the post-burst period is higher than the maximal burst surprise
                            m = m+ind_max_S-1; % TW add -1
                        end
                    else
                        m = length(ISI)-n; % TW: don't think this can ever be reached?
                    end
                    
                    if n+m > length(ISI)
                        m = length(ISI)-n; % TW: also seems impossible to reach, but maybe with earlier off-by-1 errors
                    end
                    % at this point, the end of the index of the end of the burst
                    % is settled to be n+m (TW: n+m-1, because n:n+m-1 is m
                    % ISI's is m+1 spikes as desired)
                    
                    
                    % ******** 4. selection step - test whether adjusting the front end of the burst enhances the surprise value ******************
                    % TW: this section only checks removing one spike at a
                    % time until surprise is lowered, but original
                    % algorithm must check all possible starts (more like
                    % step 3)
                    % (editorializing here - this original algorithm seems
                    % like it could undersample bursts by accidentally
                    % finding two bursts in step 3 and removing the
                    % original burst in step 4 if the 2nd is more
                    % surprising than the first, but that's how it is. Also
                    % better than the incorrect step 4 in the original code
                    % which keeps the two bursts as a single burst)
                    
                    S = surprise_new3(mean_FR,ISI(n:n+m-1),sampling_rate); % surprise of unshortened train
                    for p = 1: m+1-min_length_of_burst % try all possible bursts, with smallest being min_length_of_burst spikes all at end
                        S(p+1) = surprise_new3(mean_FR,ISI(n+p:n+m-1),sampling_rate); % TW: m -> m-1
                    end
                    [max_S,ind_max_S] = max(S);
                    n = n+ind_max_S-1; % move start up
                    m = m-ind_max_S+1; % adjust length accordingly
                    
%                     if n > 1
%                         o = 1;
%                         while ((m-o > min_length_of_burst) && (surprise_new3(mean_FR,ISI(n+o:n+m-1),sampling_rate) >= surprise_new3(mean_FR,ISI(n+o-1:n+m-1),sampling_rate))) % TW: m -> m-1
%                             o = o+1;
%                         end
%                         % TW: as above, inc_flag lines were added to
%                         % correct for indexing off-by-1's but do so
%                         % erroneously
%                         o = o - 1;          % reducing o by one to correct for the addition that resulted in the end of the while loop
%                         n = n+o;            % adjust the beginning of the burst
%                         m = m-o;            % adjust the length of the burst
%                     end
                    
                    % at this point, the beginning of the burst is settled to be n, and the length is m+1 ***
                    if (m+1 >= min_length_of_burst) && max_S > surprise_cutoff, % TW: m -> m-1, and switched to max_S from step 4, no need to recalculate
                        burst_num = burst_num + 1;
                        data.burst.begin{f}{u}(burst_num) = n;
                        data.burst.num_spikes{f}{u}(burst_num) = m+1;
                        data.burst.surprise{f}{u}(burst_num) = max_S; % TW: again no need to recalculate
                        data.burst.rate{f}{u}(burst_num) = length(ISI(n:n+m-1))/(sum(ISI(n:n+m-1))/sampling_rate); % TW: m -> m-1
                        data.burst.max_rate{f}{u}(burst_num) = sampling_rate/min(ISI(n:n+m-1)); % TW: m -> m-1
                        data.burst.baseline_rate{f}{u}(burst_num) = mean_FR;
                        data.burst.burst_duration{f}{u}(burst_num) = sum(ISI(n:n+m-1));
                    end
                    
                    n = n+m; % TW: because n immediately incremented at start of loop, erroneously forcing bursts to be non-adjacent                                                      % adjust ISI pointer to the ISI following the burst
                    
                end
            end
        end
        
        
        % ****** Store burst parameters in the output array
        
        data.burst.num_bursts{f}(u) = length(data.burst.begin{f}{u});
        C = [];
        for n = 1:data.burst.num_bursts{f}(u)
            C(n) = data.burst.num_spikes{f}{u}(n);
        end
        data.burst.mean_spikes_per_burst{f}(u) = mean(C);
        data.burst.median_spikes_per_burst{f}(u) = median(C);
        data.burst.total_spikes_in_bursts{f}(u) = sum(C);
        
        for n = 1:data.burst.num_bursts{f}(u)
            C(n) = data.burst.rate{f}{u}(n);
        end
        data.burst.mean_intra_burst_frequency{f}(u) = mean(C);
        data.burst.median_intra_burst_frequency{f}(u) = median(C);
        data.burst.proportion_time_in_bursts{f}(u) = (data.burst.total_spikes_in_bursts{f}(u)/data.burst.mean_intra_burst_frequency{f}(u))/(sum(ISI(beg_idx+1:length(ISI)))/sampling_rate);
        data.burst.proportion_spikes_in_bursts{f}(u) = data.burst.total_spikes_in_bursts{f}(u)/length(ISI(beg_idx+1:length(ISI)));
    end
end
% ******* Exit ***********

%**********************************************************************
function [burst,deceleration] = surprise_new3(r,data,sampling_rate)
% calculates surprise index.  Parameters are ...
% r = comparison firing rate (spikes per second)
% data = ISI data to be included in the burst
% sampling_rate = sampling rate of ISI measurements

T = sum(data)/sampling_rate;
num_spikes = length(data);

p = poisscdf(num_spikes,r*T);

switch p
    case 0                                      % for very small p, a value of 10exp-100 is assumed
        burst = 0;
        if nargout > 1,deceleration = 100;end
    case 1                                      % for very high p, a value of 1-10exp-100 is assumed
        burst = 100;
        if nargout > 1,deceleration = 0;end
    otherwise
        burst = -log(1-p);
        if nargout > 1,deceleration = -log10(p);end
end


%************************************************************************