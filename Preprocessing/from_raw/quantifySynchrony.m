 function [ data ] = quantifySynchrony( data )
% Tim Whalen March 2018, last edited October 2018
% Calculates cross-correlogram and confidence intervals for synchrony
% normalized to chance levels of synchrony
% Xcorr is computed over mall windows, normalized, hen averaged together to
% allow comparisons across pairs and average out nonstationarities in
% firing rate
% Baseline is calculated for lags > cut. One unit's spike train is zeroed
% to the left and right such that the same number of real time pointsare 
% used at each lag, allowing for a stable baseline (i.e., not artificially
% decreasing due to normal zero-padding).

% data must include: ts, T, spkchans_space, movet_raw, movey_raw, files (for title if toPlot==1)
%                    If using minRate, will try to use firing rates during rest calculated by
%                    restMoveVals - if this has not been run, rates is also required.

% To specify parameters, add the sub-struct "sync" to your data struct 
% containing any of the following you wish to define (all undefined will 
% use defaults):
%   bin = length of bins in sec (default = .01). (note blen-slen)/2 (i.e.
%         maxlag) must be 2)
%   slen = time length of train not zeroed in sec. Maxlag is (blen-slen)/2, 
%          so blen > slen required (default = 4, default with blen leads to 
%          default maxlag = 4)
%   blen = time length of unzeroed train in sec, i.e. maximum time length
%          of information considered at once. (default = 12)
%   step = time step for moving window in sec (default = 4, step > slen is 
%          not recommended and will throw away data in zeroed spike train 
%   cut = time lag from zero to consider baseline activity. Note that
%         baseline will be averaged from cut to maxlag, so sufficiently
%         long maxlag will mask oscillations (default = 0.5)
%   sides = number of bins left and right of zero to report (default = 2)
%   minRate = minimum rate for unit to be used in Hz (default = 5)
%   toPlot: 1 to plot cross-correlograms, 0 otherwise. (never plots
%            autocorrs)
% Note: because len sets consistent length of comparison, a longer input
% train only affects how accurate the baseline distribution (and confidence
% interval) is, but should not affect power of test
%
% Output: Fields will be added to the sync sub-struct of data; see line 65.

% Get parameters or set to defaults
names = {'bin', 'slen', 'blen', 'step', 'cut', 'sides', 'minRate', 'toPlot'};
defaults = {.01, 4, 12, 4, 0.5, 2, 5, 0};
[data, bin, slen, blen, step, cut, sides, minRate, toPlot] = extractDataInputs(data,'sync',names,defaults);
data.sync.rest = struct();

% Get necessary data - will throw error if missing
ts_all = data.ts;
Tall = data.T;
chans_all = data.spkchans_space; % channel numbers labeled ventral to dorsal
movet_all = data.movet_raw;
movey_all = data.movey_raw;
if toPlot
    files = data.('files');
end
if minRate ~=Inf
    try
        rates_all = data.rest.restrate;
        disp('Using rates during rest')
    catch
        rates_all = data.rates;
        disp('Unable to find rates during rest - using rates for whole train instead.')
    end
end

% Initialize outputs
xcorr_avg = cell(size(ts_all)); % averaged and normalized xcorr
xcorr_avg_rest = cell(size(ts_all)); % as above, but only for windows with no movement
around0_norm = cell(size(ts_all)); % normalized xcorr from -side:sides bins
around0_norm_rest = cell(size(ts_all));
pval0_mat_rest = cell(size(ts_all)); % bootstrapped p-values at lag zero for each pair
pval0_list_rest = cell(size(ts_all)); % as above, but as a list instead of matrix
sync_list = cell(size(ts_all)); % list of synchrony values for all pairs inach file
sync_list_rest = cell(size(ts_all)); % indicates which 2 units were used for each value in sync_list
sync_list_units = cell(size(ts_all));
sync_list_units_rest = cell(size(ts_all));
sync_indices = cell(size(ts_all)); % each unit's average normalized lag-zero xcorr for all pairs it is in
sync_indices_rest = cell(size(ts_all));

binm = 1/bin; % #bins/sec
maxlag_sec = (blen-slen)/2; % sec
maxlag = round(maxlag_sec*binm); % maxlag as @bins - MUST BE INT even before rounding; rounding just converts from double
stepbins = round(maxlag*2)+1; % #xcorr bins in one step
binsperstep = round(step*binm); % #bins moved in one step
bbins = round(binm*blen);   % #signal bins for unzeroed train in one step
if mod(bbins,2)==1 % if odd
    bbins = bbins-1;
end

% for zeroing
nkeepbins = round(slen*binm); % #signal bins not to zero in second train
if mod(nkeepbins,2)==1
    nkeepbins = nkeepbins-1;
end
binstart = bbins/2 - nkeepbins/2; % zero 2nd train before here
binend = bbins/2 + nkeepbins/2; % zero 2nd train after here

ncutbins = round(cut*binm); % #bins on either side of zero that aren't baseline
cutleft = (stepbins-1)/2 - ncutbins;
cutright = (stepbins-1)/2 + ncutbins;

if (blen-slen)/step==2 % this restriction makes things much easier
    resttime = zeros(size(ts_all));
else
    disp('warning: window size (blen-slen) is not step*2 - calculating time spent in rest skipped')
end

for f = 1:length(ts_all)
    movey = movey_all{f};
    if isempty(movey)
        continue
    end
    movet = movet_all{f};
    T = Tall(f);
    totalbins = floor(T*binm);
    nsteps = floor(T/step - blen/step) + 1;

    ts = ts_all{f};
    chans = chans_all{f};
    rates = rates_all{f};
    nu = length(ts); % # units
    if (nu==1) && (isempty(ts{1}))
        continue
    end
    xcorr_avg{f} = cell(nu);
    
    xcorr_avg_rest{f} = cell(nu);
    around0_norm_rest{f} = zeros(nu,nu,sides*2+1);
    around0_norm{f} = zeros(nu,nu,sides*2+1);
    pval0_mat_rest{f} = zeros(nu,nu);
    
    for i = 1:nu
        for j = 1:nu % do all and average later, since zeroing edges of one cell means xcorr is not symmetric
            if (chans(i) == chans(j) || (rates(i) < minRate || rates(j) < minRate)) % also skips i==j
                around0_norm{f}(i,j,:) = zeros(1,1,2*sides+1)+NaN;
                around0_norm_rest{f}(i,j,:) = zeros(1,1,2*sides+1)+NaN;
                continue
            end
            
            deltifull = zeros(1,totalbins);
            deltjfull = deltifull;
            deltifull( round(binm*ts{i})+1) = 1;
            deltjfull( round(binm*ts{j})+1) = 1;
            xcos = zeros(nsteps,stepbins);
            mving = zeros(nsteps,1);
            
            for s = 1:nsteps
                starti = (s-1)*binsperstep + 1; % start index;
                endi = (s-1)*binsperstep + bbins; % end index;
                delti = deltifull(starti:endi);
                deltj = deltjfull(starti:endi); % test: size should be bbins
                deltj_zeroed = deltj;
                deltj_zeroed([1:binstart-1 1+binend:bbins]) = 0; % test: sum ~=0 should be nkeepbins)
                xco = xcorr(delti,deltj_zeroed,maxlag);
                xb = xco([1:cutleft cutright:stepbins]);
                xb_avg = mean(xb);
                if xb_avg == 0 % if no spikes during interval (probably due to signal disrutpion)
                    xcos(s,:) = zeros(size(xco));
                else
                    xcos(s,:) = xco/xb_avg; % normalize baseline to 1; now x0 is multiple of baseline
                end
                % convert real times to mv indices
                mvs = find(movet>(starti-1)*bin,1); % start, note 1/bin is FS
                mve = find(movet>(endi-1)*bin,1)-1; % end
                mving(s) = sum(abs(movey(mvs:mve))>=2)>0; % 1 if any part of larger segment has mvmt 
            end
            xcorr_avg{f}{i,j} = mean(xcos,1); % avg all steps
            xcorr_avg_rest{f}{i,j} = mean(xcos(mving==0,:),1);
            around0_norm{f}(i,j,:) = xcorr_avg{f}{i,j}(maxlag+1-sides:maxlag+1+sides);
            around0_norm_rest{f}(i,j,:) = xcorr_avg_rest{f}{i,j}(maxlag+1-sides:maxlag+1+sides);
            xc_sort = sort(xcorr_avg_rest{f}{i,j}([1:cutleft cutright:stepbins]));
            nlower = find(xc_sort>around0_norm_rest{f}(i,j,sides+1),1)-2; % -1 for indexing, -2 because want upper bound
            if isempty(nlower)
                nlower = length(xc_sort)-1;
            end
            if nlower<0
                nlower = 0;
            end
            pval0_mat_rest{f}(i,j) = (length(xc_sort)-nlower)/length(xc_sort);
            
            % from mving, get length of time in rest used for sync
            if (blen-slen)/step==2 % condition makes this way easier to compute
                % Intuitively we want to compute rest time by walking through each
                % step adding window length for each 0 and step length for each
                % subsequent adjacent zero. This opaquely does that.
                mving2 = [1;mving;1];
                starts = sum(diff(mving2)==-1);
                resttime(f) = 2*step*starts + step*(sum(mving2==0)-starts);
            end
            
            if toPlot && i>j
                figure
                plot(-maxlag*bin:bin:maxlag*bin,(xcorr_avg_rest{f}{i,j}+xcorr_avg_rest{f}{j,i})/2,'k');
                title([data.files{f} ' Unit ' int2str(i) ' vs. ' int2str(j) ' (pval = ' num2str(pval0_mat_rest{f}(i,j)) ')'])%, T = ' num2str(resttime(f)) ' sec)'])
                xlabel('Time (sec)')
                ylabel('XCorr')
            end
        end
    end
    [ sync_list{f}, sync_list_units{f}, sync_indices{f} ] = triAverageAndIndex(around0_norm{f},pval0_mat_rest{f});
    [ sync_list_rest{f}, sync_list_units_rest{f}, sync_indices_rest{f}, pval0_list_rest{f} ] = triAverageAndIndex(around0_norm_rest{f},pval0_mat_rest{f});
    
end

data.sync.sync_mat = around0_norm;
data.sync.rest.sync_mat = around0_norm_rest;
data.sync.sync_list = sync_list;
data.sync.rest.sync_list = sync_list_rest;
data.sync.sync_list_units = sync_list_units;
data.sync.rest.sync_list_units = sync_list_units_rest;
data.sync.sync_indices = sync_indices;
data.sync.rest.sync_indices = sync_indices_rest;
data.sync.rest.pval0_mat = pval0_mat_rest;
data.sync.rest.pval0_list = pval0_list_rest;
if (blen-slen) == 2*step
    data.sync.rest.resttime = resttime;
end
end

% Subroutine to average across the diagonal of a matrix, and get some
% useful values from that. Separate snce we do it for full and only rest.
function [ triavg, whichunits, indices, pvals ] = triAverageAndIndex(mat,pval_mat)
    if length(size(mat))==2 % switch to 3D
        mat2 = zeros(size(mat,1),size(mat,2),1);
        mat2(:,:,1) = mat;
        mat = mat2;
    end
    n = size(mat,1);
    m = size(mat,3);
    mid = (m+1)/2; % element at lag 0
    arenan = isnan(triu(mat(:,:,1)));
    numnan = sum(arenan(:))-n; % number of NaN's not on diagonal
    whichunits = zeros((n*(n-1)/2)-numnan,2);
    triavg = zeros((n*(n-1)/2)-numnan,size(mat,3));
    pvals = zeros((n*(n-1)/2)-numnan,1);
    count = 1;
    all_indices = zeros(n,n-1);
    for i = 1:n
        for j = 1+i:n
            if isnan(mat(i,j,1))
                continue
            end
            triavg(count,:) = (mat(i,j,:)+mat(j,i,:))/2;
            pvals(count) = (pval_mat(i,j)+pval_mat(j,i))/2;
            whichunits(count,:) = [i j];
            all_indices(i,j-1) = triavg(count,mid);
            all_indices(j,i) = triavg(count,mid);
            count = count+1;
        end
    end
    indices = mean(all_indices,2);
end