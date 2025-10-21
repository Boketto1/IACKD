function align_eeg_leap(eeg_seg_mat, lm_mat, out_dir, tol_ms, skip_pen_ms)
% Monotonic DP alignment: Sequential pairing of durations between EEG (14->66) and LEAP (start->hit)
%
% Usage:
% align_eeg_leap('...\segments_14_to_1000001.mat', '...\leap_trials_100Hz.mat', '...\out', 600, 2000)

if nargin < 4 || isempty(tol_ms),      tol_ms = 600;     end
if nargin < 5 || isempty(skip_pen_ms), skip_pen_ms = 2000; end
if ~exist(out_dir,'dir'), mkdir(out_dir); end

E = load(eeg_seg_mat);   % Need: segs_blc, times, trial_info, RESAMPLE_HZ, and keep_mask_* or keep_mask
S = load(lm_mat);        % Need: LM structure (incluting duration_sec, t, x_mm/y_mm/ball_x_pix/ball_color)

% ---------Get EEG trial list---------
reqE = {'segs_blc','times','trial_info','RESAMPLE_HZ'};
for k=1:numel(reqE), assert(isfield(E,reqE{k}), 'EEG missing: %s', reqE{k}); end

% keep mask selection
if isfield(E,'keep_mask_segments')
    keep_seg = logical(E.keep_mask_segments(:));
elseif isfield(E,'keep_mask')
    keep_seg = logical(E.keep_mask(:));
else
    warning('keep_mask_segments" and "keep_mask" not found, default to retain all.');
    keep_seg = true(numel(E.segs_blc),1);
end

% The mapping from seg to info (in the order of trial_info-valid)
valid_info_idx = find([E.trial_info.valid]==1);
assert(~isempty(valid_info_idx), 'There is no entry for valid==1 in trial_info');
assert(numel(valid_info_idx) >= sum(keep_seg), 'The valid_info_idx is less than the number of segments. Please check the upstream export.');

seg_idx_list = find(keep_seg);      
info_idx_list = valid_info_idx(seg_idx_list);

% EEG duration (ms)
dur_eeg_ms = arrayfun(@(i) E.trial_info(i).dur_ms, info_idx_list).';
nan_mask = isnan(dur_eeg_ms);
if any(nan_mask)
    for ii = find(nan_mask).'
        si = seg_idx_list(ii);
        dur_eeg_ms(ii) = E.times{si}(end)*1000;
    end
end

% ---------Retrieve the LEAP trial list (only valid==true)---------
assert(isfield(S,'LM'),'LEAP: Lack of LM structure');
LMall = S.LM(:);
Lmask = true(size(LMall));
if isfield(LMall,'valid')
    Lmask = logical([LMall.valid]);
end
LM = LMall(Lmask);
assert(~isempty(LM),'LEAP: No valid trial');

lm_trial_id = arrayfun(@(u) u.trial, LM).';
dur_lm_ms   = arrayfun(@(u) u.duration_sec, LM).' * 1000;

% --------- DP sequential pairing (cost = |duration difference|; skip cost = skip_pen_ms)---------
DE = dur_eeg_ms(:);  NE = numel(DE);
DL = dur_lm_ms(:);   NL = numel(DL);

% Dynamic programming matrix
INF = 1e18;
C = INF*ones(NE+1, NL+1);
C(1,1) = 0;
BT = zeros(NE+1, NL+1, 'uint8');  % 1=diag(match), 2=up(skip EEG), 3=left(skip LM)

for i = 0:NE
    for j = 0:NL
        if i<NE && j<NL
            cost_match = C(i+1,j+1) + abs(DE(i+1) - DL(j+1));
            if cost_match < C(i+2,j+2)
                C(i+2,j+2) = cost_match; BT(i+2,j+2) = 1;
            end
        end
        if i<NE
            cost_skipE = C(i+1,j+1) + skip_pen_ms;
            if cost_skipE < C(i+2,j+1)
                C(i+2,j+1) = cost_skipE; BT(i+2,j+1) = 2;
            end
        end
        if j<NL
            cost_skipL = C(i+1,j+1) + skip_pen_ms;
            if cost_skipL < C(i+1,j+2)
                C(i+1,j+2) = cost_skipL; BT(i+1,j+2) = 3;
            end
        end
    end
end

% backtracking pair
i = NE; j = NL;
pairs = [];
while i>0 || j>0
    b = BT(i+1, j+1);
    if b==1
        pairs(end+1,:) = [i j]; %#ok<AGROW>
        i=i-1; j=j-1;
    elseif b==2
        i=i-1;
    elseif b==3
        j=j-1;
    else
        if i>0, i=i-1; elseif j>0, j=j-1; end
    end
end
pairs = flipud(pairs);

% ---------Generate output structures Pairs and summary---------
K = size(pairs,1);
Pairs = struct('eeg_seg_idx',{},'lm_trial',{},'dur_eeg_ms',{},'dur_lm_ms',{}, ...
               'dur_diff_ms',{},'match_ok',{},'strategy',{}, ...
               't_eeg',{},'eeg_data',{},'t_lm',{},'x_mm',{},'y_mm',{}, ...
               'ball_x_pix',{},'ball_color',{});
rows = cell(K,1);

for k = 1:K
    ii = pairs(k,1);  
    jj = pairs(k,2);  

    seg_i  = seg_idx_list(ii);  
    info_i = info_idx_list(ii); 
    lm_j   = jj;

    de = DE(ii);
    dl = DL(jj);
    diff_ms = dl - de;

    ok = abs(diff_ms) <= tol_ms;

    Pairs(end+1).eeg_seg_idx = info_i;                      %#ok<AGROW>
    Pairs(end  ).lm_trial    = lm_trial_id(lm_j);
    Pairs(end  ).dur_eeg_ms  = de;
    Pairs(end  ).dur_lm_ms   = dl;
    Pairs(end  ).dur_diff_ms = diff_ms;
    Pairs(end  ).match_ok    = ok;
    Pairs(end  ).strategy    = sprintf('dp_tol_%d_skip_%d', tol_ms, skip_pen_ms);

    Pairs(end  ).t_eeg       = E.times{seg_i};
    Pairs(end  ).eeg_data    = E.segs_blc{seg_i};
    Pairs(end  ).t_lm        = LM(lm_j).t;
    Pairs(end  ).x_mm        = LM(lm_j).x_mm;
    Pairs(end  ).y_mm        = LM(lm_j).y_mm;
    if isfield(LM,'ball_x_pix'),  Pairs(end).ball_x_pix  = LM(lm_j).ball_x_pix;  end
    if isfield(LM,'ball_color'),  Pairs(end).ball_color  = LM(lm_j).ball_color;  end

    rows{k} = {info_i, lm_trial_id(lm_j), de, dl, diff_ms, ok, Pairs(end).strategy};
end

% --------- Save ---------
if isempty(rows)
    T = cell2table(cell(0,7), 'VariableNames', ...
        {'eeg_seg_idx','lm_trial','dur_eeg_ms','dur_lm_ms','dur_diff_ms','match_ok','strategy'});
else
    T = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'eeg_seg_idx','lm_trial','dur_eeg_ms','dur_lm_ms','dur_diff_ms','match_ok','strategy'});
end

writetable(T, fullfile(out_dir,'eeg_leap_pairs_summary.csv'));
save(fullfile(out_dir,'eeg_leap_pairs.mat'), 'Pairs', 'tol_ms', 'skip_pen_ms', '-v7.3');

fprintf('[Align] EEG kept=%d, Leap kept=%d\n', numel(seg_idx_list), numel(LM));
fprintf('[Align] DP matches: %d pair(s), tol=%d ms, skip_pen=%d ms\n', size(pairs,1), tol_ms, skip_pen_ms);
fprintf('[Align] Wrote: eeg_leap_pairs.mat & eeg_leap_pairs_summary.csv\n');
end
