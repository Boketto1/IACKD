function calibrate_pairs_with_events(eeg_seg_mat, lm_mat, pairs_mat, out_dir, tol_ms, opts)
% calibration: global_time ≈ a + b * eeg_time
% Anchor points: EEG 14 ↔ LM first frame (all pairs); EEG 66 ↔ LM first hit_time (only pairs with hits)
% Robustness: First, calculate b0 using the duration ratio with hits, then perform LS fitting and iterate to remove outliers
if nargin < 5 || isempty(tol_ms), tol_ms = 600; end
if nargin < 6, opts = struct; end
if ~isfield(opts,'max_outlier_ms'), opts.max_outlier_ms = 1500; end
if ~isfield(opts,'max_iter'),       opts.max_iter       = 3;    end
if ~isfield(opts,'min_pairs'),      opts.min_pairs      = 8;    end
if ~exist(out_dir,'dir'), mkdir(out_dir); end

E = load(eeg_seg_mat);
L = load(lm_mat);
P = load(pairs_mat);
assert(isfield(E,'trial_info') && isfield(E,'RESAMPLE_HZ'), 'EEG 缺 trial_info/RESAMPLE_HZ');
assert(isfield(L,'LM'), 'LM 缺 LM 结构');
assert(isfield(P,'Pairs') && ~isempty(P.Pairs), 'Pairs 为空');

fs    = double(E.RESAMPLE_HZ);
Pairs = P.Pairs(:);
LM    = L.LM(:);
lm_id = arrayfun(@(u) u.trial, LM);

% collect anchor points
X = []; Y = []; tag = [];            % tag: 1=start, 2=end(with hit)
dur_eeg = []; dur_lm = [];           % Only collect data with hits for b0
n_with_hit = 0;

for k = 1:numel(Pairs)
    ei  = Pairs(k).eeg_seg_idx;
    tid = Pairs(k).lm_trial;
    li  = find(lm_id==tid, 1);
    if isempty(li), continue; end

    % EEG anchors
    t14 = double(E.trial_info(ei).t14_samp)/fs;
    t66 = double(E.trial_info(ei).tend_samp)/fs;

    % LM start & hit
    gt = double(LM(li).global_time(:));
    if isempty(gt), continue; end
    t0 = gt(1);

    thit = NaN;
    if isfield(LM,'hit_time') && ~isempty(LM(li).hit_time)
        thit = double(LM(li).hit_time);
    elseif isfield(LM,'hit') && ~isempty(LM(li).hit)
        h = LM(li).hit(:); hi = find(h~=0,1,'first'); if ~isempty(hi), thit = gt(hi); end
    end

    % Starting anchor point
    if isfinite(t14) && isfinite(t0)
        X(end+1,1) = t14; Y(end+1,1) = t0; tag(end+1,1) = 1; %#ok<AGROW>
    end

    % Ending anchor point 
    if isfinite(t66) && isfinite(thit)
        X(end+1,1) = t66; Y(end+1,1) = thit; tag(end+1,1) = 2; %#ok<AGROW>
        dur_eeg(end+1,1) = (t66 - t14);                           %#ok<AGROW>
        dur_lm(end+1,1)  = (thit - t0);                           %#ok<AGROW>
        n_with_hit = n_with_hit + 1;
    end
end

assert(numel(X) >= opts.min_pairs, 'Insufficient event anchor points available for calibration (%d<%d)。', numel(X), opts.min_pairs);

ratio = dur_lm ./ dur_eeg; ratio = ratio(isfinite(ratio) & ratio>0);
if isempty(ratio), b0 = 1.0; else, b0 = median(ratio); end
a0 = median( Y(tag==1) - b0 * X(tag==1) ); 

% LS intends to merge and iteratively remove outliers
a = a0; b = b0;
for it = 1:opts.max_iter
    res_ms = (Y - (a + b*X)) * 1000;
    inlier = abs(res_ms) <= opts.max_outlier_ms;
    if sum(inlier) < opts.min_pairs, break; end
    A = [ones(sum(inlier),1) X(inlier)];
    coef = A \ Y(inlier);
    a = coef(1); b = coef(2);
    if all(inlier), break; end
end

res_ms   = (Y - (a + b*X)) * 1000;
med_res  = median(abs(res_ms),'omitnan');
in_tol   = mean(abs(res_ms) <= tol_ms);

stats = struct();
stats.n_pairs_total    = numel(Pairs);
stats.n_anchors        = numel(X);
stats.n_with_hit_pairs = n_with_hit;
stats.med_abs_res_ms   = med_res;
stats.in_tol_rate      = in_tol;
stats.tol_ms           = tol_ms;
stats.max_outlier_ms   = opts.max_outlier_ms;

out_mat = fullfile(out_dir, 'eeg_leap_pairs_events_from_dp.mat');
save(out_mat, 'a','b','stats');
fprintf('[CAL-HIT] global ≈ %.6f + %.9f * eeg  (|res|_med≈%.1f ms; anchors=%d; with_hit_pairs=%d)\n', ...
    a, b, med_res, numel(X), n_with_hit);
end
