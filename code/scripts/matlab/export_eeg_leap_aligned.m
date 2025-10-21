function export_eeg_leap_aligned(eeg_seg_mat, lm_mat, pairs_mat, calib_mat, out_mat, opts)
% Map LEAP to EEG timeline (-200ms -> 66). No "hit" required.
% Adds move_direct as a single string per trial (majority over in-range).
if nargin<6, opts = struct; end
if ~isfield(opts,'tol_ms'),   opts.tol_ms   = 600;  end   % only check start residual now
if ~isfield(opts,'min_cov'),  opts.min_cov  = 0.75; end
if ~isfield(opts,'keep_all'), opts.keep_all = false;end
if ~isfield(opts,'verbose'),  opts.verbose  = true; end

E = load(eeg_seg_mat);   L = load(lm_mat);  P = load(pairs_mat);  C = load(calib_mat);
assert(isfield(E,'trial_info') && isfield(E,'RESAMPLE_HZ') && isfield(E,'times') && isfield(E,'segs_blc'));
assert(isfield(L,'LM'));  assert(isfield(P,'Pairs') && ~isempty(P.Pairs));
assert(isfield(C,'a') && isfield(C,'b'));

a = C.a; b = C.b; fs = E.RESAMPLE_HZ;
Pairs = P.Pairs(:); LM = L.LM(:); lm_id = arrayfun(@(u) u.trial, LM);

OUT = struct([]); kept=0; skipped=0;

for k = 1:numel(Pairs)
    ei  = Pairs(k).eeg_seg_idx;
    tid = Pairs(k).lm_trial;
    li  = find(lm_id==tid, 1);
    if isempty(li)
        skipped = skipped+1; if opts.verbose, fprintf('[SKIP] LM trial %d not found.\n', tid); end
        continue;
    end

    % EEG -> Absolute seconds -> Mapping to LM global time
    t14_abs   = E.trial_info(ei).t14_samp / fs;
    t_eeg_rel = Pairs(k).t_eeg(:);
    t_eeg_abs = t14_abs + t_eeg_rel;
    t_target  = a + b * t_eeg_abs;

    % LM support range
    gt = LM(li).global_time(:); tmin=gt(1); tmax=gt(end);
    inrange = (t_target >= tmin) & (t_target <= tmax);
    cov = mean(inrange);

    % Hand 3D data (interpolated only for inrange; the rest are NaN)
    x = LM(li).x_mm(:); y = LM(li).y_mm(:);
    haveZ = isfield(LM,'z_mm') && ~isempty(LM(li).z_mm);
    if haveZ, z = LM(li).z_mm(:); else, z = nan(size(x)); end
    bx = LM(li).ball_x_pix(:);

    xa=nan(size(t_target)); ya=xa; za=xa; bxa=xa;
    idx=find(inrange);
    if ~isempty(idx)
        xa(idx)  = interp1(gt, x,  t_target(idx), 'linear');
        ya(idx)  = interp1(gt, y,  t_target(idx), 'linear');
        if haveZ, za(idx)= interp1(gt, z,  t_target(idx), 'linear'); end
        bxa(idx) = interp1(gt, bx, t_target(idx), 'nearest');
    end

    % Color and Direction
    bc_trial = majority_label_over_range(LM(li), 'ball_color', gt, t_target, inrange);
    md_trial = majority_label_over_range(LM(li), 'move_direct', gt, t_target, inrange);

    % Residual (only starting point)
    t0_lm = gt(1);
    start_res_ms = ((a + b*t14_abs) - t0_lm) * 1000;

    ok_res = (abs(start_res_ms) <= opts.tol_ms);
    ok_cov = (cov >= opts.min_cov);
    ok = ok_res & ok_cov;

    if ~ok && ~opts.keep_all
        skipped = skipped+1;
        if opts.verbose
            fprintf('[DROP] ei=%d tid=%d  cov=%.2f  start=%.1fms  (no end_res)\n', ei, tid, cov, start_res_ms);
        end
        continue;
    end

    % write OUT
    OUT(end+1).eeg_seg_idx = ei;                        %#ok<AGROW>
    OUT(end  ).lm_trial    = tid;
    OUT(end  ).t_ms        = round(t_eeg_rel*1000);
    OUT(end  ).EEG         = Pairs(k).eeg_data;

    OUT(end  ).x_mm        = xa(:).';
    OUT(end  ).y_mm        = ya(:).';
    OUT(end  ).z_mm        = za(:).';

    OUT(end  ).ball_x_pix  = bxa(:).';
    OUT(end  ).ball_color  = bc_trial;  
    OUT(end  ).move_direct = md_trial;   

    OUT(end  ).a           = a;
    OUT(end  ).b           = b;
    OUT(end  ).start_res_ms= start_res_ms;
    OUT(end  ).cov         = cov;

    kept = kept+1;
end

save(out_mat, 'OUT', '-v7.3');
fprintf('[EXPORT] kept=%d, dropped=%d  (tol=%d ms, min_cov=%.2f) -> %s\n', kept, skipped, opts.tol_ms, opts.min_cov, out_mat);
end

% ---- helpers ----
function s = majority_label_over_range(LM, field, gt, t_target, inrange)
% Robust: try in-range voting; if empty -> use *_label fallback; final fallback 'unknown'
% field: 'ball_color' or 'move_direct'
fallback = '';
lab_field = field;
if strcmp(field,'move_direct')
    lab_field = 'move_direct_label';
end
if isfield(LM, lab_field) && ~isempty(LM.(lab_field))
    if iscell(LM.(lab_field))
        fallback = char(LM.(lab_field){1});
    else
        fallback = char(LM.(lab_field));
    end
end

if ~isfield(LM, field) || isempty(LM.(field))
    s = iif_empty(fallback, 'unknown'); return;
end

src = LM.(field);
if iscell(src), src = string(src(:)); else, src = string(src(:)); end
src = strtrim(src); src(src=="" | ismissing(src)) = "unknown";

if any(inrange)
    idx_nn = interp1(gt, (1:numel(src))', t_target(inrange), 'nearest');
    vals = src(idx_nn);
else
    vals = string.empty(0,1);
end

if isempty(vals)
    s = iif_empty(fallback, 'unknown');
else
    s = char(mode(categorical(vals)));
end
end

function out = iif_empty(s, alt)
if isempty(s)
    out = alt;
else
    out = s;
end
end
