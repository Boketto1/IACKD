function process_leap_dual_csv(csv_ball, csv_hand, out_mat, fs_out, gap_ms, bp_band_hz, filt_order)
% Fuse left_ball_1 (trial, ball_x_pix, ball_color, move_direct, hit) + left_leap_1 (hand 3D)
% Output LM(k) with:
%   trial, t, global_time
%   x_mm(_raw), y_mm(_raw), z_mm(_raw)
%   ball_x_pix, ball_color (per-sample)
%   move_direct_seq (per-sample), move_direct_label (single label)
%   hit (per-sample logical), hit_time (abs sec or NaN)
%   duration_sec, valid

if nargin < 4 || isempty(fs_out),     fs_out = 100;      end
if nargin < 5 || isempty(gap_ms),     gap_ms  = 250;     end
if nargin < 6 || isempty(bp_band_hz), bp_band_hz = [0.02 4]; end
if nargin < 7 || isempty(filt_order), filt_order = 4;    end

Tb = readtable(csv_ball);
Th = readtable(csv_hand);
Tb.Properties.VariableNames = lower(strrep(strtrim(Tb.Properties.VariableNames),' ','_'));
Th.Properties.VariableNames = lower(strrep(strtrim(Th.Properties.VariableNames),' ','_'));

need_b = {'global_time','trial','ball_x_pix','ball_color'};
for i=1:numel(need_b), assert(ismember(need_b{i}, Tb.Properties.VariableNames), 'ball CSVmissing columns: %s', need_b{i}); end

% move_direct
md_candidates = {'move_direct','movedirect','move_direction','movement_direction','movement','direction','move_dir','dir'};
md_name = '';
for i=1:numel(md_candidates)
    if ismember(md_candidates{i}, Tb.Properties.VariableNames)
        md_name = md_candidates{i}; break;
    end
end
% hit
hit_candidates = {'hit','is_hit','contact','touched','touch','hit_flag'};
hit_name = '';
for i=1:numel(hit_candidates)
    if ismember(hit_candidates{i}, Tb.Properties.VariableNames)
        hit_name = hit_candidates{i}; break;
    end
end

% hand 3D
need_h = {'timestamp','x_mm','y_mm','z_mm'};
for i=1:numel(need_h), assert(ismember(need_h{i}, Th.Properties.VariableNames), 'hand CSVmissing columns: %s', need_h{i}); end

% time is converted from minutes to seconds
gt_ball = normalize_time_sec(Tb.global_time);
ts_hand = normalize_time_sec(Th.timestamp);

% Manually estimated linear mapping (please adjust here if necessary)
a = 3.040141; 
b = 0.999999894;
ts_hand = (ts_hand - a) / b;

% data column
trial_b = Tb.trial;
bx_b    = Tb.ball_x_pix;
bc_b    = string(Tb.ball_color);

% original move_direct
if ~isempty(md_name)
    md_raw = string(Tb.(md_name));
else
    md_raw = strings(height(Tb),1);
end
md_raw = strtrim(md_raw);
md_raw(md_raw=="" | ismissing(md_raw)) = "unknown";

% original hit
if ~isempty(hit_name)
    hit_raw = Tb.(hit_name);
    hit_bool = coerce_hit_to_logical(hit_raw);
else
    warning('The hit column was not found in the ball CSV file, so all trials are considered as no hits');
    hit_bool = false(height(Tb),1);
end

x_h = Th.x_mm;  y_h = Th.y_mm;  z_h = Th.z_mm;

% ==================Filter design: supports bandpass or lowpass==================
fs = fs_out;
fc = double(bp_band_hz(:)');           % Allow [lo hi] or scalar fc
if numel(fc)==2 && all(isfinite(fc)) && fc(1)>0
    % ------ band pass ------
    Wn_bp = fc/(fs/2);
    Wn_bp(1)=max(Wn_bp(1),1e-6); 
    Wn_bp(2)=min(Wn_bp(2),0.999);
    if Wn_bp(1) >= Wn_bp(2)
        Wn_lp = min(fc(2)/(fs/2),0.999);
        [b1,a1] = butter(filt_order, Wn_lp, 'low');
        mode_filt = 'low';
        Wn_for_fallback = Wn_lp; 
    else
        [b1,a1] = butter(filt_order, Wn_bp, 'bandpass');
        mode_filt = 'band';
        Wn_for_fallback = Wn_bp; 
    end
elseif numel(fc)==1 && isfinite(fc) && fc>0
    % ------ low-pass ------
    Wn_lp = min(fc/(fs/2),0.999);
    [b1,a1] = butter(filt_order, Wn_lp, 'low');
    mode_filt = 'low';
    Wn_for_fallback = Wn_lp;  
else
    error('bp_band_hz parameter is invalid: pass [lo hi] for bandpass, or pass scalar fc for lowpass.');
end
% ===================================================================

dt = 1/fs;
trials = unique(trial_b);
LM = struct([]);
for k = 1:numel(trials)
    mskb = (trial_b == trials(k)); if ~any(mskb), continue; end
    gt_k = gt_ball(mskb); bx_k = bx_b(mskb); bc_k = bc_b(mskb); md_k = md_raw(mskb); h_k = hit_bool(mskb);
    [gt_k,ord] = sort(gt_k);
    bx_k=bx_k(ord); bc_k=bc_k(ord); md_k=md_k(ord); h_k=h_k(ord);

    % Start and end: Start point = first frame; 
    % End point = first hit time (if no hit, then last frame)
    t0   = gt_k(1);
    hIdx = find(h_k,1,'first');
    if ~isempty(hIdx), t_end = gt_k(hIdx); else, t_end = gt_k(end); end
    dur_sec = t_end - t0; if ~(isfinite(dur_sec) && dur_sec>0), continue; end

    % Grid (absolute seconds)
    tgrid = (t0:dt:t_end)';

    % Nearest neighbor of ball/label to grid
    bx_r  = interp1(gt_k, bx_k, tgrid, 'nearest', 'extrap');
    bc_r  = interp1_categorical(gt_k, bc_k, tgrid);   % cellstr
    md_r  = interp1_categorical(gt_k, md_k, tgrid);   % cellstr

    % Hit-by-hit sampling (all are 1 from the first hit)
    hit_r = false(size(tgrid));
    if ~isempty(hIdx)
        hit_r(tgrid >= gt_k(hIdx)) = true;
        hit_time = gt_k(hIdx);
    else
        hit_time = NaN;
    end

    % Primary direction (majority voting)
    mdk_clean = md_k; mdk_clean(mdk_clean=="" | ismissing(mdk_clean)) = "unknown";
    md_label = char(mode(categorical(mdk_clean)));

    % Hand data window (+/-100 ms)
    pad = 0.100;
    inWin = (ts_hand >= (t0 - pad)) & (ts_hand <= (t_end + pad));
    ts_k  = ts_hand(inWin); xk = x_h(inWin); yk = y_h(inWin); zk = z_h(inWin);

    cover_ok = (~isempty(ts_k)) && (min(ts_k) <= t0 + 0.05) && (max(ts_k) >= t_end - 0.05);
    gaps_bad = false; if numel(ts_k)>=2, gaps_bad = any(diff(ts_k) > gap_ms/1000); end

    % First, perform linear interpolation to a uniform grid to obtain the raw trajectory
    if isempty(ts_k)
        x_raw = nan(size(tgrid)); y_raw = x_raw; z_raw = x_raw;
    else
        x_raw = interp1(ts_k, xk, tgrid, 'linear', 'extrap');
        y_raw = interp1(ts_k, yk, tgrid, 'linear', 'extrap');
        z_raw = interp1(ts_k, zk, tgrid, 'linear', 'extrap');
    end

    % ==================Robust filtering==================
    nfact1 = 3*(max(length(a1),length(b1))-1);
    minN1  = nfact1 + 1;
    Nsig   = numel(x_raw);

    if all(isfinite(x_raw)) && Nsig >= minN1
        x_f = filtfilt(b1,a1,double(x_raw));
        y_f = filtfilt(b1,a1,double(y_raw));
        z_f = filtfilt(b1,a1,double(z_raw));
    elseif all(isfinite(x_raw))
        if strcmp(mode_filt,'low')
            [b2,a2] = butter(2, Wn_for_fallback, 'low');
        else
            [b2,a2] = butter(2, Wn_for_fallback, 'bandpass');
        end
        nfact2 = 3*(max(length(a2),length(b2))-1);
        minN2  = nfact2 + 1;
        if Nsig >= minN2
            x_f = filtfilt(b2,a2,double(x_raw));
            y_f = filtfilt(b2,a2,double(y_raw));
            z_f = filtfilt(b2,a2,double(z_raw));
        else
            x_f = x_raw; y_f = y_raw; z_f = z_raw;
            warning('Trial too short for filtfilt (N=%d). Using unfiltered hand signals.', Nsig);
        end
    else
        x_f = x_raw; y_f = y_raw; z_f = z_raw;
    end
    % =============================================================================================

    % write LM
    LM(end+1).trial            = trials(k);                 %#ok<AGROW>
    LM(end  ).t                = tgrid - t0;
    LM(end  ).global_time      = tgrid;

    LM(end  ).x_mm_raw         = x_raw;  LM(end).y_mm_raw  = y_raw;  LM(end).z_mm_raw = z_raw;
    LM(end  ).x_mm             = x_f;    LM(end).y_mm      = y_f;    LM(end).z_mm     = z_f;

    LM(end  ).ball_x_pix       = bx_r;
    LM(end  ).ball_color       = bc_r;

    LM(end  ).move_direct_seq  = md_r;                      % per-sample
    LM(end  ).move_direct_label= md_label;                  % single label

    LM(end  ).hit              = hit_r;                     % per-sample 0/1
    LM(end  ).hit_time         = hit_time;                  % abs sec (NaN if none)

    LM(end  ).duration_sec     = dur_sec;
    LM(end  ).valid            = logical(cover_ok && ~gaps_bad);
end

% statistical printing
n_with_hit = sum(~isnan([LM.hit_time]));
fprintf('[LM] trials=%d, valid=%d, with_hit=%d\n', numel(LM), sum([LM.valid]), n_with_hit);

save(out_mat, 'LM', '-v7.3');
fprintf('[LM] Saved: %s\n', out_mat);
end

% ---------- helpers ----------
function t = normalize_time_sec(t)
t = double(t(:)); mx = max(t);
if mx > 1e8, t = t/1e6; elseif mx > 1e5, t = t/1e3; end
end
function v = coerce_hit_to_logical(x)
if islogical(x), v = x(:);
elseif isnumeric(x), v = x(:) ~= 0;
elseif iscellstr(x) || isstring(x)
    s = lower(strtrim(string(x(:))));
    v = ismember(s, ["1","true","t","yes","y"]);
else, error('The data type of the "hit" column cannot be identified');
end
v = logical(v(:));
end
function C = interp1_categorical(t_src, s_src, t_dst)
s = string(s_src(:));
s = strtrim(s); s(s=="" | ismissing(s)) = "unknown";
idx = interp1(t_src(:), (1:numel(s))', t_dst(:), 'nearest', 'extrap');
C = cellstr(s(idx));
end
