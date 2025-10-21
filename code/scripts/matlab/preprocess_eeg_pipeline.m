function preprocess_premotor_55_14(CDT_FILE, CSV_BALL, OUT_MAT, opts)
% 预处理 + 提取预运动(55→14)片段 + 读取 CSV 标签（color_ball / move_direct）
% 仅保存一个你指定路径的 .mat 文件（OUT_MAT）
%
% 必选入参：
%   CDT_FILE  : .cdt 文件路径
%   CSV_BALL  : left_ball_1.csv 路径（包含 color_ball / move_direct / trial / global_time）
%   OUT_MAT   : 输出的 .mat 路径（若目录不存在会自动创建）
%
% 可选参数(opts)（有默认）：
%   .resample_hz = 200
%   .bp_lo_hz    = 0.1
%   .bp_hi_hz    = 50
%   .del_chans   = {'HEO','VEOG','TRIGGER'}
%   .iclabel_thr = 0.90      % Eye/Muscle ≥此概率才移除
%   .reref_mode  = 'average'
%   .pre_ms      = 1000      % 固定段长度（以14为终点的 1s）
%   .tol_ms_lo   = 800       % 接受 55→14 时长下限
%   .tol_ms_hi   = 1200      % 接受 55→14 时长上限

if nargin < 4, opts = struct; end
opts = setdef(opts, 'resample_hz', 200);
opts = setdef(opts, 'bp_lo_hz',    0.1);
opts = setdef(opts, 'bp_hi_hz',    45);
opts = setdef(opts, 'del_chans',   {'M1','M2','HEOG','VEOG','TRIGGER'});
opts = setdef(opts, 'iclabel_thr', 0.90);
opts = setdef(opts, 'reref_mode',  'average');
opts = setdef(opts, 'pre_ms',      1000);
opts = setdef(opts, 'tol_ms_lo',   800);
opts = setdef(opts, 'tol_ms_hi',   1200);

% 确保输出目录存在
outdir = fileparts(OUT_MAT);
if ~isempty(outdir) && ~exist(outdir,'dir'), mkdir(outdir); end

% ---------------- 0) 载入 EEG ----------------
fprintf('[0] Loading: %s\n', CDT_FILE);
[~,~,ext] = fileparts(CDT_FILE);
switch lower(ext)
    case '.set', EEG = pop_loadset(CDT_FILE);
    case {'.cdt','.dap','.dat','.rs3'}
        EEG = pop_loadcurry(CDT_FILE, 'CurryLocations','on', 'CurryEvents','on');
    otherwise, error('Unsupported EEG file: %s', ext);
end
EEG = eeg_checkset(EEG);

% ---------------- 1) 电极位置 -----------------
fprintf('[1] Set channel locations...\n');
EEG = try_set_chanlocs(EEG, '');

% ---------------- 2) 降采样 200 Hz -------------
fprintf('[2] Resample to %d Hz...\n', opts.resample_hz);
EEG = pop_resample(EEG, opts.resample_hz);

% ---------------- 3) 0.1–50 Hz 带通 ------------
fprintf('[3] Bandpass %.1f–%.1f Hz...\n', opts.bp_lo_hz, opts.bp_hi_hz);
EEG = pop_eegfiltnew(EEG, opts.bp_lo_hz, opts.bp_hi_hz);

% ---------------- 4) 删通道 --------------------
fprintf('[4] Drop chans: %s\n', strjoin(opts.del_chans, ', '));
EEG = delete_if_present(EEG, opts.del_chans);

% ---------------- 5) ICA+ICLabel（Eye/Muscle）--
fprintf('[5] ICA + ICLabel (Eye/Muscle >= %.2f)\n', opts.iclabel_thr);
EEG = run_ica_remove_eye_muscle(EEG, opts.iclabel_thr);

% ---------------- 6) 重参考 --------------------
fprintf('[6] Re-reference: %s\n', opts.reref_mode);
EEG = do_reref(EEG, opts.reref_mode);

% ---------------- 7) 提取 55→14 片段 ------------
fprintf('[7] Extract 55→14 (var) & [-%d..0]ms (1s) segments...\n', opts.pre_ms);
[segs_var, times_var, segs_1s, times_1s, trial_info, keep_mask_1s] = ...
    extract_premotor_segments(EEG, 55, 14, opts.pre_ms, [opts.tol_ms_lo opts.tol_ms_hi]);
% —— 统一 1s 片段长度为 200 点（C×200）
[segs_1s, times_1s] = force_len_fixed(segs_1s, times_1s, 200, EEG.nbchan, EEG.srate);

% ---------------- 8) 读 CSV 标签 ----------------
fprintf('[8] Read labels from CSV: %s\n', CSV_BALL);
labels = read_ball_labels(CSV_BALL);

% 以“出现顺序”对齐标签与 EEG 段
n_eeg = numel(segs_var);
n_lab = height(labels);
m = min(n_eeg, n_lab);
labels = labels(1:m, :);
if n_eeg > m
    pad = table( (m+1:n_eeg)', repmat("",n_eeg-m,1), repmat("",n_eeg-m,1), ...
        'VariableNames', {'trial_csv','color_ball','move_direct'});
    labels = [labels; pad];
end

% ---------------- 保存到你指定的 .mat ----------------
meta = struct();
meta.srate       = EEG.srate;
meta.chanlabels  = {EEG.chanlocs.labels}';
meta.nbchan      = EEG.nbchan;
meta.preprocess  = opts;

save(OUT_MAT, 'segs_var','times_var','segs_1s','times_1s', ...
    'trial_info','keep_mask_1s','labels','meta','-v7.3');
fprintf('Saved: %s\n', OUT_MAT);
fprintf('Done.\n');
end

% ====== 辅助函数 ======
function s = setdef(s, f, v), if ~isfield(s, f), s.(f) = v; end, end

function EEG = try_set_chanlocs(EEG, CHANLOC_FILE)
try
    if ~isempty(CHANLOC_FILE) && exist(CHANLOC_FILE,'file')
        EEG = pop_chanedit(EEG, 'lookup', CHANLOC_FILE); return;
    end
    eeglabroot = fileparts(which('eeglab.m'));
    cands = { ...
        fullfile(eeglabroot,'plugins','dipfit','standard_BEM','elec','standard_1005.elc'), ...
        fullfile(eeglabroot,'plugins','dipfit','standard_BEM','elec','standard_1020.elc'), ...
        fullfile(eeglabroot,'sample_locs','standard-10-5-cap385.elp')};
    for i=1:numel(cands)
        if exist(cands{i},'file'), EEG = pop_chanedit(EEG, 'lookup', cands{i}); return; end
    end
catch ME
    warning('Set chanlocs failed: %s', ME.message);
end
end

function EEG = delete_if_present(EEG, chans)
labels = {EEG.chanlocs.labels};
mask = ismember(labels, chans);
if any(mask)
    EEG = pop_select(EEG, 'nochannel', labels(mask));
end
end

function EEG = run_ica_remove_eye_muscle(EEG, thr)
EEG = pop_runica(EEG, 'icatype','runica','extended',1,'interrupt','off');
if exist('pop_iclabel','file')==2
    EEG = pop_iclabel(EEG, 'default');
    if isfield(EEG,'etc') && isfield(EEG.etc,'ic_classification') && ...
       isfield(EEG.etc.ic_classification,'ICLabel')
        P = EEG.etc.ic_classification.ICLabel.classifications;  % [IC x 7]
        % 列: Brain(1) Muscle(2) Eye(3) Heart(4) Line(5) Chan(6) Other(7)
        rm = find(P(:,2)>=thr | P(:,3)>=thr);   % 仅删肌电&眼电
        if ~isempty(rm)
            fprintf('    Removing %d ICs (Eye/Muscle >= %.2f)\n', numel(rm), thr);
            EEG = pop_subcomp(EEG, rm, 0);
        else
            fprintf('    No IC meets Eye/Muscle >= %.2f\n', thr);
        end
    else
        warning('ICLabel result not found; skip removal.');
    end
else
    warning('ICLabel plugin not found; skip IC removal.');
end
end

function EEG = do_reref(EEG, mode)
try
    if ischar(mode) && strcmpi(mode,'average')
        EEG = pop_reref(EEG, []);
    elseif iscell(mode)
        lab = {EEG.chanlocs.labels};
        idx = find(ismember(lab, mode)); 
        if isempty(idx), warning('Ref chans missing; using average.'); EEG = pop_reref(EEG, []);
        else, EEG = pop_reref(EEG, idx); end
    else
        EEG = pop_reref(EEG, []);
    end
catch ME
    warning('Reref failed: %s; using average.', ME.message);
    EEG = pop_reref(EEG, []);
end
end

function [segs_var, times_var, segs_1s, times_1s, info, keep_1s] = ...
    extract_premotor_segments(EEG, code_start, code_end, pre_ms, tol_ms_pair)
% 输出：
%   segs_var/times_var : 55→14 变长段（0 对齐 14；时间单位秒）
%   segs_1s/times_1s   : 固定 1s（[-1,0]s，0 对齐 14）
%   info               : 结构体数组（i55/i14/t55/t14/dur_ms/valid）
%   keep_1s            : 满足 dur 在 [tol_lo, tol_hi]ms 的 trial 掩码

srate = EEG.srate;
types = arrayfun(@(e) type2num(e.type), EEG.event);
lats  = [EEG.event.latency];                % 样本点
idx55 = find(types==code_start);            % 55
idx14 = find(types==code_end);              % 14

segs_var = {}; times_var = {};
segs_1s  = {}; times_1s  = {};
info = struct('i55',{},'i14',{},'t55',{},'t14',{},'dur_ms',{},'valid',{});
keep_1s = [];

pre_samp = round(pre_ms/1000 * srate);

for ii = 1:numel(idx14)
    i14 = idx14(ii);
    t14 = lats(i14);

    % “相邻的上一枚 55”：在上一个 14 之后且本次 14 之前的最后一个 55
    prev14_lat = -inf; if ii>1, prev14_lat = lats(idx14(ii-1)); end
    cand = idx55(lats(idx55) > prev14_lat & lats(idx55) < t14);
    i55 = NaN; t55 = NaN; valid = false; dur_ms = NaN;
    if ~isempty(cand)
        i55 = cand(end); t55 = lats(i55);
        dur_ms = (t14 - t55)/srate*1000;
        valid = (dur_ms>=0);
    end

    % 变长 55→14
    if valid
        s1 = max(1, round(t55));
        s2 = min(size(EEG.data,2), round(t14));
        segs_var{end+1}  = EEG.data(:, s1:s2);
        t0               = ((s1:s2) - t14)/srate;   % 0 对齐 14
        times_var{end+1} = t0;
    else
        segs_var{end+1}  = [];
        times_var{end+1} = [];
    end

    % 固定 1s（仅当 dur 在容差内）
    ok1s = valid && (dur_ms >= tol_ms_pair(1)) && (dur_ms <= tol_ms_pair(2));
    if ok1s
        s1f = round(t14 - pre_samp);
        s2f = round(t14);
        if s1f >= 1 && s2f <= size(EEG.data,2)
            segs_1s{end+1}  = EEG.data(:, s1f:s2f);
            times_1s{end+1} = ((s1f:s2f) - t14)/srate;
        else
            segs_1s{end+1}  = [];
            times_1s{end+1} = [];
            ok1s = false;
        end
    else
        segs_1s{end+1}  = [];
        times_1s{end+1} = [];
    end

    info(end+1).i55     = i55;     %#ok<AGROW>
    info(end  ).i14     = i14;
    info(end  ).t55     = t55;
    info(end  ).t14     = t14;
    info(end  ).dur_ms  = dur_ms;
    info(end  ).valid   = valid;
    keep_1s(end+1,1)    = ok1s;    %#ok<AGROW>
end
end

function v = type2num(x)
if isnumeric(x), v = x; return; end
if ischar(x) || isstring(x)
    y = str2double(x); if ~isnan(y), v = y; return; end
end
v = NaN;
end

function T = read_ball_labels(csvfile)
Tb = readtable(csvfile);
Tb.Properties.VariableNames = lower(strrep(strtrim(Tb.Properties.VariableNames),' ','_'));
allnames = Tb.Properties.VariableNames;

trial_col = pick_col(allnames, {'trial','trial_id','trialindex','trial_idx'});
assert(~isempty(trial_col), 'CSV 缺少 trial 列');
color_col = pick_col(allnames, {'color_ball','ball_color','color','colour'});
assert(~isempty(color_col), 'CSV 缺少 color_ball / ball_color 列');
move_col  = pick_col(allnames, {'move_direct','movedirect','move_direction','movement_direction','movement','direction','move_dir','dir'});
assert(~isempty(move_col), 'CSV 缺少 move_direct 列');

% —— 清洗：去掉首尾引号、空白并小写 —— 
clean = @(s) regexprep(lower(strtrim(string(s))), '^["'']+|["'']+$', '');

trial = Tb.(trial_col);
color = clean(Tb.(color_col));
move  = clean(Tb.(move_col));

u = unique(trial, 'stable');
trial_csv  = zeros(numel(u),1);
color_lab  = strings(numel(u),1);
move_lab   = strings(numel(u),1);
label_code = -1*ones(numel(u),1,'int8');

for k = 1:numel(u)
    m = (trial == u(k));
    c_k = char(mode(categorical(color(m))));
    d_k = char(mode(categorical(move(m))));
    color_lab(k) = string(c_k);
    move_lab(k)  = string(d_k);
    trial_csv(k) = u(k);
    label_code(k)= map_label_code(d_k, c_k);  % 映射到 0/1/2/3/-1
end

T = table(trial_csv, color_lab, move_lab, label_code, ...
          'VariableNames', {'trial_csv','color_ball','move_direct','label_code'});

    function name = pick_col(names, aliases)
        name = '';
        for ii = 1:numel(aliases)
            hit = find(strcmpi(names, aliases{ii}), 1, 'first');
            if ~isempty(hit), name = names{hit}; return; end
        end
    end

    function code = map_label_code(move_dir, color_ball)
        md = clean(move_dir); cb = clean(color_ball);
        if md=="left"  && cb=="red",    code = int8(0); return; end
        if md=="right" && cb=="red",    code = int8(1); return; end
        if md=="left"  && cb=="yellow", code = int8(2); return; end
        if md=="right" && cb=="yellow", code = int8(3); return; end
        code = int8(-1);
    end
end
function [Sfix, Tfix] = force_len_fixed(Sin, Tin, target_len, nbchan, srate)
% 将 cell 数组中的每段 EEG 统一成 C×target_len：
% - 末尾裁剪（T>target_len）
% - 末尾 0 填充（T<target_len）
% 同步裁剪/补齐 times（线性延拓），单位秒
n = numel(Sin);
Sfix = cell(n,1);
Tfix = cell(n,1);

for i = 1:n
    Xi = Sin{i};
    ti = Tin{i};

    if isempty(Xi)
        Sfix{i} = zeros(nbchan, target_len);
        if ~isempty(srate) && isfinite(srate)
            dt = 1./double(srate);
            Tfix{i} = ((-target_len+1):0) * dt;   % 长度=target_len，结束于 0
        else
            Tfix{i} = zeros(1, target_len);
        end
        continue;
    end

    [C,T] = size(Xi);
    if C ~= nbchan
        % 若通道数异常，尽量适配当前段的通道数
        nb = C;
    else
        nb = nbchan;
    end

    % —— 数据裁剪/填充 —— 
    if T >= target_len
        Sfix{i} = Xi(:, 1:target_len);            % 裁掉末尾（最靠近 0ms 的点）
        t_out    = safe_times_trim(ti, target_len);
    else
        Xo = zeros(nb, target_len);
        Xo(:,1:T) = Xi;
        Sfix{i} = Xo;
        t_out    = safe_times_pad(ti, target_len, srate);
    end
    Tfix{i} = t_out(:)';  % 行向量
end
end

function t_out = safe_times_trim(t_in, target_len)
% 裁剪 times 到前 target_len 个样本；若 times 为空则回退为 0 序列
if isempty(t_in)
    t_out = zeros(1, target_len);
else
    t_in = t_in(:)';
    if numel(t_in) >= target_len
        t_out = t_in(1:target_len);
    else
        % 不该发生，但兜底：不足则线性外推
        dt = median(diff(t_in));
        last = t_in(end);
        need = target_len - numel(t_in);
        t_extra = last + dt*(1:need);
        t_out = [t_in, t_extra];
    end
end
end

function t_out = safe_times_pad(t_in, target_len, srate)
% 将 times 末尾按等步长补足到 target_len
if isempty(t_in)
    if ~isempty(srate) && isfinite(srate)
        dt = 1./double(srate);
        t_out = ((-target_len+1):0) * dt;
    else
        t_out = zeros(1, target_len);
    end
    return;
end
t_in = t_in(:)';
T = numel(t_in);
if T == 1
    dt = (1./double(srate));
else
    dt = median(diff(t_in));
    if ~isfinite(dt) || dt<=0
        dt = (1./double(srate));
    end
end
need = target_len - T;
if need <= 0
    t_out = t_in(1:target_len);
else
    t_extra = t_in(end) + dt*(1:need);
    t_out = [t_in, t_extra];
end
end
