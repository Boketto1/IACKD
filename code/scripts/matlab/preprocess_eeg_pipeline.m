function preprocess_eeg_pipeline()
% ================== CONFIG (edit to your paths) ==================
INPUT_FILE    = 'F:\SC_data\a\sub1 zdb\Acquisition 08l.cdt';
OUTPUT_DIR    = 'D:\EEG_data\SC_data\s02';
CHANLOC_FILE  = '';   % optional custom loc file, leave empty to use EEGLAB templates
DEL_CHANS     = {'HEO','VEO','TRIGGER','HEOG','VEOG'};

RESAMPLE_HZ   = 100;         % resample target Hz
BP_LO_HZ      = 0.1;         % bandpass low cut
BP_HI_HZ      = 30;          % bandpass high cut
ICA_EYE_PTHR  = 0.90;        % remove ICs with Eye prob >= this
REREF_MODE    = 'average';   % 'average' or {'M1','M2'}

START_CODE    = 14;          % event code start
END_CODE      = 1000001;     % event code end (=66)
BASE_MS       = [-200 0];    % baseline window relative to 14 (ms)
MAX_GAP_MS    = 8000;        % max allowed 14->END (ms). use inf to disable
REQUIRE_FULL_BASELINE = true;

AMP_THRESH    = 350;         % absolute amplitude threshold (uV)
% ================================================================

if ~exist(OUTPUT_DIR,'dir'), mkdir(OUTPUT_DIR); end

%% [0] Load
fprintf('[0] Loading dataset: %s\n', INPUT_FILE);
[~,~,ext] = fileparts(INPUT_FILE);
switch lower(ext)
    case '.set'
        EEG = pop_loadset(INPUT_FILE);
    case {'.cdt', '.dap', '.dat', '.rs3'}
        EEG = pop_loadcurry(INPUT_FILE, 'CurryLocations','on', 'CurryEvents','on');
    otherwise
        error('Unsupported input file: %s', ext);
end
EEG = eeg_checkset(EEG);

%% [1] Channel locations
fprintf('[1] Setting channel locations...\n');
EEG = try_set_chanlocs(EEG, CHANLOC_FILE);
EEG = eeg_checkset(EEG);

%% [2] Resample
fprintf('[2] Resampling to %d Hz...\n', RESAMPLE_HZ);
EEG = pop_resample(EEG, RESAMPLE_HZ);
EEG = eeg_checkset(EEG);

%% [3] Bandpass
fprintf('[3] Bandpass %.1f-%.1f Hz (zero-phase FIR)...\n', BP_LO_HZ, BP_HI_HZ);
EEG = pop_eegfiltnew(EEG, BP_LO_HZ, BP_HI_HZ);
EEG = eeg_checkset(EEG);

%% [4] Delete channels
fprintf('[4] Deleting channels if present: %s\n', strjoin(DEL_CHANS, ', '));
EEG = delete_if_present(EEG, DEL_CHANS);
EEG = eeg_checkset(EEG);

%% [5] ICA + ICLabel (remove ocular only)
fprintf('[5] ICA + ICLabel (remove ocular >= %.2f)...\n', ICA_EYE_PTHR);
EEG = run_ica_and_remove_eyes(EEG, ICA_EYE_PTHR);
EEG = eeg_checkset(EEG);

%% [6] Re-reference
fprintf('[6] Re-referencing (%s)...\n', REREF_MODE);
EEG = do_reref(EEG, REREF_MODE);
EEG = eeg_checkset(EEG);

%% Save continuous preprocessed set
preproc_set = fullfile(OUTPUT_DIR, 'EEG_preproc.set');
EEG.setname = 'EEG_preproc';
pop_saveset(EEG, 'filename', preproc_set);
fprintf('    Saved: %s\n', preproc_set);

%% [7] Extract segments: adjacent start (14) to end (END_CODE), with pre baseline
fprintf('[7] Extract segments (adjacent 14 -> %d, with %d ms pre)...\n', END_CODE, abs(BASE_MS(1)));
[segs, times, trial_info] = extract_segments_varlen(EEG, START_CODE, END_CODE, BASE_MS(1), MAX_GAP_MS, REQUIRE_FULL_BASELINE);

% map: which trial_info rows produced segments (their indices, in order)
valid_info_idx = find([trial_info.valid] == 1);
ntr = numel(segs);
if numel(valid_info_idx) ~= ntr
    warning('Consistency: valid info entries (%d) != number of segments (%d).', numel(valid_info_idx), ntr);
end

%% [8] Baseline correction
fprintf('[8] Baseline correcting [%d %d] ms rel. to 14...\n', BASE_MS(1), BASE_MS(2));
[segs_blc, reject_baseline] = baseline_correct_segments(segs, times, BASE_MS);

%% [9] Absolute amplitude reject (any sample exceeds +/- threshold)
fprintf('[9] Reject if ANY sample in whole segment exceeds +/- %d uV...\n', AMP_THRESH);
ntr = numel(segs_blc);
max_abs_all = zeros(ntr,1);
for i = 1:ntr
    Xi = segs_blc{i};
    if isempty(Xi)
        max_abs_all(i) = Inf;
    else
        max_abs_all(i) = max(abs(Xi(:)));
    end
end
keep_amp  = (max_abs_all <= AMP_THRESH);
keep_base = (~reject_baseline(:));

% final keep at segment level (length = ntr)
keep_mask_segments = keep_amp(:) & keep_base(:);

% project back to info level (length = numel(trial_info))
keep_mask_info = false(numel(trial_info),1);
m = min(ntr, numel(valid_info_idx));
if m > 0
    keep_mask_info(valid_info_idx(1:m)) = keep_mask_segments(1:m);
end

fprintf('    Segments: %d | kept by amp: %d | kept by baseline: %d | FINAL kept: %d\n', ...
    ntr, sum(keep_amp), sum(keep_base), sum(keep_mask_segments));

%% Save segments and metadata
out_mat = fullfile(OUTPUT_DIR, 'segments_14_to_1000001.mat');
save(out_mat, 'segs','times','segs_blc','trial_info', ...
              'keep_mask_segments','keep_mask_info', ...
              'AMP_THRESH','BASE_MS','START_CODE','END_CODE','RESAMPLE_HZ','-v7.3');
fprintf('    Saved segments: %s\n', out_mat);

%% Export CSV aligned to trial_info (one row per END)
out_csv = fullfile(OUTPUT_DIR, 'valid_trials.csv');
T = make_info_level_table(trial_info, keep_mask_info);
writetable(T, out_csv);
fprintf('    Saved trial list: %s\n', out_csv);

fprintf('All done.\n');
end

% ================= helpers =================

function EEG = try_set_chanlocs(EEG, CHANLOC_FILE)
try
    if ~isempty(CHANLOC_FILE) && exist(CHANLOC_FILE,'file')
        EEG = pop_chanedit(EEG, 'lookup', CHANLOC_FILE);
        return;
    end
    eeglabroot = fileparts(which('eeglab.m'));
    candidates = { ...
        fullfile(eeglabroot,'plugins','dipfit','standard_BEM','elec','standard_1005.elc'), ...
        fullfile(eeglabroot,'plugins','dipfit','standard_BEM','elec','standard_1020.elc'), ...
        fullfile(eeglabroot,'sample_locs','standard-10-5-cap385.elp') ...
    };
    for i = 1:numel(candidates)
        if exist(candidates{i},'file')
            EEG = pop_chanedit(EEG, 'lookup', candidates{i});
            return;
        end
    end
    warning('Channel location template not found. Skip.');
catch ME
    warning('Failed to set channel locations: %s', ME.message);
end
end

function EEG = delete_if_present(EEG, chans)
labels = {EEG.chanlocs.labels};
mask = ismember(labels, chans);
if any(mask)
    EEG = pop_select(EEG, 'nochannel', labels(mask));
else
    fprintf('    (None of %s present; skip)\n', strjoin(chans, ', '));
end
end

function EEG = run_ica_and_remove_eyes(EEG, pthr)
EEG = pop_runica(EEG, 'icatype','runica','extended',1,'interrupt','off');
if exist('pop_iclabel','file') == 2
    EEG = pop_iclabel(EEG, 'default');
    ok = isfield(EEG,'etc') && isfield(EEG.etc,'ic_classification') && ...
         isfield(EEG.etc.ic_classification,'ICLabel');
    if ok
        probs = EEG.etc.ic_classification.ICLabel.classifications;
        if size(probs,2) >= 3
            rm = find(probs(:,3) >= pthr);  % Eye class
            if ~isempty(rm)
                fprintf('    Removing %d ocular IC(s) (>= %.2f)\n', numel(rm), pthr);
                EEG = pop_subcomp(EEG, rm, 0);
            else
                fprintf('    No ocular IC >= %.2f found.\n', pthr);
            end
        end
    else
        warning('ICLabel result not found; skip IC-based removal.');
    end
else
    warning('ICLabel plugin not found; skip ocular removal.');
end
end

function EEG = do_reref(EEG, mode)
try
    if ischar(mode) && strcmpi(mode,'average')
        EEG = pop_reref(EEG, []);
    elseif iscell(mode)
        labels = {EEG.chanlocs.labels};
        idx = find(ismember(labels, mode));
        if isempty(idx)
            warning('Re-ref channels %s not found. Using average.', strjoin(mode, ', '));
            EEG = pop_reref(EEG, []);
        else
            EEG = pop_reref(EEG, idx);
        end
    else
        EEG = pop_reref(EEG, []);
    end
catch ME
    warning('Re-reference failed (%s). Using average.', ME.message);
    EEG = pop_reref(EEG, []);
end
end

function [segs, times, info] = extract_segments_varlen(EEG, code_start, code_end, pre_ms, max_gap_ms, require_full_baseline)
% For each END, take the nearest preceding START that is between previous END and current END.
pre_ms = abs(pre_ms);
srate = EEG.srate;
types = arrayfun(@(e) type2num(e.type), EEG.event);
lats  = [EEG.event.latency];

idx_s = find(types == code_start);
idx_e = find(types == code_end);

segs = {}; times = {};
info = struct('i_start',{},'i_end',{},'t14_samp',{},'tend_samp',{},'dur_ms',{},'valid',{});
if isempty(idx_e) || isempty(idx_s), return; end

pre_samp = round(pre_ms/1000 * srate);
[~, order_e] = sort(lats(idx_e));
idx_e = idx_e(order_e);

last_end_lat = -inf;
for ee = 1:numel(idx_e)
    i_end = idx_e(ee);
    t_end = lats(i_end);
    cand_s = idx_s(lats(idx_s) > last_end_lat & lats(idx_s) < t_end);

    valid  = true;
    i_start = NaN; t_start = NaN; dt_ms = NaN;

    if isempty(cand_s)
        valid = false;
    else
        i_start = cand_s(end);
        t_start = lats(i_start);
        dt_ms   = (t_end - t_start)/srate*1000;
        if ~(dt_ms >= 0 && dt_ms <= max_gap_ms), valid = false; end
    end

    if valid
        s1 = round(t_start - pre_samp);
        s2 = round(t_end);
        if s1 < 1
            if require_full_baseline, valid = false; else, s1 = 1; end
        end
        if s2 > size(EEG.data,2), valid = false; end

        if valid
            seg = EEG.data(:, s1:s2);
            t0  = ((s1:s2) - t_start) / srate;  % seconds; 0 at event 14
            segs{end+1}  = seg;                 %#ok<AGROW>
            times{end+1} = t0;                  %#ok<AGROW>
        end
    end

    info(end+1) = struct( ...
        'i_start',  i_start, ...
        'i_end',    i_end, ...
        't14_samp', t_start, ...
        'tend_samp',t_end, ...
        'dur_ms',   dt_ms, ...
        'valid',    logical(valid));

    last_end_lat = t_end;
end
end

function v = type2num(x)
if isnumeric(x), v = x; return; end
if ischar(x) || isstring(x)
    y = str2double(x);
    if ~isnan(y), v = y; return; end
end
v = NaN;
end

function [segs_blc, rej] = baseline_correct_segments(segs, times, base_ms)
n = numel(segs);
segs_blc = segs;
rej = false(n,1);
for i = 1:n
    t = times{i};
    s = segs{i};
    if numel(t) < 2 || isempty(s), rej(i) = true; continue; end
    dt_ms = median(diff(t))*1000;
    tol   = max(1, 0.6*abs(dt_ms));
    idx = find( t*1000 >= (base_ms(1)-tol) & t*1000 <= (base_ms(2)+tol) );
    if numel(idx) < 3, rej(i) = true; continue; end
    base = mean(s(:, idx), 2);
    segs_blc{i} = s - base;
end
end

function T = make_info_level_table(info, keep_mask_info)
n = numel(info);
seg_idx_for_info = nan(n,1);
valid_idx = find([info.valid] == 1);
seg_idx_for_info(valid_idx) = (1:numel(valid_idx))';

i14 = arrayfun(@(x) x.i_start,   info)';
iend= arrayfun(@(x) x.i_end,     info)';
t14 = arrayfun(@(x) x.t14_samp,  info)';
te  = arrayfun(@(x) x.tend_samp, info)';
dur = arrayfun(@(x) x.dur_ms,    info)';

T = table((1:n)', seg_idx_for_info, i14, iend, t14, te, dur, logical([info.valid]'), logical(keep_mask_info(:)), ...
    'VariableNames', {'info_row','seg_idx','ev14_index','evEnd_index','t14_sample','tEnd_sample','dur_ms','valid','keep'});
end
