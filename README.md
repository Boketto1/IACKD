# IACKD
The Intention–Action Conflict EEG–Hand Kinematics Dataset (IACKD) is a joint resource for studying congruent and incongruent intention–action conditions during unimanual control.
# Intention–Action Conflict EEG–Hand Kinematics Dataset (IACKD)

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Platform](https://img.shields.io/badge/Platform-Win%20%7C%20macOS%20%7C%20Linux-blue.svg)](#)
[![Data DOI](https://img.shields.io/badge/DOI-<DOI>-orange.svg)](<DOI_link>)

Synchronized **EEG** and **3-D hand kinematics** collected under **congruent** and **incongruent** intention–action conditions during unimanual control.  
15 participants, 7,040 trials (6,720 valid). Designed for **intention decoding**, **visuomotor conflict** studies, and **continuous trajectory decoding**.

---

## 1) Get the data

- 🌐 **Public repository (with DOI):** <DOI_link>  
- 📦 **Mirror / sample subset (optional):** <sample_zip_link>  
- 📁 **Folder layout preview:** see [`Data Records`](#data-records-summary)

> Data license: **CC BY 4.0**. Please cite (see **Cite this work**).

---

## 2) Repository structure
<pre>
.
├─ Datasets/s3/preprocessed data/# small sample files for quick tests
├─ code/scripts/                 # MATLAB & Python analysis scripts
│  ├─ matlab/
│  │  ├─ preprocess_eeg_pipeline.m
│  │  ├─ preprocess_premotor_55_14.m
│  │  ├─ process_leap_dual_csv.m
│  │  ├─ align_eeg_leap.m
│  │  ├─ calibrate_pairs_with_events.m
│  │  └─ export_eeg_leap_aligned.m
│  └─ python/
│     ├─ rp4pre_move.py        # RP / ERD-ERSsync 
│     ├─ res_completion_time.py  # residuals & completion-time summaries     
│     └─ ersp_itc_trajectory.py  # ERSP / ITC / trajectories
├─ code/output                 # validation figures (optional)
└─ README.md
</pre>



## 3) Quick start

### MATLAB
### MATLAB example for pre_move
```matlab
preprocess_premotor_55_14( ...
  'D:\datasets\raw data\s1\eeg\Acquisition 01l.cdt', ...           % Path to the .cdt EEG file
  'D:\datasets\raw data\s1\leap motion\left_ball_1.csv ', ...       % CSV file containing color_ball / move_direct information
  'D:\output\pre_eeg_L1.mat' ...                                    % Output path for the resulting .mat file
);
```


### MATLAB example for move_execution
```matlab
% Part 1
% EEG segment preprocessing
preprocess_eeg_pipeline % Modify the internal path to the raw EEG file: "D:\datasets\raw data\s1\eeg\Acquisition 01l.cdt"
```
```matlab
% Part 2
% Align the time axes of *_ball_*.csv and *_leap_*.csv files
ball_csv = 'D:\datasets\raw data\s1\leap motion\left_ball_1.csv';
leap_csv = 'D:\datasets\raw data\s1\leap motion\left_leap_1.csv';

M = exact_offset_from_xy(ball_csv, leap_csv, ...
    'tol_mm', 0, ...
    'normalize_time', false, ...
    'by_trial', true, ...       % Recommended to enable this to avoid cross-trial mismatches
    'show_plot', true, ...
    'bin_ms', 5);
% Edit lines 51–52 in process_leap_dual_csv5.m and replace a and b with the values output by exact_offset_from_xy.m
```
```matlab
% Part 3
% Merge two Leap Motion streams and generate a unified file
process_leap_dual_csv( ...
  'D:\datasets\raw data\s1\leap motion\left_ball_1.csv', ...
  'D:\datasets\raw data\s1\leap motion\left_leap_1.csv', ...
  'D:\output\leap\leap_L1.mat', ... % Output path for the preprocessed and aligned Leap Motion data
  100, ...       % Output sampling rate (Hz)
  250, ...       % Threshold for severe frame loss in hand tracking (ms)
  4, ...         % Band-pass filter cutoff for hand trajectory (Hz)
  4);            % 4th-order Butterworth filter
```
```matlab
% Part 4
% Align EEG data with Leap Motion data
align_eeg_leap( ...
  'D:\output\eeg\eeg_L1.mat', ...
  'D:\output\leap\leap_L1.mat', ...
  'D:\output\out', ... % Output folder for align_eeg_leap; generates eeg_leap_pairs.mat and eeg_leap_pairs_summary.csv
  300, ...  % Tolerance for duration difference between EEG and kinematic data (ms)
  2000);    % Penalty for skipping trials in dynamic programming (controls matching count)
```
```matlab
% Part 5
% Recalibrate the aligned pairs and compute residual tolerance
calibrate_pairs_with_events2( ...
  'D:\output\eeg\eeg_L1.mat', ... 
  'D:\output\leap\leap_L1.mat', ...
  'D:\output\out\eeg_leap_pairs.mat', ...
  'D:\output\out', ... % Output directory for this script; generates eeg_leap_pairs_events_from_dp.mat
  600); % Residual tolerance for event alignment (ms)
  ```
```matlab
% Part 6
% Export the final aligned multimodal dataset
export_eeg_leap_aligned( ...
  'D:\output\eeg\eeg_L1.mat', ...
  'D:\output\leap\leap_L1.mat', ...
  'D:\output\out\eeg_leap_pairs.mat', ...
  'D:\output\out\eeg_leap_pairs_events_from_dp.mat', ...
  'D:\output\final_out\eeg_leap_L1.mat');
  ```
### python
```python
import scipy.io as sio
m = sio.loadmat('eeg_leap_L1.mat', squeeze_me=True, struct_as_record=False)
out = m['OUT'][0]
print(out['EEG'].shape)      # (nChan, T) at 100 Hz
print(out['x_mm'].shape)     # (T,) 3-D trajectory also at 100 Hz
```
## 4) Environment / dependencies
MATLAB R2023b (tested) + EEGLAB v2025

Optional toolboxes: Signal Processing

Python ≥ 3.9, packages: numpy, mat73, scipy, matplotlib, mne (optional)

Set a working directory with write permission for exported .mat/.csv files.

## 5) Data Records summary
Raw EEG (CURRY): .cdt + sidecars .ceo, .dpa (1024 Hz; events 14/66)

Leap Motion CSV: *_ball_*.csv (labels & screen x), *_leap_*.csv (x_mm,y_mm,z_mm @170 Hz)

Pre-move EEG preprocessed data/pre_move/

segs_1s, times_1s (fixed −1.0→0.0 s)

segs_var, times_var (full 55→14, variable length)

Aligned EEG–hand preprocessed data/move_execution/

OUT struct per trial: EEG, x_mm/y_mm/z_mm, ball_x_pix, t_ms, labels, (a,b), start_res_ms, cov

Scripts preprocessed data/script/ + README.md (how to run)

Units: EEG (µV), position (mm), screen (px), time (ms/s). Missing coverage → NaN.

## 6) Reproducibility & notes
Event alignment uses affine mapping t_EEG → t_Leap with parameters (a,b).

Coverage ratio cov is computed over the valid Leap interval; out-of-coverage samples are NaN.

Trials with |start_res_ms| > 30 ms can be excluded (recommended).

Use segs_1s/times_1s for equal-length models; use segs_var/times_var to keep the full pre-move context.
## 7) License
Data: Creative Commons CC BY 4.0

Code: MIT (unless otherwise stated in file headers)
## 8) Contact & issues
GitHub Issues: please open a ticket with steps to reproduce and logs.
## 9) Acknowledgements
Thanks to all participants and lab members who supported data collection and validation.
