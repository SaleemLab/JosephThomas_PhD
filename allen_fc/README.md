# Allen Visual Coding Neuropixels, functional-connectivity sessions: per-session .mat files

One file per session, `session_<id>.mat`. Each holds every sorted unit (noise clusters removed) on every probe, with its spike times for the whole recording. Alongside are the stimulus tables, optotagging pulses, running and eye tracking. Everything is read straight from the Allen NWB file, and nothing is filtered, binned or re-zeroed. LFP is not included.

**Where things are**

The data are on the server, and everything else is in this repository folder.

| | |
|---|---|
| `Y:\ibn-vision\DATA\ExternalDatasets\ABI_VisualCoding_Neuropixels\exports_for_joseph\` | the 26 session files, `session_<id>.mat` (about 14 GB in total). Nothing else is there |
| `psth_tutorial.m` | a step-by-step tutorial: start here |
| `example_v1_speed_state.m` | V1 population PSTH by dot speed, running vs stationary, for every session, with per-session and across-session figures |
| `helpers/` | the small functions both scripts use |
| `figures/` | figures saved by the scripts |
| `export_summary.csv` | one row per session: unit counts, V1 probe, eye tracking, invalid spans, dot-motion protocol, file size |

Open `psth_tutorial.m` and run it one section at a time. Both scripts add `helpers/` to the MATLAB path themselves. Their `data_dir` setting points at the server folder; change it if the drive is mapped to a different letter on your machine.

```matlab
data_dir = 'Y:\ibn-vision\DATA\ExternalDatasets\ABI_VisualCoding_Neuropixels\exports_for_joseph';
S = load(fullfile(data_dir, 'session_779839471.mat'));
units = struct2table(S.units);       % one row per unit
st = S.spike_times{5};               % spike times (s) of units row 5
fl = struct2table(S.stim.flashes);   % one row per flash
```

**Conventions**

- All times are in seconds on the Allen master clock, the same clock for spikes, stimuli, running and eye tracking. Nothing needs aligning.
- `S.spike_times{i}` belongs to row `i` of `S.units`, and so does `S.waveform(i, :)`.
- Units are sorted by probe, then from the probe's top (pia) to its tip.
- Fields named `*local_index` are Allen's 0-based indices. Everything else follows MATLAB.
- Text columns are cell arrays of char. `struct2table` turns any struct of columns here into a table.

**Size (pilot session 779839471, one of the larger ones):** 2,235 units, 908 of them QC-pass, 87 M spikes. 550 MB on disk and 0.8 GB in memory. `load` takes about 8 s in MATLAB R2025a.

## Variables

### `units`: one row per unit

| Field | Meaning |
|---|---|
| `unit_id` | Allen unit id (unique across the whole dataset) |
| `probe`, `probe_id` | probe name (`probeA`–`probeF`) and Allen probe id |
| `area` | CCF structure acronym of the unit's peak channel, e.g. `VISp`, `LGd`, `CA1`. `VIS` = visual cortex not assigned to an area; `''` = no label (outside the brain) |
| `qc_pass` | `true` if `isi_violations < 0.5`, `amplitude_cutoff < 0.1` and `presence_ratio > 0.9` (`meta.qc`). The usual starting filter |
| `is_rs` | regular-spiking (`waveform_duration > 0.4` ms); `false` = fast-spiking |
| `n_spikes` | number of spikes in `spike_times{i}` |
| `firing_rate` | mean rate over the session (Hz, Allen) |
| `peak_channel_id`, `peak_channel_local_index`, `peak_channel_valid` | channel where the unit is largest: Allen channel id, 0-based index on its probe (0–383), and `false` if Allen marked that channel invalid |
| `probe_vertical_position_um`, `probe_horizontal_position_um` | position of the peak channel on the shank (µm; vertical is measured from the tip) |
| `ccf_ap_um`, `ccf_dv_um`, `ccf_lr_um` | Allen CCF coordinates of the peak channel (µm). NaN for the few invalid channels missing from Allen's channel table |
| `v1_depth_from_pia_um`, `v1_depth_norm`, `in_v1_column` | V1 probe only (`meta.v1_probe`): depth along the cortical normal from the pial surface (µm), the same as a fraction of cortical thickness (0 = pia, 1 = white matter), and `0 ≤ v1_depth_norm ≤ 1`. NaN on other probes, and in sessions without our column geometry. Estimated from the CCF registration and the LFP surface, so treat it as approximate |
| `isi_violations`, `amplitude_cutoff`, `presence_ratio`, `snr`, `isolation_distance`, `l_ratio`, `d_prime`, `nn_hit_rate`, `nn_miss_rate`, `silhouette_score`, `max_drift`, `cumulative_drift` | Allen's spike-sorting quality metrics (definitions in the AllenSDK documentation, "Visual Coding Neuropixels" quality metrics) |
| `amplitude`, `waveform_duration`, `waveform_halfwidth`, `PT_ratio`, `repolarization_slope`, `recovery_slope`, `spread`, `velocity_above`, `velocity_below` | Allen's waveform metrics (amplitude in µV, durations in ms, spread in µm) |
| `cluster_id`, `local_index` | Kilosort cluster id, and the unit's 0-based index within its probe's sort |

### `spike_times`: n_units × 1 cell

Each cell is a column vector of spike times (s). Spike times cover the whole recording, stimuli and optotagging included.

### `waveform`: n_units × 82

Mean waveform on the peak channel (µV), sampled at 30 kHz (2.7 ms). Allen chooses the peak channel from the spike-sorting template. For about 5% of units the largest mean waveform sits a few channels away from it.

### `stim.<name>`: one struct of columns per stimulus, one row per presentation

Every table has `start_time`, `stop_time`, `stimulus_name`, and (apart from `spontaneous`) `stimulus_block`. The parameter columns are numeric where every value is a single number. Vector-valued settings such as `size = '[250.0, 250.0]'` stay as text.

| Table | What it is | Main condition columns |
|---|---|---|
| `spontaneous` | grey screen between stimuli. The longest interval (~30 min) is `meta.spont_block_s` | (none) |
| `flashes` | 250 ms full-field flashes | `color` (1 = white, −1 = black) |
| `drifting_gratings_contrast` | full-field drifting gratings at many contrasts | `orientation`, `contrast` |
| `drifting_gratings_75_repeats` | 4 orientations × 2 contrasts, 75 repeats each, 2 s on | `orientation`, `contrast` |
| `gabors` | small Gabor patches for receptive-field mapping | `x_position`, `y_position`, `orientation` |
| `natural_movie_one_more_repeats`, `natural_movie_one_shuffled` | 30 s natural movie, intact and with frames shuffled. **One row per frame** (900 frames per repeat) | `frame` (0–899) |
| `dot_motion` | random-dot motion | `Dir` (direction, deg), `Speed` (the stimulus software's value), `Speed_deg_s` (the speed shown, deg/s), `coherence` |
| `dot_motion_nonstandard` | only in 766640955 and 767871931: the dots used a different protocol (see below) | `Dir`, `Speed` |

**Dot speeds.** In 24 of the 26 sessions the dots move in 4 directions (`Dir` −45, 0, 45, 90) at 7 speeds. Their `Speed` values map to deg/s as follows, and `Speed_deg_s` holds the deg/s value. `meta.dot_speed_deg_s` holds the same table.

| `Speed` | 0 | 0.0005 | 0.001 | 0.005 | 0.01 | 0.02 | 0.04 |
|---|---|---|---|---|---|---|---|
| deg/s | 0 | 16 | 32 | 64 | 128 | 256 | 512 |

The raw values are not proportional to deg/s, so the two are matched by rank. The first two FC sessions, 766640955 and 767871931, instead showed 8 directions (0–315) at 4 speeds (0.001–0.02), plus 15 trials with neither (NaN). Their table is renamed `dot_motion_nonstandard`, so code that asks for `dot_motion` skips them rather than pooling them by mistake, and `meta.dot_motion_note` says why. Their deg/s speeds are not known, so they have no `Speed_deg_s`.

### `epochs`: one row per stimulus block

`stimulus_name`, `stimulus_block`, `start_time` (first presentation), `stop_time` (last presentation), `n_presentations`. Each spontaneous interval is its own row, with `stimulus_block` NaN. This is the quickest way to see the session's timeline.

### `opto`: optotagging light pulses, at the end of the session

`start_time`, `stop_time`, `duration` (s), `level` (Allen's light level), `stimulus_name` (`pulse`, `fast_pulses`, `raised_cosine`), `condition` (description). These only tag cells in the Cre lines (Sst, Pvalb, Vip; see `meta.genotype`).

### `invalid`: time spans Allen flagged as bad

`start_time`, `stop_time`, `probe_id` (−1 = all probes). Spikes on that probe within these spans should be ignored. Most sessions have none.

### `running`, `eye`

- `running.t`, `running.speed_cm_s`: wheel speed (sampled irregularly, ~30–50 Hz).
- `eye.t`, `eye.pupil_area`, `eye.eye_area`, `eye.gaze_x_cm`, `eye.gaze_y_cm`, `eye.gaze_azimuth_deg`, `eye.gaze_altitude_deg`: Allen's filtered eye tracking (30 Hz). Gaze is the position on the screen (cm) and the same position as an angle (deg). Values are NaN during blinks and bad fits (17% of frames in the pilot).
- Pupil and eye area are in **arbitrary units**. The NWB labels them pixels², but the values are around 1e-3. Compare them within a session only; `sqrt(pupil_area / pi)` is proportional to diameter.
- Two sessions have no eye tracking, and their `eye` is an empty struct (`meta.has_eye_tracking` is `false`).

### `meta`

`session_id`, `session_type`, `subject_id`, `genotype`, `sex`, `age_in_days`, `species`, `session_start_time`, `spont_block_s` (the longest spontaneous interval, [start stop]), `v1_probe`, `has_eye_tracking`, `n_units`, `n_qc_pass`, `qc` (thresholds behind `units.qc_pass`), `rs_fs_threshold_ms`, `dot_motion_note` (whether this session's dot motion is the standard protocol), `dot_speed_deg_s` (`Speed` → deg/s), `probes` (name, id, unit count and sampling rates of each probe), `source`, `exported`, `export_code`.

## PSTHs step by step: `psth_tutorial.m`

**Start here.** `psth_tutorial.m` is a script to run one section at a time. Each section is a short snippet you can copy into your own code. The settings at the top are the session file, stimulus, window, bin width, smoothing, single-cell normalisation, the state window, a minimum trial count and the areas. The default areas are the six visual cortical areas (`VISp`, `VISl`, `VISal`, `VISrl`, `VISam`, `VISpm`) and the two visual thalamic nuclei (`LGd`, `LP`).

Not every session has cells in every area. Section 1 prints the unit count of each chosen area. Section 4 stops with a message if its area has none. Section 11 keeps an empty area's field, but with a NaN MUA (not zeros) and a `SingleCell` of size M × T × 0. Section 12 skips the session. A session without the chosen stimulus stops in section 2 with a message.

The sections, in order:

1. Load a session and list its stimuli and units per area.
2. List the trials of one stimulus and its conditions.
3. Label each trial running, stationary or neither, and count the trials of each condition in each state.
4. Bin one area's spikes into trials × bins × units, then smooth.
5. Get the PSTH (mean ± SEM over trials) of one unit and of the area's MUA.
   - 5b. Normalise each cell (sqrt, z-score, baseline-subtracted, soft-norm), then average the normalised units into a population PSTH in which every unit counts equally.
6. Split by stimulus condition.
7. Split by one parameter (e.g. orientation), pooling over the rest.
8. Split by behavioural state. A PSTH from fewer than `min_trials` trials stays NaN, and states can be pooled over conditions.
   - 8b. Split by state and by specific stimulus parameters together. You get a grid with one panel per value of `row_param` × `col_param` (e.g. dot `Dir` × `Speed`), running and stationary overlaid, and the trial counts in each title. Set `row_param = ''` to pool over it. A tuning curve follows: the response while the stimulus is on, against `col_param`, for each state.
9. Split by cell type: RS vs FS, V1 depth, and a rough optotag test (Cre lines only).
10. Make cross-validation folds within each condition, then compute split-half reliability per unit. Odd/even and first/second-half subsets are also shown.
11. Build the agreed layout: `Session.<stimulus>.TrialsStim`, `.TrialID`, `.Area.<name>.MUA` / `.SingleCell`.
12. Compute the same PSTH for every session in the folder.

It uses six small functions in `helpers/`, each short enough to read in a minute:

| Function | What it does |
|---|---|
| `[trials, conditions] = trial_table(S, stim_name)` | One row per trial: `onset`, `duration`, every parameter that varies, `condition`, `block`. Natural movies give one trial per 30 s repeat |
| `[state, speed] = trial_state(S, onsets, window, opt)` | 1 = running, 2 = stationary, 0 = neither. Running means speed > 3 cm/s for ≥ 80% of the window; stationary means \|speed\| < 1 cm/s for ≥ 80% (all adjustable). The speed is smoothed over 1 s first. Returns mean speed and both fractions for your own rule |
| `counts = bin_spikes(S.spike_times(rows), onsets, edges)` | spike counts, trials × bins × units, with `edges` in s from onset |
| `y = smooth_time(x, bin_s, method, width_s)` | smoothing along time (dimension 2): `'none'`, `'boxcar'`, `'gaussian'` (SD `width_s`), or `'causal'` (the past half of the Gaussian, so responses don't start early). The kernel is rescaled at the window edges |
| `y = normalise_cells(x, method, opt)` | Normalises each unit of trials × bins × units on its own. `'none'`; `'sqrt'`; `'zscore'`, with mean and SD from `opt.ref_bins` (default all bins, or e.g. `t < 0` for the baseline); `'baseline'` (subtracts the `ref_bins` mean); `'softnorm'`, which is x / (range + `opt.c`, 5 spikes/s) with the range taken over the condition-averaged PSTH when `opt.condition` is given. In the agreed layout the MUA is summed before normalising, and `SingleCell` is normalised with `cell_norm` |
| `B = trial_block(S, stim_name, areas, opt)` | All of the above in one call. Returns one stimulus of one session in the agreed format: `B.TrialsStim`, `B.Conditions`, `B.TrialID`, `B.Area.<name>.MUA` / `.SingleCell` / `.units`, `B.Trials`, `B.t`. Soft-norm by default. It also sets trials that overlap an Allen invalid span to NaN, for the units on that probe. An area with no units gets a NaN MUA and a `SingleCell` of size M × T × 0 |

On the pilot session, the whole tutorial takes about 70 s for `dot_motion` and 3.5 minutes for the 75-repeat gratings in MATLAB R2025a, most of it in section 11 (all eight areas). Binning one area (V1, 130 QC units, 600 trials × 300 bins) takes a few seconds and needs 94 MB. The pilot mouse hardly ran during the gratings, with 0–3 running trials per condition. Section 8 shows how `min_trials` handles that.

## Notes

- **Source:** Allen Brain Observatory, Visual Coding Neuropixels (Siegle et al. 2021, *Nature* 592:86–92), functional-connectivity sessions. Please cite it, and see the Allen terms of use.
- **Excluded:** units Allen labels `noise` are left out. Allen's own default filter is stricter on presence ratio (> 0.95) than `qc_pass`.
- **Spikes in invalid spans** are still in `spike_times`. Use `invalid` to mask them.
- **Missing data:** LFP, raw wheel rotation, per-spike amplitudes and full multichannel waveforms are not in these files. They are in the NWB files and can be exported if needed.
