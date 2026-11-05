%% PSTHs from the Allen FC session files: a step-by-step tutorial
% Run one section at a time (Ctrl+Enter in the editor). Each section is a short snippet that only needs
% the variables made by the sections above it, so copy whichever ones you need into your own code.
%
% The session files are on the server (data_dir below). It uses these small functions, in helpers/:
%   trial_table      one row per trial: onset, stimulus parameters, condition number
%   trial_state      running / stationary / neither on each trial
%   bin_spikes       spike counts, trials x time bins x units
%   smooth_time      smoothing along time
%   normalise_cells  z-score, soft-norm, sqrt, ... of each unit
%   trial_block      all of the above in one call, returning the agreed format (section 11)
% Every one of them is short: open it to see exactly what it does.
%
% For the V1 population PSTH by dot speed, running vs stationary, across every session, run
% example_v1_speed_state.m instead: it is the finished analysis, built on trial_block.

%% Settings: change these

here = fileparts(mfilename('fullpath'));   % this script's folder ('' = MATLAB's current folder)
addpath(fullfile(here, 'helpers'));       % the small functions this uses
data_dir = 'Y:\ibn-vision\DATA\ExternalDatasets\ABI_VisualCoding_Neuropixels\exports_for_joseph';   % the session_<id>.mat files (on the server): change if your drive letter differs
file = 'session_779839471.mat';
stim_name = 'dot_motion';  % any field of S.stim (section 1 lists them)
window = [-0.2 1.8];         % s relative to stimulus onset
bin_s = 0.01;                % bin width (s); the window should hold a whole number of bins
smooth_method = 'gaussian';  % 'none' | 'boxcar' | 'gaussian' | 'causal'
smooth_s = 0.035;             % Gaussian SD, or boxcar width (s)
cell_norm = 'softnorm';      % 'none' | 'sqrt' | 'zscore' | 'baseline' | 'softnorm' (section 5b)
norm_ref = 'all';            % z-score / baseline reference: 'all' bins, or 'pre' (bins before onset)
softnorm_c = 5;              % soft-norm constant (spikes/s)
state_window = [-0.2 1.8];        % s relative to onset over which running / stationary is judged
min_trials = 5;              % a PSTH from fewer trials than this is left as NaN
qc_only = true;              % only units that pass Allen's quality thresholds (units.qc_pass)
areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', ...   % the 6 visual cortical areas
         'LGd', 'LP'};                                              % and the 2 visual thalamic nuclei
% Not every session has cells in every area. Section 1 prints how many each one has, an area with none is
% skipped, and in section 11 its MUA is NaN rather than 0. unique(units.area) lists every area in a session.

%% 1. Load one session

S = load(fullfile(data_dir, file));
units = struct2table(S.units);               % one row per unit; row i goes with S.spike_times{i}
fprintf('Session %d, %s, %d units (%d pass QC)\n', S.meta.session_id, S.meta.genotype, ...
    height(units), sum(units.qc_pass));
disp(fieldnames(S.stim))                     % the stimuli in this session
disp(groupcounts(units(units.qc_pass, :), 'area'))   % units per area, every area in the session

% the chosen areas, including those with no cells in this session
n_units = zeros(numel(areas), 1);
for a = 1:numel(areas)
    n_units(a) = sum(strcmp(units.area, areas{a}) & (units.qc_pass | ~qc_only));
end
disp(table(areas', n_units, 'VariableNames', {'area', 'n_units'}))

%% 2. Trials of one stimulus, and its conditions

% dot_motion: two sessions (766640955, 767871931) showed different dot speeds and directions. Their table is
% S.stim.dot_motion_nonstandard instead, so they never get pooled with the rest; S.meta.dot_motion_note says why.
% In the other 24, trials.Speed_deg_s is the speed in deg/s (0, 16, 32, 64, 128, 256, 512).
if ~isfield(S.stim, stim_name)
    error('Session %d has no %s: pick one from section 1. %s', S.meta.session_id, stim_name, S.meta.dot_motion_note);
end
[trials, conditions] = trial_table(S, stim_name);
disp(conditions)                             % one row per parameter combination, with n_trials
fprintf('First trials: \n')
disp(head(trials))                           % one row per trial

%% 3. Behavioural state on each trial: 1 = running, 2 = stationary, 0 = neither

[trials.state, speed] = trial_state(S, trials.onset, state_window);
trials = [trials, speed];

% How many trials of each condition fall in each state. Some mice hardly run, so expect zeros.
conditions.n_running = accumarray(trials.condition, trials.state == 1);
conditions.n_stationary = accumarray(trials.condition, trials.state == 2);
disp(conditions)

% Other thresholds:  trial_state(S, trials.onset, state_window, struct('run_cm_s', 5, 'run_frac', 0.5))
% Or your own rule:  trials.state = 2 - (trials.mean_speed_cm_s > 3);

%% 4. Bin the spikes of one area

area = 'VISp';
rows = find(strcmp(units.area, area) & (units.qc_pass | ~qc_only));   % rows of units / S.spike_times
if isempty(rows)
    error('No %s units in session %d: pick another area from section 1.', area, S.meta.session_id);
end
edges = window(1):bin_s:window(2);           % bin edges, s from onset
t = edges(1:end - 1) + bin_s / 2;            % bin centres, for plotting

counts = bin_spikes(S.spike_times(rows), trials.onset, edges);     % trials x bins x units

% Allen flagged a few time spans as bad (3 sessions). Set the trials they overlap to NaN, for the units on
% that probe (probe_id -1 = every probe). Most sessions have none, so this loop usually does nothing.
for j = 1:numel(S.invalid.start_time)
    bad = trials.onset + window(2) > S.invalid.start_time(j) & trials.onset + window(1) < S.invalid.stop_time(j);
    on_probe = S.invalid.probe_id(j) == -1 | units.probe_id(rows) == S.invalid.probe_id(j);
    counts(bad, :, on_probe) = NaN;
end

rate = smooth_time(counts / bin_s, bin_s, smooth_method, smooth_s);  % spikes/s
fprintf('%s: %d trials x %d bins x %d units\n', area, size(rate));

%% 5. Mean and SEM over trials, for one unit and for the whole area (MUA)

% rate(trials, bins, units): average over dimension 1 to get a PSTH.
unit = 1;                                    % one unit: a column of units(rows(unit), :)
psth_unit = mean(rate(:, :, unit), 1, 'omitnan');   % omitnan: skips trials set to NaN above

MUA = sum(rate, 3);                          % trials x bins: summed over the area's units
psth = mean(MUA, 1, 'omitnan');
sem = std(MUA, 0, 1, 'omitnan') / sqrt(sum(~isnan(MUA(:, 1))));

figure;
subplot(1, 2, 1); plot(t, psth_unit, 'k');
xlabel('s from onset'); ylabel('spikes/s'); title(sprintf('unit %d', units.unit_id(rows(unit)))); box off;
subplot(1, 2, 2); hold on;
fill([t, fliplr(t)], [psth - sem, fliplr(psth + sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t, psth, 'k');
xlabel('s from onset'); ylabel('spikes/s'); title(sprintf('%s MUA, all trials, mean +- SEM', area)); box off;

%% 5b. Normalise each cell

% normalise_cells works on trials x bins x units and treats each unit on its own.
if strcmp(norm_ref, 'pre')
    ref_bins = t < 0;                        % mean and SD from the bins before stimulus onset
else
    ref_bins = true(size(t));                % from every bin of every trial
end
rate_z = normalise_cells(rate, 'zscore', struct('ref_bins', ref_bins));
rate_sqrt = normalise_cells(rate, 'sqrt');
% sqrt after smoothing makes the first and last bins dip a little: they average fewer bins, so they are
% noisier, and the sqrt of a noisier rate is lower on average. To avoid it, sqrt the counts first:
%   rate_sqrt = smooth_time(sqrt(counts / bin_s), bin_s, smooth_method, smooth_s);
rate_base = normalise_cells(rate, 'baseline', struct('ref_bins', t < 0));   % change from pre-stimulus, spikes/s
rate_soft = normalise_cells(rate, 'softnorm', struct('c', softnorm_c, 'condition', trials.condition));
rate_norm = normalise_cells(rate, cell_norm, struct('ref_bins', ref_bins, 'c', softnorm_c, ...
    'condition', trials.condition));         % the one chosen in Settings

% Population PSTH in which every unit counts equally: average the normalised units.
% (The MUA, sum(rate, 3), is dominated by the units that fire most.)
figure;
versions = {rate, 'none (spikes/s)'; rate_sqrt, 'sqrt'; rate_z, 'z-score'; rate_soft, 'soft-norm'};
for i = 1:4
    subplot(1, 4, i);
    plot(t, mean(mean(versions{i, 1}, 3, 'omitnan'), 1, 'omitnan'), 'k');
    xlabel('s from onset'); title(versions{i, 2}); box off;
end

%% 6. Split by stimulus condition
K = height(conditions);
psth_cond = nan(K, numel(t));                % condition x bins
for k = 1:K
    psth_cond(k, :) = mean(MUA(trials.condition == k, :), 1, 'omitnan');
end

figure; plot(t, psth_cond');
legend(compose('condition %d', 1:K), 'Location', 'eastoutside');   % the conditions table says what each is
xlabel('s from onset'); ylabel('spikes/s'); title(sprintf('%s MUA by condition', area)); box off;

%% 7. Split by one parameter, pooling over the others

% Any column of trials works as a selector, e.g. orientation (gratings), Dir or Speed (dot_motion), color (flashes).
param = 'Speed_deg_s';                       % dot_motion: Speed in deg/s (Speed is the raw stimulus value)
values = unique(trials.(param));
figure; hold on;
for i = 1:numel(values)
    sel = trials.(param) == values(i);
    plot(t, mean(MUA(sel, :), 1, 'omitnan'), 'DisplayName', sprintf('%s %g (%d trials)', param, values(i), sum(sel)));
end
legend('Location', 'eastoutside'); xlabel('s from onset'); ylabel('spikes/s'); box off;

% Two parameters at once:  sel = trials.orientation == 90 & trials.contrast == 0.8;

%% 8. Split by behavioural state, with a minimum number of trials

% Each condition and state gets a PSTH only if it has at least min_trials trials; otherwise it stays NaN.
psth_state = nan(K, numel(t), 2);            % condition x bins x state (1 = running, 2 = stationary)
n_state = zeros(K, 2);
for k = 1:K
    for s = 1:2
        sel = trials.condition == k & trials.state == s;
        n_state(k, s) = sum(sel);
        if n_state(k, s) >= min_trials
            psth_state(k, :, s) = mean(MUA(sel, :), 1, 'omitnan');
        end
    end
end
disp(array2table(n_state, 'VariableNames', {'n_running', 'n_stationary'}))

% Pooled over conditions, which helps when each condition has few running trials
figure; hold on;
labels = {'running', 'stationary'};
for s = 1:2
    sel = trials.state == s;
    if sum(sel) >= min_trials
        plot(t, mean(MUA(sel, :), 1, 'omitnan'), 'DisplayName', sprintf('%s (%d trials)', labels{s}, sum(sel)));
    end
end
legend('Location', 'eastoutside'); xlabel('s from onset'); ylabel('spikes/s');
title(sprintf('%s MUA, all conditions', area)); box off;

%% 8b. Split by behavioural state AND specific stimulus parameters (e.g. dot direction x speed)

% One panel per combination of row_param and col_param values, running and stationary overlaid.
% Set row_param = '' to pool over it (e.g. every direction together, one panel per speed).
col_param = 'Speed_deg_s';                   % any column of trials: 'Speed_deg_s', 'Dir', 'orientation', 'contrast', ...
row_param = 'Dir';
labels = {'running', 'stationary'};
colours = {[0.85 0.33 0.10], [0 0.45 0.74]};

col_values = unique(trials.(col_param));
if isempty(row_param)
    row_values = NaN;
else
    row_values = unique(trials.(row_param));
end
nr = numel(row_values);
nc = numel(col_values);
psth_vs = nan(nr, nc, numel(t), 2);          % row value x col value x bins x state
n_vs = zeros(nr, nc, 2);                     % trials behind each PSTH

figure;
for r = 1:nr
    for c = 1:nc
        is_stim = trials.(col_param) == col_values(c);
        if ~isempty(row_param)
            is_stim = is_stim & trials.(row_param) == row_values(r);
        end
        subplot(nr, nc, (r - 1) * nc + c); hold on;
        for s = 1:2
            sel = is_stim & trials.state == s;
            n_vs(r, c, s) = sum(sel);
            if n_vs(r, c, s) >= min_trials
                psth_vs(r, c, :, s) = mean(MUA(sel, :), 1, 'omitnan');
                plot(t, squeeze(psth_vs(r, c, :, s)), 'Color', colours{s});
            end
        end
        title({sprintf('%s %g', col_param, col_values(c)), ...
            sprintf('run %d, still %d', n_vs(r, c, 1), n_vs(r, c, 2))}, 'FontWeight', 'normal', 'FontSize', 7);
        if c == 1 && ~isempty(row_param)
            ylabel(sprintf('%s %g', row_param, row_values(r)));
        end
        box off;
    end
end
sgtitle(sprintf('%s MUA: running (orange) vs stationary (blue), panels with < %d trials left empty', ...
    area, min_trials));

% Tuning curve: mean response while the stimulus is on, against col_param, for each state (pooled over row_param)
on = t >= 0 & t < median(trials.duration);
response = mean(MUA(:, on), 2);              % one number per trial
figure; hold on;
for s = 1:2
    m = nan(nc, 1);
    e = nan(nc, 1);
    for c = 1:nc
        sel = trials.(col_param) == col_values(c) & trials.state == s;
        if sum(sel) >= min_trials
            m(c) = mean(response(sel), 'omitnan');
            e(c) = std(response(sel), 'omitnan') / sqrt(sum(~isnan(response(sel))));
        end
    end
    errorbar(1:nc, m, e, 'o-', 'Color', colours{s}, 'DisplayName', labels{s});
end
xticks(1:nc); xticklabels(string(col_values)); xlabel(col_param);
ylabel(sprintf('%s MUA while the stimulus is on (spikes/s)', area)); legend; box off;

%% 9. Split by cell type

% Regular-spiking (mostly excitatory) vs fast-spiking (mostly PV interneurons), from spike width
is_rs = units.is_rs(rows);
figure; hold on;
plot(t, mean(mean(rate(:, :, is_rs), 3), 1, 'omitnan'), 'DisplayName', sprintf('RS, %d units', sum(is_rs)));
plot(t, mean(mean(rate(:, :, ~is_rs), 3), 1, 'omitnan'), 'DisplayName', sprintf('FS, %d units', sum(~is_rs)));
legend; xlabel('s from onset'); ylabel('spikes/s per unit'); title(area); box off;

% Cortical depth on the V1 probe (0 = pia, 1 = white matter; NaN elsewhere), e.g. superficial vs deep
depth = units.v1_depth_norm(rows);
superficial = depth < 0.4;
deep = depth >= 0.6;

% Any other column of units can split cells the same way, e.g. units.firing_rate(rows) > 5

%% 9b. Optotagged cells (Sst, Pvalb and Vip Cre lines only)

% S.opto lists the light pulses at the end of the session. A rough tag: a unit that fires much more in the
% 10 ms of light than in the 10 ms before it. Check the result by eye before relying on it.
opto = struct2table(S.opto);
pulses = opto.start_time(strcmp(opto.stimulus_name, 'pulse') & opto.duration == 0.01);
c = bin_spikes(S.spike_times(rows), pulses, [-0.01 0 0.01]);   % pulses x [before, during] x units
before = squeeze(mean(c(:, 1, :), 1));       % mean spikes per pulse, 10 ms before the light
during = squeeze(mean(c(:, 2, :), 1));       % and during it
tagged = during > 0.5 & during > 5 * before;
fprintf('%s: %d of %d units look light-driven (%s)\n', area, sum(tagged), numel(tagged), S.meta.genotype);

%% 10. Trial subsets for cross-validation

% Each condition's trials are shuffled and dealt into folds, so every fold has every condition.
n_folds = 2;
rng(0);                                      % the same folds every time you run this
trials.fold = zeros(height(trials), 1);
for k = 1:K
    idx = find(trials.condition == k);
    idx = idx(randperm(numel(idx)));
    trials.fold(idx) = mod(0:numel(idx) - 1, n_folds) + 1;
end

% PSTHs from fold 1 and fold 2: condition x bins x units
psth_fold = nan(K, numel(t), numel(rows), n_folds);
for f = 1:n_folds
    for k = 1:K
        sel = trials.condition == k & trials.fold == f;
        psth_fold(k, :, :, f) = mean(rate(sel, :, :), 1, 'omitnan');
    end
end

% Split-half reliability of each unit: correlation between its two halves, all conditions in a row
a = reshape(psth_fold(:, :, :, 1), [], numel(rows));
b = reshape(psth_fold(:, :, :, 2), [], numel(rows));
reliability = nan(numel(rows), 1);
for u = 1:numel(rows)
    r = corrcoef(a(:, u), b(:, u));
    reliability(u) = r(1, 2);
end
figure; histogram(reliability, -0.2:0.05:1);
xlabel('split-half correlation of the PSTH'); ylabel('units'); title(area); box off;

% Simpler subsets: odd and even trials, or the first and second half of the session
odd = mod(1:height(trials), 2)' == 1;
first_half = trials.onset < median(trials.onset);

% Folds within a state: run the loop above on trials.condition == k & trials.state == 2

%% 11. The agreed layout: Session.<stimulus>.TrialsStim / TrialID / Area.<name>.MUA / SingleCell
% trial_block.m does sections 2-5b in one call, for any list of areas, and returns the agreed format.
% example_v1_speed_state.m uses it on every session.
opt = struct('window', window, 'bin_s', bin_s, 'smooth_method', smooth_method, 'smooth_s', smooth_s, ...
    'cell_norm', cell_norm, 'softnorm_c', softnorm_c, 'norm_ref', norm_ref, ...
    'state_window', state_window, 'qc_only', qc_only);
B = trial_block(S, stim_name, areas, opt);
%   B.TrialsStim               M x 1, a row of B.Conditions
%   B.TrialID                  M x 1, 1 = running, 2 = stationary, 0 = neither
%   B.Area.VISp.MUA            M x T, summed spikes/s
%   B.Area.VISp.SingleCell     M x T x N, normalised with cell_norm
%   B.Area.VISp.units          which unit each page of SingleCell is
Session = struct('session_id', S.meta.session_id, 'genotype', S.meta.genotype);
Session.(stim_name) = B;

% Selecting trials uses only these fields, e.g. running trials of one stimulus parameter:
%   speed = B.Conditions.Speed_deg_s(B.TrialsStim);
%   pop = mean(mean(B.Area.VISp.SingleCell(speed == 128 & B.TrialID == 1, :, :), 3), 1, 'omitnan');

% For example, the running / stationary ratio of the V1 MUA:
%   ratio = mean(B.Area.VISp.MUA(B.TrialID == 1, :), 'omitnan') ./ mean(B.Area.VISp.MUA(B.TrialID == 2, :), 'omitnan');
% save(sprintf('session_%d_%s.mat', Session.session_id, stim_name), '-struct', 'Session', '-v7.3');



%% 12. The same PSTH for every session in a folder

files = dir(fullfile(data_dir, 'session_*.mat'));
files = files(~cellfun(@isempty, regexp({files.name}, '^session_\d+\.mat$')));   % not saved outputs
psth_all = {};
for i = 1:numel(files)
    Si = load(fullfile(files(i).folder, files(i).name), 'units', 'spike_times', 'stim', 'meta');
    if ~isfield(Si.stim, stim_name)          % e.g. the two sessions with nonstandard dot motion
        fprintf('%s: no %s, skipped\n', files(i).name, stim_name);
        continue
    end
    Ui = struct2table(Si.units);
    r = find(strcmp(Ui.area, area) & Ui.qc_pass);
    if isempty(r)
        fprintf('%s: no %s units, skipped\n', files(i).name, area);
        continue
    end
    tr = trial_table(Si, stim_name);
    X = smooth_time(bin_spikes(Si.spike_times(r), tr.onset, edges) / bin_s, bin_s, smooth_method, smooth_s);
    psth_all{end + 1} = mean(sum(X, 3), 1);  %#ok<SAGROW>
    fprintf('%s: %d %s units\n', files(i).name, numel(r), area);
end
