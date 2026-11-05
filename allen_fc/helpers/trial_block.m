function B = trial_block(S, stim_name, areas, opt)
%TRIAL_BLOCK  One stimulus of one session in the agreed format.
%
%   B = trial_block(S, 'dot_motion', {'VISp'})
%   B = trial_block(S, 'dot_motion', {'VISp', 'LGd'}, struct('bin_s', 0.02, 'cell_norm', 'none'))
%
%   B.TrialsStim               M x 1   condition of each trial: a row of B.Conditions
%   B.Conditions               K rows  the stimulus parameters of each condition, e.g. Dir, Speed_deg_s
%   B.TrialID                  M x 1   1 = running, 2 = stationary, 0 = neither
%   B.Area.<name>.MUA          M x T   spikes/s summed over the area's units (before normalising)
%   B.Area.<name>.SingleCell   M x T x N, each unit normalised with opt.cell_norm
%   B.Area.<name>.units        N rows  which unit each page of SingleCell is (rows of S.units)
%   B.Trials                   M rows  everything else per trial: onset, parameters, wheel speed
%   B.t                        1 x T   bin centres, s from stimulus onset
%   B.opt                      the options used
%
%   opt fields (defaults in brackets):
%     window         [-0.5 2.5]   s relative to onset
%     bin_s          [0.01]       bin width (s)
%     smooth_method  ['gaussian'] 'none' | 'boxcar' | 'gaussian' | 'causal' (smooth_time.m)
%     smooth_s       [0.02]       Gaussian SD or boxcar width (s)
%     cell_norm      ['softnorm'] 'none' | 'sqrt' | 'zscore' | 'baseline' | 'softnorm' (normalise_cells.m)
%     softnorm_c     [5]          soft-norm constant (spikes/s)
%     norm_ref       ['all']      z-score / baseline reference bins: 'all', or 'pre' (before onset)
%     state_window   [[]]         s relative to onset over which running is judged; [] = window
%     state_opt      [struct()]   thresholds for trial_state.m, e.g. struct('run_cm_s', 5)
%     qc_only        [true]       only units with units.qc_pass
%
%   An area with no units in this session keeps its field, with MUA all NaN and SingleCell M x T x 0.
%   Trials overlapping an Allen invalid span (S.invalid) are NaN for the units on that probe.

o = struct('window', [-0.5 2.5], 'bin_s', 0.01, 'smooth_method', 'gaussian', 'smooth_s', 0.02, ...
    'cell_norm', 'softnorm', 'softnorm_c', 5, 'norm_ref', 'all', 'state_window', [], ...
    'state_opt', struct(), 'qc_only', true);
if nargin == 4
    for f = fieldnames(opt)'
        o.(f{1}) = opt.(f{1});
    end
end
if ~isfield(S.stim, stim_name)
    error('Session %d has no %s.', S.meta.session_id, stim_name);
end

% trials, conditions and behavioural state
[trials, conditions] = trial_table(S, stim_name);
state_window = o.state_window;
if isempty(state_window)
    state_window = o.window;
end
[trials.state, speed] = trial_state(S, trials.onset, state_window, o.state_opt);
trials = [trials, speed];

edges = o.window(1):o.bin_s:o.window(2);
t = edges(1:end - 1) + o.bin_s / 2;
if strcmp(o.norm_ref, 'pre')
    ref_bins = t < 0;
else
    ref_bins = true(size(t));
end

B.session_id = S.meta.session_id;
B.stimulus = stim_name;
B.TrialsStim = trials.condition;
B.Conditions = conditions;
B.TrialID = trials.state;
B.Trials = trials;
B.t = t;
B.opt = o;

units = struct2table(S.units);
M = height(trials);
for a = 1:numel(areas)
    r = find(strcmp(units.area, areas{a}) & (units.qc_pass | ~o.qc_only));
    if isempty(r)                                % no cells here: keep the field, but make it obviously empty
        B.Area.(areas{a}).MUA = nan(M, numel(t), 'single');
        B.Area.(areas{a}).SingleCell = zeros(M, numel(t), 0, 'single');
        B.Area.(areas{a}).units = units([], :);
        continue
    end
    counts = bin_spikes(S.spike_times(r), trials.onset, edges);

    % Allen invalid spans: NaN on the trials they overlap, for the units on that probe (-1 = every probe)
    for j = 1:numel(S.invalid.start_time)
        bad = trials.onset + o.window(2) > S.invalid.start_time(j) & ...
              trials.onset + o.window(1) < S.invalid.stop_time(j);
        on_probe = S.invalid.probe_id(j) == -1 | units.probe_id(r) == S.invalid.probe_id(j);
        counts(bad, :, on_probe) = NaN;
    end

    X = smooth_time(counts / o.bin_s, o.bin_s, o.smooth_method, o.smooth_s);   % spikes/s
    B.Area.(areas{a}).MUA = sum(X, 3);
    B.Area.(areas{a}).SingleCell = normalise_cells(X, o.cell_norm, ...
        struct('ref_bins', ref_bins, 'c', o.softnorm_c, 'condition', trials.condition));
    B.Area.(areas{a}).units = units(r, :);
end
end
