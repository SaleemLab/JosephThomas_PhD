%% V1 population PSTH for each dot speed, running vs stationary: per session and across sessions
% For each session, this builds the agreed format with trial_block.m (TrialsStim, TrialID,
% Area.<name>.SingleCell) and then uses only those fields. Each unit is soft-normed, so every unit counts
% equally in the population PSTH, and every session counts equally in the average across sessions.
%
% Press Run. It loads every session_<id>.mat in data_dir (on the server), one at a time, and takes about
% 5-10 minutes for all 26. The helper functions it uses are in helpers/. Skipped: the two sessions with a different dot-motion
% protocol (they have no stim.dot_motion) and sessions with no units in the chosen area.

%% Settings
here = fileparts(mfilename('fullpath'));   % this script's folder ('' = MATLAB's current folder)
addpath(fullfile(here, 'helpers'));       % the small functions this uses
data_dir = 'Y:\ibn-vision\DATA\ExternalDatasets\ABI_VisualCoding_Neuropixels\exports_for_joseph';   % the session_<id>.mat files (on the server): change if your drive letter differs
fig_dir = fullfile(here, 'figures', 'v1_speed_state');   % figures are saved here, next to this script; '' = leave them open
area = 'VISp';
stim_name = 'dot_motion';
speeds = [0 16 32 64 128 256 512];      % dot speeds, deg/s (B.Conditions.Speed_deg_s)
min_trials = 5;                         % a session's PSTH from fewer trials than this is left out (NaN)

opt = struct();                         % passed to trial_block; anything not set takes its default
opt.window = [-0.2 1.8];                % s from onset (the dots are on for 1 s)
opt.bin_s = 0.01;
opt.smooth_method = 'gaussian';
opt.smooth_s = 0.02;
opt.cell_norm = 'softnorm';             % x / (range + c), per unit
opt.softnorm_c = 5;                     % spikes/s
opt.state_window = [0 1];               % running / stationary judged while the dots are on
opt.qc_only = true;

%% Every session
files = dir(fullfile(data_dir, 'session_*.mat'));
files = files(~cellfun(@isempty, regexp({files.name}, '^session_\d+\.mat$')));
if ~isempty(fig_dir) && ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
nt = round(diff(opt.window) / opt.bin_s);
psth = nan(numel(files), numel(speeds), nt, 2);  % session x speed x time x state (1 = running, 2 = stationary)
n_trials = zeros(numel(files), numel(speeds), 2);
session_ids = nan(numel(files), 1);

for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name), 'units', 'spike_times', 'stim', 'running', 'meta', 'invalid');
    session_ids(i) = S.meta.session_id;
    if ~isfield(S.stim, stim_name)
        fprintf('%d: no %s, skipped\n', session_ids(i), stim_name);
        continue
    end
    B = trial_block(S, stim_name, {area}, opt);  % the agreed format
    clear S
    if isempty(B.Area.(area).units)
        fprintf('%d: no %s units, skipped\n', session_ids(i), area);
        continue
    end

    speed = B.Conditions.Speed_deg_s(B.TrialsStim);          % M x 1, the dot speed of each trial
    pop = mean(B.Area.(area).SingleCell, 3, 'omitnan');       % M x T, population: the mean over units
    for k = 1:numel(speeds)
        for s = 1:2
            sel = speed == speeds(k) & B.TrialID == s;
            n_trials(i, k, s) = sum(sel);
            if n_trials(i, k, s) >= min_trials
                psth(i, k, :, s) = mean(pop(sel, :), 1, 'omitnan');
            end
        end
    end
    fprintf('%d: %d %s units, running trials per speed %s\n', session_ids(i), ...
        height(B.Area.(area).units), area, mat2str(n_trials(i, :, 1)));

    % one figure for this session
    title_str = sprintf('Session %d, %s, %d units (soft-normed)', session_ids(i), area, height(B.Area.(area).units));
    plot_speed_state(B.t, speeds, squeeze(psth(i, :, :, :)), [], squeeze(n_trials(i, :, :)), 'trials', title_str);
    if ~isempty(fig_dir)
        exportgraphics(gcf, fullfile(fig_dir, sprintf('session_%d.png', session_ids(i))));
        close(gcf);
    end
end
t = B.t;

%% Across sessions: the mean of the session PSTHs, +- SEM over sessions
n_sessions = squeeze(sum(~isnan(psth(:, :, 1, :)), 1));     % speed x state
grand = squeeze(mean(psth, 1, 'omitnan'));                  % speed x time x state
sem = squeeze(std(psth, 0, 1, 'omitnan')) ./ sqrt(permute(n_sessions, [1 3 2]));
plot_speed_state(t, speeds, grand, sem, n_sessions, 'sessions', ...
    sprintf('%s across sessions, mean +- SEM over sessions (soft-normed)', area));
if ~isempty(fig_dir)
    exportgraphics(gcf, fullfile(fig_dir, 'across_sessions.png'));
end

% Left in the workspace: psth (session x speed x time x state), n_trials, session_ids, t, speeds


%% Plotting: one panel per speed, running (orange) and stationary (blue)
function plot_speed_state(t, speeds, P, E, n, n_word, title_str)
% P, E: speed x time x state (E = [] for no error band); n: speed x state counts shown in each title
colours = {[0.85 0.33 0.10], [0 0.45 0.74]};
figure('Position', [50 50 1500 320]);
ax = gobjects(1, numel(speeds));
for k = 1:numel(speeds)
    ax(k) = subplot(1, numel(speeds), k); hold on;
    for s = 1:2
        y = squeeze(P(k, :, s));
        if all(isnan(y))
            continue
        end
        if ~isempty(E)
            e = squeeze(E(k, :, s));
            fill([t, fliplr(t)], [y - e, fliplr(y + e)], colours{s}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        plot(t, y, 'Color', colours{s}, 'LineWidth', 1.2);
    end
    xline(0, ':'); xline(1, ':');                % dots on, dots off
    title({sprintf('%g deg/s', speeds(k)), sprintf('run %d, still %d %s', n(k, 1), n(k, 2), n_word)}, ...
        'FontWeight', 'normal');
    xlabel('s from onset'); box off;
    if k == 1
        ylabel('population rate (soft-normed)');
    end
end
linkaxes(ax, 'xy');
sgtitle([title_str, ': running (orange), stationary (blue)']);
end
