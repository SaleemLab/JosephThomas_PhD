function [trials, conditions] = trial_table(S, stim_name)
%TRIAL_TABLE  One row per trial of a stimulus: onset, duration, parameters and condition number.
%
%   [trials, conditions] = trial_table(S, 'drifting_gratings_75_repeats')
%
%   trials.onset       stimulus onset (s, the same clock as the spike times)
%   trials.duration    how long the stimulus was on (s)
%   trials.<param>     every parameter that changes across this stimulus's trials,
%                      e.g. orientation and contrast for gratings, Dir and Speed for dot_motion
%   trials.condition   which row of conditions this trial belongs to
%   trials.block       Allen stimulus block (some stimuli are shown in more than one block)
%   conditions         one row per combination of the parameters, with n_trials
%
%   Natural movies: one trial per 30 s repeat (S.stim has one row per frame).

T = struct2table(S.stim.(stim_name));
if startsWith(stim_name, 'natural_movie')
    first = find(T.frame == 0);                  % first frame of each repeat
    last = [first(2:end) - 1; height(T)];        % last frame of each repeat
    T.stop_time(first) = T.stop_time(last);
    T = T(first, :);
end

trials = table(T.start_time, T.stop_time - T.start_time, 'VariableNames', {'onset', 'duration'});

% keep the parameters that change from trial to trial
skip = {'start_time', 'stop_time', 'stimulus_block', 'stimulus_index', 'stimulus_name', 'frame'};
params = setdiff(T.Properties.VariableNames, skip, 'stable');
for p = params
    v = T.(p{1});
    if isnumeric(v)
        varies = numel(unique(v(~isnan(v)))) + any(isnan(v)) > 1;   % NaN counts as one value
    else
        varies = numel(unique(v)) > 1;
    end
    if varies
        trials.(p{1}) = v;
    end
end
params = trials.Properties.VariableNames(3:end);

% one condition per combination of those parameters
if isempty(params)
    trials.condition = ones(height(trials), 1);
    conditions = table(1, height(trials), 'VariableNames', {'condition', 'n_trials'});
else
    % unique() treats every NaN as different, so give NaN a stand-in value (Inf) while grouping. Two sessions
    % have 15 dot_motion trials with no speed or direction (NaN), and they become one condition this way.
    key = trials(:, params);
    for p = params
        if isnumeric(key.(p{1}))
            key.(p{1})(isnan(key.(p{1}))) = Inf;
        end
    end
    [~, first, trials.condition] = unique(key);
    conditions = trials(first, params);
    conditions.n_trials = accumarray(trials.condition, 1);
    conditions = [table((1:height(conditions))', 'VariableNames', {'condition'}), conditions];
end

if ismember('stimulus_block', T.Properties.VariableNames)
    trials.block = T.stimulus_block;
end
end
