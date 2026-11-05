function [state, speed] = trial_state(S, onsets, window, opt)
%TRIAL_STATE  Running (1), stationary (2) or neither (0) on each trial, from the wheel speed.
%
%   state = trial_state(S, trials.onset, [0 2])
%   [state, speed] = trial_state(S, trials.onset, [0 2], struct('run_cm_s', 5))
%
%   window   s relative to onset over which the speed is judged
%   opt      any of these fields (defaults in brackets):
%              run_cm_s    [3]    running threshold (cm/s)
%              run_frac    [0.8]  running if the speed is above run_cm_s for at least this fraction of the window
%              still_cm_s  [1]    stationary threshold (cm/s)
%              still_frac  [0.8]  stationary if |speed| is below still_cm_s for at least this fraction
%              smooth_s    [1]    moving mean over the speed first: the raw wheel signal jitters by about
%                                 +-2 cm/s even when the mouse is still
%   state    M x 1: 1 = running, 2 = stationary, 0 = neither
%   speed    M x 3 table: mean_speed_cm_s, frac_running, frac_still, for applying your own rule instead

o = struct('run_cm_s', 3, 'run_frac', 0.8, 'still_cm_s', 1, 'still_frac', 0.8, 'smooth_s', 1);
if nargin == 4
    for f = fieldnames(opt)'
        o.(f{1}) = opt.(f{1});
    end
end

% the wheel is sampled irregularly: put it on a regular 50 Hz grid, then smooth
dt = 0.02;
[t, i] = unique(S.running.t(:));
v = S.running.speed_cm_s(i);
ok = isfinite(v);
tg = (t(find(ok, 1)):dt:t(find(ok, 1, 'last')))';
vg = interp1(t(ok), v(ok), tg);
vg = movmean(vg, max(1, round(o.smooth_s / dt)));

M = numel(onsets);
state = zeros(M, 1);
mean_speed_cm_s = nan(M, 1);
frac_running = nan(M, 1);
frac_still = nan(M, 1);
for m = 1:M
    in = tg >= onsets(m) + window(1) & tg < onsets(m) + window(2);
    if ~any(in)
        continue
    end
    mean_speed_cm_s(m) = mean(vg(in));
    frac_running(m) = mean(vg(in) > o.run_cm_s);
    frac_still(m) = mean(abs(vg(in)) < o.still_cm_s);
    if frac_running(m) >= o.run_frac
        state(m) = 1;
    elseif frac_still(m) >= o.still_frac
        state(m) = 2;
    end
end
speed = table(mean_speed_cm_s, frac_running, frac_still);
end
