function y = normalise_cells(x, method, opt)
%NORMALISE_CELLS  Normalise each unit of trials x time x units data on its own.
%
%   y = normalise_cells(rate, 'zscore')
%   y = normalise_cells(rate, 'zscore', struct('ref_bins', t < 0))     mean and SD from the baseline only
%   y = normalise_cells(rate, 'softnorm', struct('c', 5, 'condition', trials.condition))
%
%   method   'none'      unchanged
%            'sqrt'      square root: makes the trial-to-trial variance more similar across firing rates
%            'zscore'    (x - mean) / SD, both taken over all trials and the bins in ref_bins
%            'baseline'  x - mean over ref_bins: the change from baseline, in the data's units (e.g. spikes/s)
%            'softnorm'  x / (range + c): range is the max - min of the unit's trial-averaged PSTH, taken
%                        within each condition if opt.condition is given (Churchland et al. 2012).
%                        c stops weakly firing units from being blown up
%   opt      any of these fields (defaults in brackets):
%              ref_bins   logical 1 x T, the bins that mean and SD come from           [every bin]
%              c          soft-norm constant, in the data's units                         [5, i.e. 5 spikes/s]
%              condition  M x 1 condition of each trial, for softnorm                     [all trials together]
%
%   A unit with SD 0 (no spikes in the reference bins) gives NaN for 'zscore'.
%   Ratios such as running / stationary only make sense with 'none', 'sqrt' or 'softnorm':
%   z-scored and baseline-subtracted values can be zero or negative.

o = struct('ref_bins', true(1, size(x, 2)), 'c', 5, 'condition', ones(size(x, 1), 1));
if nargin == 3
    for f = fieldnames(opt)'
        o.(f{1}) = opt.(f{1});
    end
end

switch method
    case 'none'
        y = x;
    case 'sqrt'
        y = sqrt(x);
    case {'zscore', 'baseline'}
        ref = x(:, o.ref_bins, :);
        mu = mean(ref, [1 2], 'omitnan');       % 1 x 1 x N: one mean per unit
        y = x - mu;
        if strcmp(method, 'zscore')
            sd = std(ref, 0, [1 2], 'omitnan');
            sd(sd == 0) = NaN;
            y = y ./ sd;
        end
    case 'softnorm'
        conds = unique(o.condition);
        P = nan(numel(conds), size(x, 2), size(x, 3), 'like', x);   % condition x time x unit PSTHs
        for k = 1:numel(conds)
            P(k, :, :) = mean(x(o.condition == conds(k), :, :), 1, 'omitnan');
        end
        range = max(P, [], [1 2]) - min(P, [], [1 2]);              % 1 x 1 x N
        y = x ./ (range + o.c);
    otherwise
        error('normalise_cells: method must be none, sqrt, zscore, baseline or softnorm');
end
end
