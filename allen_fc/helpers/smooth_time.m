function y = smooth_time(x, bin_s, method, width_s)
%SMOOTH_TIME  Smooth binned data along time (dimension 2), e.g. trials x time x units.
%
%   y = smooth_time(rate, 0.01, 'gaussian', 0.02)    Gaussian with a 20 ms SD, for 10 ms bins
%
%   method   'none'      no smoothing
%            'boxcar'    moving mean over width_s
%            'gaussian'  Gaussian kernel with SD width_s, cut at +-3 SD
%            'causal'    the past half of that Gaussian only: each bin averages itself and earlier bins,
%                        so a response never appears to start before it really does
%
%   At the edges of the window the kernel is rescaled over the bins that exist, so the first and last
%   bins are not pulled towards zero. Smooth single trials first and then average, or average first:
%   the mean comes out the same.

if strcmp(method, 'none')
    y = x;
    return
end
w = width_s / bin_s;                             % width in bins
switch method
    case 'boxcar'
        h = floor(round(w) / 2);                 % an odd number of bins keeps the kernel centred
        k = ones(1, 2 * h + 1);
    case {'gaussian', 'causal'}
        h = ceil(3 * w);
        lag = -h:h;
        k = exp(-lag .^ 2 / (2 * w ^ 2));
        if strcmp(method, 'causal')
            k(lag < 0) = 0;                      % y(t) = sum over lag of k(lag) * x(t - lag)
        end
    otherwise
        error('smooth_time: method must be none, boxcar, gaussian or causal');
end
y = convn(x, k, 'same') ./ convn(ones(1, size(x, 2)), k, 'same');
end
