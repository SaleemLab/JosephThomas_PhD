function counts = bin_spikes(spike_times, onsets, edges)
%BIN_SPIKES  Spike counts per trial, time bin and unit.
%
%   counts = bin_spikes(S.spike_times(rows), trials.onset, -0.5:0.01:2.5)
%
%   spike_times  cell array with one vector of spike times (s) per unit, e.g. S.spike_times(rows)
%   onsets       M x 1 trial onsets (s)
%   edges        bin edges in s relative to onset (T + 1 edges give T bins)
%   counts       M x T x N spike counts: trials x time bins x units (single precision, exact for counts)
%
%   Divide by the bin width for a rate in Hz: counts / bin_s.

E = onsets(:) + edges(:)';                       % M x (T+1) bin edges in session time
[e, ~, j] = unique(E(:));                        % every edge once, in order
counts = zeros(numel(onsets), numel(edges) - 1, numel(spike_times), 'single');
for k = 1:numel(spike_times)
    % number of spikes before each edge; the count in a bin is the difference between its two edges
    n_before = cumsum(histcounts(spike_times{k}, [-Inf; e; Inf]));
    n_before = reshape(n_before(j), size(E));
    counts(:, :, k) = diff(n_before, 1, 2);
end
end
