function [dips, bounds, dip_sig] = findDips(signal, dip_len, low_sig, high_signal_prc) 
%FINDDIPS finds dips in a set of signals.
%[DIPS, BOUNDS, DIP_SIG] = FINDDIPS(SIGNAL, dip_len, LOW_SIG, HIGH_SIGNAL_PRC)
% A dip is a local minimum with low enough signal (see below) surrounded by
% peaks with high enough signal.
% Local maxima must have signal > max(low_sig, p), where p is the
% high_signal_prc-th percentile of the signal in that region. Local minima
% must have signal < p. This removes very low local maxima and very high
% local minima. We don't require local minima to be < low_sig, since dips between
% very high peaks could have somewhat high signal.

nsig = size(signal, 1);
dips = cell(nsig, 1); % Dip start and end positions relative to the start of the corresponding region
bounds = cell(nsig, 1); % Positions of peaks surrounding each dip
dip_sig = cell(nsig, 1); % Mean signal at dips and surrounding peaks

for s = 1:nsig
    tmps = signal{s};
    high_signal_cut = prctile(tmps, high_signal_prc);
    [xmax, imax, xmin, imin] = extrema(tmps); % Find local maxima and minima
    sel_max = xmax > max(high_signal_cut, low_sig); % Select local maxima with high enough signal
    xmax = xmax(sel_max);
    imax = imax(sel_max);
    sel_min = xmin < high_signal_cut; 
    xmin = xmin(sel_min);
    imin = imin(sel_min);
    
    % Each dip must be surrounded by peaks, so if there's < 2 peaks,
    % there's no dip.
    if length(xmax) > 1
        nmin = length(imin);
        tmp_dips = zeros(nmin, 2);
        tmp_bounds = zeros(nmin, 2);
        tmp_sig = zeros(nmin, 3);
        
        % extrema returns maxima and minima sorted by height not by
        % position
        [imax, idx] = sort(imax, 'ascend');
        xmax = xmax(idx);
        [imin, idx] = sort(imin, 'ascend');
        xmin = xmin(idx);
        
        % For each local minimum (i.e. potential dip)
        for i = 1:nmin
            right_pos = find(imax > imin(i), 1, 'first'); % Peak immediately to the right of local minimum
            left_pos = find(imax < imin(i), 1, 'last'); % Peak immediately to the left of the local minimum
            if ~isempty(right_pos) && ~isempty(left_pos)
                % The dip region is the region with signal at most 
                % (min(signal_at_peaks)/4 higher than the signal at the
                % local minimum. This region might contain inside it local
                % 'bumps' (eg. small local maxima that were removed).
                cutoff = xmin(i) + (min(xmax(right_pos), xmax(left_pos)) - xmin(i)) / 4;
                dip_region = find(tmps(imax(left_pos):imax(right_pos)) < cutoff);
                
                reg_len = imax(right_pos) - imax(left_pos); % distance between neighboring peaks
                
                % Make sure the depleted region is not too short or too
                % long and has significantly lower signal than the adjacent
                % peaks.
                if ~isempty(dip_region) &&  reg_len > dip_len(1) && reg_len < dip_len(2) && ...
                        xmin(i) < min(xmax(right_pos), xmax(left_pos)) / 2
                    tmp_bounds(i, :) = [imax(left_pos), imax(right_pos)];
                    tmp_dips(i, :) = imax(left_pos) - 1 + [dip_region(1), dip_region(end)];
                    tmp_sig(i, :) = [nanmean(tmps(tmp_dips(i, 1):tmp_dips(i, 2))), ...
                        xmax(left_pos), xmax(right_pos)];
                end
            end
        end
        
        % Remove local minima that were not dips
        sel = all(tmp_dips > 0, 2);
        tmp_dips = tmp_dips(sel, :);
        tmp_bounds = tmp_bounds(sel, :);
        tmp_sig = tmp_sig(sel, :);
        
        % Due to the removal of some local maxima, we might end up with
        % duplicate dip regions.
        [bounds{s}, ui, ~] = unique(tmp_bounds, 'rows');
        dips{s} = tmp_dips(ui, :);
        dip_sig{s} = tmp_sig(ui, :);
    end
end
end