function writeDips(intervalData, step, dips, bounds, dip_sig, outfile)

nsig = size(intervalData, 1);
f = fopen(outfile, 'w');

for s = 1:nsig,
    tmp_dips = dips{s};
    for d = 1:size(tmp_dips, 1)
        % when step > 1, then index i with respect to the interval start is
        % in fact index step * (i - 1) + 1 with respect to the interval start.
        % This corresponds to absolute position
        % interval.start + ( step * (i - 1) + 1 ) - 1 = 
        % interval.start + step * (i - 1)
        % assuming the first position is numbered 1.
        start = intervalData.start(s) + step * (tmp_dips(d, 1) - 1);
        stop = intervalData.start(s) + step * (tmp_dips(d, 2) - 1);
        bound_pos = sprintf('\t%d', intervalData.start(s) + step * (bounds{s}(d, :) - 1));
        sig = sprintf('\t%.5f', dip_sig{s}(d, :));
        fprintf(f, '%s\t%d\t%d%s%s\n', char(intervalData.chr(s)), start, stop, bound_pos, sig);
    end
end
fclose(f);
end