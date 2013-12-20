clearvars

MAYAROOT = '/media/fusion10/work/chromatinVariation';
signal_dir = fullfile(MAYAROOT, 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/extractSignal/');
outdir = signal_dir; %fullfile(signal_dir, 'avgSig');
if ~isdir(outdir)
    mkdir(outdir);
end

signal_files = dir(fullfile(signal_dir, 'SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC.mat'));
write_mean = false;

for i = 1:length(signal_files)
    iname = char(signal_files(i).name);
    load(fullfile(signal_dir, iname));
    outfile = fullfile(MAYAROOT, 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/extractSignal/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC_test.txt'); %fullfile(outdir, strrep(iname, '.mat', '.txt'));
    if ~exist(outfile, 'file')
        if write_mean
            sig_out = cellfun(@nanmean, signal);
            dlmwrite(outfile, avg);
        else
            dlmwrite(outfile, signal, 'delimiter', ',');
        end
    end
end