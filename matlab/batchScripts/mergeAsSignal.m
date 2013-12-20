clearvars

% Merges extracted signal around the same set of SNPs for different
% individuals.

MAYAROOT = '/media/fusion10/work/chromatinVariation';
signal_dir = fullfile(MAYAROOT, 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/extractSignal');
outdir = fullfile(MAYAROOT, 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/extractSignal');
if ~isdir(outdir)
    mkdir(outdir);
end
outfile = fullfile(outdir, 'SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC.mat');

f = fopen(fullfile(MAYAROOT, 'rawdata/metadata/chromatinVariation_combrep_names_H3K27AC.tab'), 'r');
C = textscan(f, '%s%s%s%s');
fclose(f);
cell_line = unique(C{3}); % Cell lines to consider

% Signal files to read
signal_files = dir(fullfile(signal_dir, 'SNYDER_HG19_all_H3K27AC_AS_AT_SNYDER_HG19_*mat'));
nfiles = length({signal_files.name});

for i = 1:length(cell_line)
    % Find the signal file corresponding to this cell line (there should be
    % only one...)
    signal_file_ind = ~cellfun(@isempty, regexp({signal_files.name}, strcat('SNYDER_HG19_', cell_line(i)), 'once'));
    if sum(signal_file_ind) == 1
        dat = load(fullfile(signal_dir, char(signal_files(signal_file_ind).name)));
        [nsig, len_sig] = size(dat.signal);
        names =  strcat(cell_line(i), '_', arrayfun(@(x) {num2str(x)}, 1:nsig));
        if i == 1
            signal = zeros(nsig * nfiles, len_sig, 'single');
            all_names = cellstr(repmat('.', [nsig * nfiles, 1]));
        end
        signal(((i - 1) * nsig + 1):(i * nsig), :) = dat.signal;
        all_names(((i - 1) * nsig + 1):(i * nsig)) = names;
    end
end
intervalData = dataset();
intervalData.chr = nominal(repmat(dat.intervalData.chr, [nfiles, 1]));
intervalData.start = repmat(dat.intervalData.start, [nfiles, 1]);
intervalData.stop = repmat(dat.intervalData.stop, [nfiles, 1]);
intervalData.strand = repmat(dat.intervalData.strand, [nfiles, 1]);
intervalData.name = all_names;
save(outfile, 'intervalData', 'signal');