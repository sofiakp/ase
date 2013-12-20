clearvars

% Extract signal around AS SNPs

MAYAROOT = '/media/fusion10/work/chromatinVariation';
signal_dir = fullfile(MAYAROOT, 'rawdata/signal/combrep/fc/matFiles');
snp_dir = fullfile(MAYAROOT, 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/text');
outdir = fullfile(MAYAROOT, 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/extractSignal');
if ~isdir(outdir)
    mkdir(outdir);
end

mark = 'CTCF'
win = 500;

f = fopen(fullfile(MAYAROOT, ['rawdata/metadata/chromatinVariation_combrep_names_', mark, '.tab']), 'r');
C = textscan(f, '%s%s%s%s');
fclose(f);
sname = C{1};
iname = C{2};
cell_line = C{3};
mark = C{4};

signal_files = dir(signal_dir);

for i = 1:length(sname)
    signal_file_ind = ~cellfun(@isempty, regexp({signal_files.name}, strcat(sname(i), '_VS_', iname(i)), 'once'));
    if sum(signal_file_ind) == 1
       signal_file = fullfile(signal_dir, char(signal_files(signal_file_ind).name));
       snp_file = fullfile(snp_dir, char(strcat('SNYDER_HG19_all_', mark, '_rep.hitInd.txt')));
       outfile = fullfile(outdir, char(strcat('SNYDER_HG19_all_', mark, '_AS_AT_', sname(i), '.mat')));
       extractSignal(snp_file, signal_file, 'if', 'summit', 'o', outfile, 'ov', 'signal',...
           'us', 'false', 'sl', win, 'sr', win, 'fw', true, 'mf', 'samplerate', 'mp', 10, 'ms', 100);
    end
end

