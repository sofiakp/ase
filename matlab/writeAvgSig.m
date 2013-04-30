clear

MAYAROOT = '/net/isilon5/encode/work/akundaje/work/histoneVariation';
signal_dir = fullfile(MAYAROOT, 'rawdata/signal/combrep/extractSignal/rand/fc/avgSig/merged_Mar13/matFiles');
outdir = fullfile(MAYAROOT, 'rawdata/signal/combrep/extractSignal/rand/fc/avgSig/merged_Mar13/textFiles');
#outdir = fullfile(signal_dir, 'avgSig');
if ~isdir(outdir)
    mkdir(outdir);
end

signal_files = dir(fullfile(signal_dir, '*mat'));
disp(signal_files)

for i = 1:length(signal_files)
    iname = char(signal_files(i).name);
outfile = fullfile(outdir, strrep(iname, '.mat', '.txt'));
if ~exist(outfile, 'file')
    load(fullfile(signal_dir, iname));
    avg = cellfun(@nanmean, signal);
    dlmwrite(outfile, avg);
end
end
