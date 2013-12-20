clearvars

% Finds dips in a list of datasets

MAYAROOT = '/media/fusion10/work/chromatinVariation';
mergedir = fullfile(MAYAROOT, 'rawdata/signal/combrep/peakFiles/merged/'); % Directory with merged peaks for the same mark
peakdir = fullfile(MAYAROOT, 'rawdata/signal/combrep/peakFiles'); % Directory with peak files for each individual
% Extracted signal at the merged peaks for each individual
signaldir = fullfile(MAYAROOT, 'rawdata/signal/combrep/extractSignal/llr');
outdir = fullfile(MAYAROOT, 'rawdata/signal/combrep/dips/llr');
if ~isdir(outdir)
    mkdir(outdir);
end
tmpdir = '/media/c300/tmp';

% Files with merged peaks (all peaks from all individuals for the same
% mark) have names
% SNYDER_HG19_<mark>_merged.encodePeak.gz
%
% Files with extracted signal around the merged peaks follow the format
% SNYDER_HG19_<mark>_merged_AT_<sname>.mat
%
% Files with individual-specific peaks start with the prefix
% <sname>_VS_<iname>
f = fopen(fullfile(MAYAROOT, 'rawdata/metadata/chromatinVariation_combrep_names.tab'), 'r');
C = textscan(f, '%s%s%s%s');
fclose(f);
sname = C{1};
iname = C{2};
cell_line = C{3};
mark = C{4};

peak_files = dir(peakdir);

for i = 1:length(sname)
    % Merged peaks
    merged_file = fullfile(mergedir, strcat('SNYDER_HG19_', char(mark(i)), '_merged.encodePeak.gz'));
    % Signal around the merged peaks
    signal_file = fullfile(signaldir, strcat('SNYDER_HG19_', char(mark(i)), '_merged_AT_', char(sname(i)), '.mat'));
    peak_file_ind = ~cellfun(@isempty, regexp({peak_files.name}, strcat(sname(i), '_VS_', iname(i)), 'once'));
    
    if exist(merged_file, 'file') && exist(signal_file, 'file') && sum(peak_file_ind) == 1
        % Individual-specific peaks
        peak_file = fullfile(peakdir, char(peak_files(peak_file_ind).name));
        tmp_file = fullfile(tmpdir, strcat('TMP_', char(peak_files(peak_file_ind).name), '.txt'));
        outpref = fullfile(outdir, strcat('SNYDER_HG19_', char(mark(i)), '_merged_AT_', char(sname(i)), '_dips'));
        
        disp(['Starting ', char(sname(i))]);
        
        % Find which merged peaks overlap with individual-specific peaks.
        % Merged peaks that don't overlap (i.e. are not present in that
        % individual) are excluded from dip calling.
        system(['intersectBed -a ', merged_file, ' -b ', peak_file, ' -c  | cut -f4 > ', tmp_file]);
        ov = dlmread(tmp_file);
        delete(tmp_file);
        
        % Load the signal around the merged peaks and call dips
        load(signal_file);
        nsig = size(signal, 1);
        dips = cell(nsig, 1);
        bounds = cell(nsig, 1);
        dip_sig = cell(nsig, 1);
        low_sig =  prctile(cellfun(@nanmean, signal(ov == 0)), 90);
        [dips(ov > 0), bounds(ov > 0), dip_sig(ov > 0)] = findDips(signal(ov > 0), [4, 80], low_sig, 60);
        
        % Write the results in mat and text format
        save(strcat(outpref, '.mat'), 'intervalData', 'dips', 'bounds', 'dip_sig');
        writeDips(intervalData, 10, dips, bounds, dip_sig, strcat(outpref, '.txt'));
    end
end