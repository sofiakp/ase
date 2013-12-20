clearvars

MAYAROOT = '/media/fusion10/work/chromatinVariation';
dip_file = fullfile(MAYAROOT, 'rawdata/signal/combrep/dips/llr/bed/SNYDER_HG19_H3K27AC_merged_dips.bed');
bedpref = 'SNYDER_HG19_H3K27AC_merged_dips';
signal_dir = fullfile(MAYAROOT, 'rawdata/signal/combrep/llr/matFiles');
outdir = fullfile(MAYAROOT, 'rawdata/signal/combrep/dips/llr/avgSig');
if ~isdir(outdir)
    mkdir(outdir);
end

f = fopen(fullfile(MAYAROOT, 'rawdata/metadata/chromatinVariation_combrep_names_H3K27AC.tab'), 'r');
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
       outpref = fullfile(outdir, strcat(bedpref, '_AT_', char(sname(i)), '_avgSig'));
       tmp_file = [outpref, '_TMP.mat'];
       disp(outpref);
       
       % Average signal in dip region
       [sig, intervalData, ~] = extractSignal(dip_file, signal_file, 'if', 'bed', 'o', [outpref, '.mat'], 'ov', 'signal',...
           'us', 'false', 'sl', 0, 'sr', 0, 'mf', 'mean');
       dip_avg = zeros(size(sig, 1), 3);
       dip_avg(:, 1) = cell2mat(sig);
       
       % Average signal upstream of the dip region
       intervalData_cp = intervalData;
       intervalData.stop = intervalData_cp.start - 1;
       intervalData.start = intervalData_cp.start - 200;
       save(tmp_file, 'intervalData');
       [sig, ~, ~] = extractSignal(tmp_file, signal_file, 'if', 'mat', 'o', [outpref, '.mat'], 'ov', 'signal',...
           'us', 'false', 'sl', 0, 'sr', 0, 'fw', true, 'mf', 'mean');
       dip_avg(:, 2) = sig;
       
       % Interval data downstream of the dip region
       intervalData = intervalData_cp;
       intervalData.stop = intervalData_cp.stop + 200;
       intervalData.start = intervalData_cp.stop + 1;
       save(tmp_file, 'intervalData');
       [sig, ~, ~] = extractSignal(tmp_file, signal_file, 'if', 'mat', 'o', [outpref, '.mat'], 'ov', 'signal',...
           'us', 'false', 'sl', 0, 'sr', 0, 'fw', true, 'mf', 'mean');
       dip_avg(:, 3) = sig;
       
       % Restore the original dip regions
       intervalData = intervalData_cp;
       save([outpref, '.mat'], 'intervalData', 'dip_avg');
       dlmwrite([outpref, '.txt'], dip_avg, 'delimiter', '\t', 'precision', '%.5f');
    end
end