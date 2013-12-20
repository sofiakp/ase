clearvars

MAYAROOT = '/media/fusion10/work/chromatinVariation/';
infile = fullfile(MAYAROOT, 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/extractSignal/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC.mat');
outdir = fullfile(MAYAROOT, 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt');
[~, outpref, ~] = fileparts(infile);
outpref = [outpref, '_cagt_'];
low_sig = 2;
%%
cagt(infile, 'od', outdir, 'op', outpref, 'tt', 'asSNP', 'st', 'H3K27AC', ...
    'lowSignalCut', low_sig, 'lowVarCut', 1.0e-2, 'nanTreat', 'interpolate', ...
    'distance', 'correlation', 'avgFun', 'mean', 'k', 40, 'replicates', 10, 'maxiter', 500, ...
    'mergeDist', 0.2, 'overwrite', true);
%%
load(infile);
results = load(fullfile(outdir, [outpref, 'results.mat']));

oth = load('/media/fusion10/work/chromatinVariation/rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/extractSignal/SNYDER_HG19_GM12878_H3K27AC_AS_AT_SNYDER_HG19_GM12878_H3K4ME1.mat');
other_signals = oth.signal;
other_names = {'H3K4ME1'};
other_types = {'H3K4ME1'};

params.xrange = -500:10:500; 
params.merged = true;
params.signalType = 'H3K27AC';
makeSignalTable(fullfile(outdir, [outpref, 'clusterData_multi.txt']), results, signal, params, ...
    other_signals, other_names, other_types);