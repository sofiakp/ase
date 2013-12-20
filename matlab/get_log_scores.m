function get_log_scores()
% Computes the log-ratio of motif scores for all pairs of individuals. You
% must have scanned the motifs already and stored the scores for each
% parental haplotype of each individual in a separate file.
% The scores are combined over all subsequences that overlap the same
% region (see accum_idx).

indir = '/home/sofiakp/projects/Anshul/data/maya/temp'; % directory with fasta files and scores
inpref = 'Gm12878_allTFBS.sorted.noPol.'; % fasta files should be <inpref><indiv>.[mp]aternal.fa.
insuf = '.pouya.scores'; % Score files should be <inpref><indiv>.[mp]aternal<insuf>.mat
outdir = '/home/sofiakp/projects/Anshul/data/maya/temp/tfRegress_at_H3K27AC';

% Read all fastas and extract the individual names
cont = dir(fullfile(indir, '*.fa'));
indivs = unique(cellfun(@(x) regexprep(regexprep(x, '.(m|p)aternal.fa', ''), inpref, ''), {cont.name}, 'UniformOutput', false));
nindivs = length(indivs);
out_indivs = strrep(strrep(indivs, 'GM2', 'HG2'), 'SNYDER', 'MS1'); % Fix names

% These must be the same scores that you scanned on the fasta files
load('/home/sofiakp/projects/Anshul/data/motifs/hocomoco_pouya_stam_jaspar_pssms.mat');
sel = strcmp(pssm_source,'pouya');
pssms = pssms(sel);
% Maximum score that each PWM can achieve.
pwm_max = cellfun(@(x) sum(max(log(x), [], 1)), pssms, 'UniformOutput', true);

% This should have two columns. The first column is an index of a region
% (eg. an H3K27ac peak). The second is an index of one of the regions that
% were scanned. The scores of all scanned regions that overlap the same
% output region will be combined. So the output will have as many rows as
% unique values in the first column of accum_idx.
accum_idx = dlmread(fullfile(outdir, [inpref, 'idx_H3K27AC.txt']), '\t');
assert(size(accum_idx, 2) == 2);

for i = 1:nindivs
  scores1 = read_scores(indir, inpref, indivs{i}, insuf, accum_idx);
  assert(size(scores1, 2) == sum(sel));
  max_weight_mat = repmat(pwm_max', size(scores1, 1), 1);
  for j = (i + 1):nindivs
    disp([indivs{i}, '_vs_', indivs{j}]);
    scores2 = read_scores(indir, inpref, indivs{j}, insuf, accum_idx);
    log_scores = scores1 - scores2; % The scores are already in log space.
    % Weight the log-ratio of each motif instance by the maximum across the
    % two individuals. This will drive to 0 score differences where both
    % individuals have a really bad score.
    % The weight is 2/(1 + exp(-max_score+tot_max)), where max_score is the
    % maximum between the two individuals and tot_max is the maximum that
    % can be achieved for that PWM. The maximum value of the score
    % is 0, in which case the weight is 2/(1+exp(0)) = 1. The minimum value
    % of the score is -inf, in which case the weight -> 0.
    weights = 2 ./ (1 + exp(-max(scores1, scores2) + max_weight_mat)); 
    log_scores = log_scores .* weights; 
    outfile = fullfile(outdir, [inpref, out_indivs{i}, '_vs_', out_indivs{j}, insuf, '.txt']);
    dlmwrite(outfile, log_scores, 'delimiter', '\t', 'precision','%.3f');
  end
end
end

function accum_scores = read_scores(indir, inpref, indiv, insuf, accum_idx)
% Compute the average on the two haplotypes.
scores_struct = load(fullfile(indir, [inpref, indiv, '.maternal', insuf, '.mat']), 'scores');
scores = scores_struct.scores;
scores_struct = load(fullfile(indir, [inpref, indiv, '.paternal', insuf, '.mat']), 'scores');
scores = (scores + scores_struct.scores) / 2;

scores(scores == 0) = NaN; % Only series of Ns can give a score of 0. max(num, NaN) = num so these will be effectively ignored.

ntfs = size(scores, 2);
nseq = max(accum_idx(:, 1));
accum_scores = -inf(nseq, ntfs);
for t = 1:size(accum_idx, 1)
  accum_scores(accum_idx(t, 1), :) = max(accum_scores(accum_idx(t, 1), :), scores(accum_idx(t, 2), :));
end
accum_scores(accum_scores == 0) = NaN;

% This is too slow
%for t = 1:ntfs
%  accum_scores(:, t) = accumarray(accum_idx(:, 1), scores(accum_idx(:, 2), t), [nseq, 1], @(x) sum(x, 1), 0);
%end
end