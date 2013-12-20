clearvars
suf = '.jaspar.scores';
indir = '/home/sofiakp/projects/Anshul/data/maya/temp';
files = dir(fullfile(indir, '*.fa'));
nfiles = size(files, 1);

load('/home/sofiakp/projects/Anshul/data/motifs/hocomoco_pouya_stam_jaspar_pssms.mat');
sel = strcmp(pssm_source,'pouya');
pssms = pssms(sel);

pwms = cellfun(@(x) [log(x); zeros(1,size(x,2))], pssms, 'UniformOutput', false);
pwms_rc = cellfun(@pwmrc,pssms,'UniformOutput',false);
pwms_rc = cellfun(@(x) [log(x);zeros(1,size(x,2))], pwms_rc, 'UniformOutput', false); % pwm to pssm, and add 'N' position

for f = 1:nfiles
  filename = fullfile(indir, files(f).name);
  [~, base, ~] = fileparts(filename);
  disp(base);
  subseqs = loadfa(filename);
  outfile = fullfile(indir, [base, suf, '.txt']);
  scores = max(pssmscan(pwms, subseqs, 0), pssmscan(pwms_rc, subseqs, 0));
  save(fullfile(indir, [base, suf, '.mat']), 'scores');
  if f == 1
    prev_scores = scores;
    good_rows = false(size(prev_scores, 1), 1);
  else
    good_rows = good_rows | any(prev_scores ~= scores, 2);
  end
end

dlmwrite(fullfile(indir, 'Gm12878_allTFBS.sorted.noPol.diffRows.txt'), good_rows);
%%
for f = 1:nfiles
  filename = fullfile(indir, files(f).name);
  [~, base, ~] = fileparts(filename);
  disp(base);
  outfile = fullfile(indir, [base, suf, '.txt']);
  load(fullfile(indir, [base, suf, '.mat']), 'scores');
  scores = scores(good_rows, :);
  dlmwrite(outfile, scores, 'delimiter', '\t', 'precision','%.3f');
end
%%
clearvars 
indir = '/home/sofiakp/projects/Anshul/data/maya/temp';
inpref = 'Gm12878_allTFBS.sorted.noPol.';
insuf = '.jaspar.scores';
good_rows = dlmread(fullfile(indir, 'Gm12878_allTFBS.sorted.noPol.diffRows.txt'));

indivs = {'GM12878', 'GM10847'};
nindivs = length(indivs);
accum_idx = [1 2; 1 3; 2 4];

for i = 1:nindivs
  scores_struct = load(fullfile(indir, [inpref, indivs{i}, '.maternal', insuf, '.mat']), 'scores');
  scores1 = scores_struct;
  scores_struct = load(fullfile(indir, [inpref, indivs{i}, '.paternal', insuf, '.mat']), 'scores');
  scores1 = (scores1 + scores_struct) / 2;
  scores1 = scores1(good_rows, :);
%   for j = 1:nindivs
%     scores_struct = load(fullfile(indir, [inpref, indivs{j}, 'maternal', insuf, '.mat']), 'scores');
%     scores2 = scores_struct;
%     scores_struct = load(fullfile(indir, [inpref, indivs{j}, 'paternal', insuf, '.mat']), 'scores');
%     scores2 = (scores2 + scores_struct) / 2;
%     scores2 = scores2(good_rows, :);
%   end
  accumarray(accum_idx(:, 1), scores(accum_idx(:, 2), :), [max(accum_idx), 1], @(x) sum(x, 1), 0);
end
