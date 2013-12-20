gerp_file = '/media/fusion10/work/chromatinVariation/rawdata/All_hg19_RS.mat';
mark = 'H3K4ME3';
indir = '/media/fusion10/work/chromatinVariation/rawdata/signal/combrep/extractSignal/fc/avgSig/plots/anova/';
pos_file = fullfile(indir, [mark, '_sig0.75Qt_matched_Fgt1_regionsWithF.txt']);
pos_outfile = fullfile(indir, [mark, '_sig0.75Qt_matched_Fgt1_gerp.mat']);
[pos_sig, pos_inter, ~] = extractSignal(pos_file, gerp_file, 'if', 'bed', 'sl', 0, 'sr', 0, 'ss', false, 'fw', false, ...
    'o', pos_outfile, 'ov', 'signal');

neg_file = fullfile(indir, [mark, '_sig0.75Qt_matched_Flt1_regionsWithF.txt']);
neg_outfile = fullfile(indir, [mark, '_sig0.75Qt_matched_Flt1_gerp.mat']);
[neg_sig, neg_inter, ~] = extractSignal(neg_file, gerp_file, 'if', 'bed', 'sl', 0, 'sr', 0, 'ss', false, 'fw', false, ...
    'o', neg_outfile, 'ov', 'signal');
%%

marks = {'H3K27AC', 'H3K4ME1', 'H3K4ME3'};

for m = 1:length(marks)
    load(fullfile(indir, [marks{m}, '_sig0.75Qt_matched_Fgt1_gerp.mat']), 'signal', 'intervalData');
    pos_sig = signal;
    pos_inter = intervalData;
    load(fullfile(indir, [marks{m}, '_sig0.75Qt_matched_Flt1_gerp.mat']), 'signal', 'intervalData');
    neg_sig = signal;
    neg_inter = intervalData;
    
    tot_pos_hits = cellfun(@(x) sum(x > 2), pos_sig);
    tot_neg_hits = cellfun(@(x) sum(x > 2), neg_sig);
    tot_pos_len = sum(pos_inter.stop - pos_inter.start + 1);
    tot_neg_len = sum(neg_inter.stop - neg_inter.start + 1);
    
    %pos_mean = cellfun(@(x) nanmedian(x), pos_sig);
    %neg_mean = cellfun(@(x) nanmedian(x), neg_sig);
    
    sum_bino_pval = binocdf(sum(tot_pos_hits > 0), length(tot_pos_hits), sum(tot_neg_hits > 0)/length(tot_neg_hits));
    tot_bino_pval = binocdf(sum(tot_pos_hits), sum(tot_pos_len), sum(tot_neg_hits)/sum(tot_neg_len));
    pos_sig_all = [pos_sig{:}];
    neg_sig_all = [neg_sig{:}];
    wilco_pval = ranksum(pos_sig_all, neg_sig_all);
    enrich = median(pos_sig_all)/median(neg_sig_all);
    fprintf('%s\t%g\t%g\t%g\t%g\t%g\t%g\n', marks{m}, sum_bino_pval, tot_bino_pval, wilco_pval, enrich, median(pos_sig_all), median(neg_sig_all));
end