function depth_matching_experiment(opts)
%DEPTH_MATCHING_EXPERIMENT  Should the deeper sample be downsampled first?
%
%   The two reference curves differ between samples mainly through alpha,
%   which is fixed by sequencing depth. That raises an obvious question: if
%   the deeper sample were thinned to match, the two curves would nearly
%   coincide, so would the DV result improve?
%
%   This script answers it by thinning sample 2 with a common binomial
%   probability p = HM1/HM2, where HM is the harmonic mean library size, and
%   comparing the resulting gene ranking against the full-depth one. It runs
%   the thinning under two independent seeds, so the ranking change caused by
%   the random thinning can be told apart from the change caused by matching
%   depth.
%
%   FINDING (bundled GSM3308547/8 data): no, do not downsample.
%
%     * alpha does converge, 1.5970 -> 1.998 against 1.9973 for sample 1, and
%       the gap between the two reference curves falls by about 60%.
%     * but two runs of the same procedure agree at only rho = 0.774, while a
%       thinned run agrees with the full-depth run at rho = 0.848. The
%       thinning injects more noise than the bias it removes.
%     * the top 100 genes are stable either way (92 to 93 shared), so the
%       damage is in the bulk of the ranking, where scores are near zero.
%
%   The useful by-product is that phi does NOT converge once depth is
%   matched: 0.152 for sample 1 against 0.204 for sample 2. With alpha equal,
%   that residual cannot be blamed on sequencing depth, so it is a real
%   global difference in biological overdispersion between the conditions.
%
%   NAME-VALUE OPTIONS
%     DataFile       .mat holding the SingleCellExperiment. Default is
%                    GSM3308547_GSM3308548.mat next to this file.
%     Seeds          random seeds for the thinning. Default [1 2]; at least
%                    two are needed for the noise floor to be measurable.
%     CompareSpline  also fit the splines, for a third point of comparison.
%                    Default true.
%
%   See also SC_ANALYTICFIT, SC_SPLINEFIT, TEST.

arguments
    opts.DataFile (1,1) string = ""
    opts.Seeds (1,:) double {mustBeInteger} = [1 2]
    opts.CompareSpline (1,1) logical = true
end

here = fileparts(mfilename("fullpath"));
if opts.DataFile == ""
    opts.DataFile = fullfile(here, "GSM3308547_GSM3308548.mat");
end
if ~isfile(opts.DataFile)
    error("depth_matching_experiment:NoData", "Data file not found: %s", opts.DataFile);
end
if numel(opts.Seeds) < 2
    error("depth_matching_experiment:NeedTwoSeeds", ...
        "At least two seeds are needed to measure the noise floor.");
end

batches = ["GSM3308547", "GSM3308548"];

% ---- raw, QC-filtered counts for the two samples ------------------------
% sce is a handle object and selectcells changes it in place, so copy first.
R = cell(1, 2);
gg = cell(1, 2);
for k = 1:2
    S = load(opts.DataFile, "sce");
    a = copy(S.sce);
    a = a.selectcells(a.c_batch_id == batches(k));
    a = a.qcfilter;
    R{k} = full(a.X);
    gg{k} = a.g;
    clear S a
end
if ~isequal(gg{1}, gg{2})
    [genes, ia, ib] = intersect(gg{1}, gg{2}, "stable");
    R{1} = R{1}(ia, :);
    R{2} = R{2}(ib, :);
else
    genes = gg{1};
end
R1 = R{1};  R2 = R{2};

L1 = sum(R1, 1).';
L2 = sum(R2, 1).';
hm1 = 1 / mean(1 ./ L1);
hm2 = 1 / mean(1 ./ L2);
p = hm1 / hm2;

fprintf('%d genes, %d and %d cells\n', numel(genes), size(R1, 2), size(R2, 2));
fprintf('harmonic mean library: %.0f (s1) vs %.0f (s2)  ->  thin s2 by p = %.4f\n\n', ...
    hm1, hm2, p);

% ---- baseline, full depth ----------------------------------------------
[dv0, sg0, par1, par2, dvGenes] = i_rundv(R1, R2, genes);
fprintf('full depth : alpha %.4f vs %.4f   phi %.4f vs %.4f\n', ...
    par1.alpha, par2.alpha, par1.phi, par2.phi);

% ---- thinned, one run per seed ------------------------------------------
nSeed = numel(opts.Seeds);
dvT = zeros(numel(dv0), nSeed);
sgT = zeros(numel(dv0), nSeed);
parT = cell(1, nSeed);
for s = 1:nSeed
    rng(opts.Seeds(s));
    R2t = i_thin(R2, p);
    [dvT(:, s), sgT(:, s), ~, parT{s}] = i_rundv(R1, R2t, genes);
    fprintf('seed %-4d  : alpha %.4f vs %.4f   phi %.4f vs %.4f\n', ...
        opts.Seeds(s), par1.alpha, parT{s}.alpha, par1.phi, parT{s}.phi);
end

% ---- how far apart are the two reference curves? ------------------------
mu = logspace(-3, log10(300), 4000).';
ycurve = @(al, ph) log1p(sqrt(al ./ mu + ph));
fprintf('\nmax |dy| between references: full depth %.4f  ->  depth matched %.4f\n', ...
    max(abs(ycurve(par1.alpha, par1.phi) - ycurve(par2.alpha, par2.phi))), ...
    max(abs(ycurve(par1.alpha, par1.phi) - ycurve(parT{1}.alpha, parT{1}.phi))));

% ---- what happened to the ranking? --------------------------------------
fprintf('\nSpearman\n');
fprintf('  seed %d vs seed %d      %.3f   <- noise floor of the thinning\n', ...
    opts.Seeds(1), opts.Seeds(2), corr(dvT(:,1), dvT(:,2), 'type', 'Spearman'));
fprintf('  seed %d vs full depth   %.3f\n', ...
    opts.Seeds(1), corr(dvT(:,1), dv0, 'type', 'Spearman'));

[~, ord0] = sort(dv0, 'descend');
[~, ord1] = sort(dvT(:,1), 'descend');
[~, ord2] = sort(dvT(:,2), 'descend');

if opts.CompareSpline
    dvS = i_splinedv(R1, R2, genes);
    fprintf('  full depth vs spline   %.3f\n', corr(dv0, dvS, 'type', 'Spearman'));
    fprintf('  seed %d vs spline       %.3f\n', ...
        opts.Seeds(1), corr(dvT(:,1), dvS, 'type', 'Spearman'));
    [~, ordS] = sort(dvS, 'descend');
end

fprintf('\ntop-100 overlap\n');
fprintf('  seed %d vs seed %d      %d\n', opts.Seeds(1), opts.Seeds(2), ...
    numel(intersect(ord1(1:100), ord2(1:100))));
fprintf('  seed %d vs full depth   %d\n', opts.Seeds(1), ...
    numel(intersect(ord1(1:100), ord0(1:100))));
if opts.CompareSpline
    fprintf('  seed %d vs spline       %d\n', opts.Seeds(1), ...
        numel(intersect(ord1(1:100), ordS(1:100))));
end

fprintf('\n%% of genes called more variable in sample 2\n');
fprintf('  full depth  %.1f%%\n', 100 * mean(sg0 > 0));
fprintf('  thinned    ');
fprintf(' %.1f%%', 100 * mean(sgT > 0, 1));
fprintf('\n');

fprintf('\ntop 10 full depth : %s\n', strjoin(dvGenes(ord0(1:10))', ', '));
fprintf('top 10 thinned    : %s\n', strjoin(dvGenes(ord1(1:10))', ', '));
end


% =========================================================================
function Rt = i_thin(R, p)
% Exact binomial thinning of a count matrix. Each read is kept with
% probability p: expand every count into its individual reads, keep at
% random, then re-aggregate. Same distribution as binornd(R, p), but about
% thirty times faster on a matrix this size.
Rt = R;
nz = Rt > 0;
v = Rt(nz);
readof = repelem((1:numel(v)).', v);
kept = rand(numel(readof), 1) < p;
Rt(nz) = accumarray(readof(kept), 1, [numel(v) 1]);
end


% =========================================================================
function [dv, sg, pa, pb, genesOut] = i_rundv(Ra, Rb, genelist)
% DV score for one pair of count matrices, using the analytic curve.
[Ta, Fa, ~, pa] = sc_analyticfit(Ra, genelist, SortIt = false);
[Tb, Fb, ~, pb] = sc_analyticfit(Rb, genelist, SortIt = false);
Ta = sortrows(Ta, 'genes');   Fa = Fa(Ta.nearidx, :);
Tb = sortrows(Tb, 'genes');   Fb = Fb(Tb.nearidx, :);
va = [Ta.lgu, Ta.lgcv, Ta.dropr] - Fa;
vb = [Tb.lgu, Tb.lgcv, Tb.dropr] - Fb;
dv = vecnorm(va - vb, 2, 2);
sg = sign(vecnorm(vb, 2, 2) - vecnorm(va, 2, 2));
genesOut = Ta.genes;   % the scores above are in THIS order, not the input order
end


% =========================================================================
function dv = i_splinedv(Ra, Rb, genelist)
% The same DV score, but from the fitted spline, for comparison.
ws = warning('off', 'all');
[Ta, ~, ~, ptsA] = sc_splinefit(sc_norm(Ra, 'type', 'libsize'), genelist, true, false);
[Tb, ~, ~, ptsB] = sc_splinefit(sc_norm(Rb, 'type', 'libsize'), genelist, true, false);
warning(ws);
Ta = sortrows(Ta, 'genes');
Tb = sortrows(Tb, 'genes');
va = [Ta.lgu, Ta.lgcv, Ta.dropr] - ptsA(Ta.nearidx, :);
vb = [Tb.lgu, Tb.lgcv, Tb.dropr] - ptsB(Tb.nearidx, :);
dv = vecnorm(va - vb, 2, 2);
end
