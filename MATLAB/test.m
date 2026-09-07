% Spline-DV demo: two ways to build the reference curve, side by side.
%
% Every gene is placed at a point (log mean, log CV, dropout rate). A
% reference curve is drawn through that cloud for each sample, and a gene's
% deviation from the curve measures how variable it is. Comparing the two
% deviation vectors gives the DV score.
%
%   METHOD 1  sc_splinefit    a smoothing spline fitted through the genes
%   METHOD 2  sc_analyticfit  a closed-form curve from a gamma-Poisson model
%                             (derived in docs/analytic-curve.md)
%
% The script runs both and compares the gene rankings they produce.
% Requires scGEAToolbox on the path.


%% 1. Load the data and split it into the two samples

load GSM3308547_GSM3308548.mat

% NOTE: sce is a handle object, so selectcells changes it in place instead of
% returning a new one. Copy it first, or filtering sample 1 also destroys
% sample 2 and both come out empty.
sce1 = copy(sce);
sce1 = sce1.selectcells(sce1.c_batch_id == "GSM3308547");
sce1 = sce1.qcfilter;

sce2 = copy(sce);
sce2 = sce2.selectcells(sce2.c_batch_id == "GSM3308548");
sce2 = sce2.qcfilter;

% Keep only the genes that survived QC in both samples.
if ~isequal(sce1.g, sce2.g)
    [g_ori, ia, ib] = intersect(sce1.g, sce2.g, 'stable');
    X1_raw = full(sce1.X(ia, :));
    X2_raw = full(sce2.X(ib, :));
else
    g_ori  = sce1.g;
    X1_raw = full(sce1.X);
    X2_raw = full(sce2.X);
end

fprintf('%d genes, %d cells in sample 1, %d cells in sample 2\n\n', ...
    numel(g_ori), size(X1_raw, 2), size(X2_raw, 2));


%% 2. METHOD 1 - the fitted spline

% sc_splinefit expects library-size normalized data.
X1_norm = sc_norm(X1_raw, 'type', 'libsize');
X2_norm = sc_norm(X2_raw, 'type', 'libsize');

[S1, ~, ~, splinepts1] = sc_splinefit(X1_norm, g_ori, true, false);
[S2, ~, ~, splinepts2] = sc_splinefit(X2_norm, g_ori, true, false);

% Put both tables in the same gene order so the rows line up.
S1 = sortrows(S1, 'genes', 'ascend');
S2 = sortrows(S2, 'genes', 'ascend');

% S.nearidx points at the closest curve point for each gene, so subtracting
% it from the gene's own position gives the deviation vector.
vS1 = [S1.lgu, S1.lgcv, S1.dropr] - splinepts1(S1.nearidx, :);
vS2 = [S2.lgu, S2.lgcv, S2.dropr] - splinepts2(S2.nearidx, :);

DiffDist_spline = vecnorm(vS1 - vS2, 2, 2);
DiffSign_spline = sign(vecnorm(vS2, 2, 2) - vecnorm(vS1, 2, 2));

% The spline only covers the range of genes actually seen, so a gene sitting
% on its very first or last point is not really on the curve at all. The
% original method drops those genes by setting their score to zero.
atEnd = S1.nearidx == 1 | S1.nearidx == max(S1.nearidx) | ...
        S2.nearidx == 1 | S2.nearidx == max(S2.nearidx);
DiffDist_spline(atEnd) = 0;

fprintf('spline:   %d genes zeroed at the ends of the curve\n', sum(atEnd));


%% 3. METHOD 2 - the analytic curve

% sc_analyticfit takes RAW counts, because it needs the library sizes to work
% out alpha and the dropout law. It normalizes internally.
[A1, foot1, curve1, par1] = sc_analyticfit(X1_raw, g_ori, SortIt = false);
[A2, foot2, curve2, par2] = sc_analyticfit(X2_raw, g_ori, SortIt = false);

fprintf('analytic: sample 1  alpha = %.4f  phi = %.4f\n', par1.alpha, par1.phi);
fprintf('analytic: sample 2  alpha = %.4f  phi = %.4f\n', par2.alpha, par2.phi);

% Same gene order as above. foot is the point on the curve closest to each
% gene, so it plays the role splinepts(nearidx,:) played in method 1.
A1 = sortrows(A1, 'genes', 'ascend');   foot1 = foot1(A1.nearidx, :);
A2 = sortrows(A2, 'genes', 'ascend');   foot2 = foot2(A2.nearidx, :);

vA1 = [A1.lgu, A1.lgcv, A1.dropr] - foot1;
vA2 = [A2.lgu, A2.lgcv, A2.dropr] - foot2;

DiffDist_analytic = vecnorm(vA1 - vA2, 2, 2);
DiffSign_analytic = sign(vecnorm(vA2, 2, 2) - vecnorm(vA1, 2, 2));

% No end correction is needed here. The analytic curve is defined for every
% mean value, so no gene can fall off the end of it.


%% 4. Collect both results in one table

assert(isequal(S1.genes, A1.genes), 'the two methods returned different genes')
genes = S1.genes;

T = table(genes, ...
    S1.lgu, S1.lgcv, S1.dropr, S2.lgu, S2.lgcv, S2.dropr, ...
    DiffDist_spline,   DiffSign_spline, ...
    DiffDist_analytic, DiffSign_analytic, ...
    'VariableNames', {'genes', ...
    'lgu_1', 'lgcv_1', 'dropr_1', 'lgu_2', 'lgcv_2', 'dropr_2', ...
    'DiffDist_spline', 'DiffSign_spline', ...
    'DiffDist_analytic', 'DiffSign_analytic'});


%% 5. Compare the two rankings

[~, rank_spline]   = sort(T.DiffDist_spline,   'descend');
[~, rank_analytic] = sort(T.DiffDist_analytic, 'descend');

fprintf('\nSpearman correlation of the two DV scores: %.3f\n', ...
    corr(T.DiffDist_spline, T.DiffDist_analytic, 'type', 'Spearman'));

for n = [100 500 1000]
    shared = numel(intersect(rank_spline(1:n), rank_analytic(1:n)));
    fprintf('top %-4d genes shared: %d\n', n, shared);
end

fprintf('\ntop 10, spline  : %s\n', strjoin(T.genes(rank_spline(1:10))',   ', '));
fprintf('top 10, analytic: %s\n\n', strjoin(T.genes(rank_analytic(1:10))', ', '));


%% 6. Draw the two curves on the gene cloud of sample 1

% Sample the analytic curve on a grid of mean values so it plots as a line.
mu = logspace(log10(expm1(min(T.lgu_1))), log10(expm1(max(T.lgu_1))), 400)';
analyticpts1 = [curve1.x(mu), curve1.y(mu), curve1.z(mu)];

figure('Color', 'w');
scatter3(T.lgu_1, T.lgcv_1, T.dropr_1, 8, [.6 .6 .65], 'filled', ...
    'MarkerFaceAlpha', 0.15);
hold on
plot3(splinepts1(:, 1), splinepts1(:, 2), splinepts1(:, 3), '-', ...
    'LineWidth', 4, 'Color', [0.85 0.33 0.10]);
plot3(analyticpts1(:, 1), analyticpts1(:, 2), analyticpts1(:, 3), '--', ...
    'LineWidth', 3, 'Color', [0.00 0.45 0.74]);
hold off
grid on
view(-37, 18)
xlabel('log(1 + mean)')
ylabel('log(1 + CV)')
zlabel('dropout rate')
title('Sample 1: reference curve, fitted vs analytic')
legend({'genes', 'spline', 'analytic'}, 'Location', 'northeast')


%% 7. Variables for test_gui

% test_gui plots the gene clouds, the two reference curves, and the
% expression profile of a selected gene. It expects these names. Swap
% splinepts for analyticpts below to view the analytic curves instead.
g = genes;

[~, loc] = ismember(genes, g_ori);
X1 = X1_norm(loc, :);
X2 = X2_norm(loc, :);

px1 = T.lgu_1;  py1 = T.lgcv_1;  pz1 = T.dropr_1;
px2 = T.lgu_2;  py2 = T.lgcv_2;  pz2 = T.dropr_2;

xyz1 = splinepts1;
xyz2 = splinepts2;

% Sort the table for reading and export; g, px1, ... stay in gene order.
T = sortrows(T, 'DiffDist_spline', 'descend');
