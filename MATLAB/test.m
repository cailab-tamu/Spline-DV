% Spline-DV differential variability, using the analytic mean / CV / dropout
% curve (SC_ANALYTICFIT) in place of the fitted spline (SC_SPLINEFIT).
%
% See docs/analytic-curve.md for the derivation. Requires scGEAToolbox.

load GSM3308547_GSM3308548.mat

% SingleCellExperiment is a handle class and selectcells mutates the object in
% place, so copy before subsetting each sample. Without the copy, the second
% selectcells filters the already-narrowed object and both samples come back
% empty.
sce1 = copy(sce);
sce1 = sce1.selectcells(sce1.c_batch_id == "GSM3308547");
sce1 = sce1.qcfilter;

sce2 = copy(sce);
sce2 = sce2.selectcells(sce2.c_batch_id == "GSM3308548");
sce2 = sce2.qcfilter;

if ~isequal(sce1.g, sce2.g)
    [g_ori, ia, ib] = intersect(sce1.g, sce2.g, 'stable');
    X1_raw = full(sce1.X(ia, :));
    X2_raw = full(sce2.X(ib, :));
else
    g_ori = sce1.g;
    X1_raw = full(sce1.X);
    X2_raw = full(sce2.X);
end

% SC_ANALYTICFIT takes raw counts: it needs the library sizes to pin down
% alpha and the dropout law, and does the library-size normalization itself.
[T1, F1, curve1, par1] = sc_analyticfit(X1_raw, g_ori, SortIt = false);
[T2, F2, curve2, par2] = sc_analyticfit(X2_raw, g_ori, SortIt = false);

fprintf('sample 1: alpha = %.4f, phi = %.4f\n', par1.alpha, par1.phi);
fprintf('sample 2: alpha = %.4f, phi = %.4f\n', par2.alpha, par2.phi);

% Align both tables on gene name. F is in the internal order, so reorder it
% through nearidx before subsetting.
T1 = sortrows(T1, 'genes', 'ascend');   F1 = F1(T1.nearidx, :);
T2 = sortrows(T2, 'genes', 'ascend');   F2 = F2(T2.nearidx, :);

if ~isequal(T1.genes, T2.genes)         % a gene may be all-zero in one sample
    [~, k1, k2] = intersect(T1.genes, T2.genes, 'stable');
    T1 = T1(k1, :);  F1 = F1(k1, :);
    T2 = T2(k2, :);  F2 = F2(k2, :);
end
g = T1.genes;

% Normalized expression, rows aligned with g, for the profile plots in
% TEST_GUI.
X1_norm = sc_norm(X1_raw, 'type', 'libsize');
X2_norm = sc_norm(X2_raw, 'type', 'libsize');
[~, loc] = ismember(g, g_ori);
X1 = X1_norm(loc, :);
X2 = X2_norm(loc, :);

px1 = T1.lgu;  py1 = T1.lgcv;  pz1 = T1.dropr;
px2 = T2.lgu;  py2 = T2.lgcv;  pz2 = T2.dropr;

% Deviation of each gene from its own condition's reference curve.
v1 = [px1, py1, pz1] - F1;
v2 = [px2, py2, pz2] - F2;

DiffDist = vecnorm(v1 - v2, 2, 2);
DiffSign = sign(vecnorm(v2, 2, 2) - vecnorm(v1, 2, 2));

% Smooth samples of the two reference curves, for plotting in TEST_GUI.
mu1 = logspace(log10(expm1(min(px1))), log10(expm1(max(px1))), 400).';
mu2 = logspace(log10(expm1(min(px2))), log10(expm1(max(px2))), 400).';
xyz1 = [curve1.x(mu1), curve1.y(mu1), curve1.z(mu1)];
xyz2 = [curve2.x(mu2), curve2.y(mu2), curve2.z(mu2)];

T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s', 'Sample1'));
T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s', 'Sample2'));

T = [T1, T2, table(DiffDist), table(DiffSign)];

% Note: the spline version zeroed out DiffDist for genes whose nearest point
% was the first or last knot, because the fitted spline only spanned the
% observed genes. The analytic curve is defined for every mean > 0, so there
% is no out-of-range case left to patch.

T = sortrows(T, "DiffDist", "descend");

fprintf('%d genes; top DV: %s\n', height(T), strjoin(T.genes_Sample1(1:5), ', '));
