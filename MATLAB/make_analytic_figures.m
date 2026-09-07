function make_analytic_figures(opts)
%MAKE_ANALYTIC_FIGURES  Regenerate the figures in docs/ for the analytic curve.
%
%   MAKE_ANALYTIC_FIGURES() rebuilds both figures used in
%   docs/analytic-curve.md and docs/Spline-DV_analytic_curve.docx from the
%   bundled GSM3308547 / GSM3308548 data:
%
%     docs/analytic_vs_spline_3d.png   Figure 1, the 3-D gene cloud with the
%                                      fitted spline and the analytic curve
%     docs/analytic_vs_spline.png      Figure 2, the same comparison one
%                                      coordinate at a time
%
%   Requires scGEAToolbox on the path (for SC_SPLINEFIT, which supplies the
%   curve the analytic form is being compared against).
%
%   NAME-VALUE OPTIONS
%     DataFile    .mat file holding the SingleCellExperiment. Default is
%                 GSM3308547_GSM3308548.mat next to this file.
%     OutDir      where the PNGs are written. Default ../docs.
%     PhiSource   how the overdispersion is chosen:
%                 "spline" (default) fits phi by minimizing the distance to
%                     the fitted spline, which is the like-for-like setting
%                     used for the RMSE table in the documentation;
%                 "genes"  uses the robust L1 fit on the gene cloud that
%                     SC_ANALYTICFIT performs, i.e. what the shipped function
%                     actually does with no spline in the loop.
%                 Both values are printed either way.
%     Resolution  export resolution in DPI. Default 320.
%
%   See also SC_ANALYTICFIT, SC_SPLINEFIT.

arguments
    opts.DataFile (1,1) string = ""
    opts.OutDir (1,1) string = ""
    opts.PhiSource (1,1) string ...
        {mustBeMember(opts.PhiSource, ["spline", "genes"])} = "spline"
    opts.Resolution (1,1) double {mustBePositive} = 320
end

here = fileparts(mfilename("fullpath"));
if opts.DataFile == ""
    opts.DataFile = fullfile(here, "GSM3308547_GSM3308548.mat");
end
if opts.OutDir == ""
    opts.OutDir = fullfile(fileparts(here), "docs");
end
if ~isfile(opts.DataFile)
    error("make_analytic_figures:NoData", "Data file not found: %s", opts.DataFile);
end
if isempty(which("sc_splinefit"))
    error("make_analytic_figures:NoToolbox", ...
        "scGEAToolbox is not on the path; sc_splinefit is required.");
end
if ~isfolder(opts.OutDir)
    mkdir(opts.OutDir);
end

cfac   = 1e4;                       % CP10K scale used by pkg.norm_libsize
batches = ["GSM3308547", "GSM3308548"];

% ---- raw, QC-filtered counts for the two samples ------------------------
% SingleCellExperiment is a handle class, so reload for each sample rather
% than deriving both from one loaded object.
R = cell(1, 2);
g = cell(1, 2);
for k = 1:2
    S = load(opts.DataFile, "sce");
    a = S.sce.selectcells(S.sce.c_batch_id == batches(k));
    a = a.qcfilter;
    R{k} = full(a.X);
    g{k} = a.g;
    clear S a
end
if ~isequal(g{1}, g{2})
    [genes, ia, ib] = intersect(g{1}, g{2}, "stable");
    R{1} = R{1}(ia, :);
    R{2} = R{2}(ib, :);
else
    genes = g{1};
end
fprintf("%d common genes; %d and %d cells\n", numel(genes), size(R{1}, 2), size(R{2}, 2));

% ---- per-sample statistics, spline curve, analytic parameters -----------
S = struct("xyz", {[], []}, "xyz1", {[], []}, "libSize", {[], []}, ...
           "alpha", {[], []}, "phi", {[], []});
for k = 1:2
    libSize = sum(R{k}, 1).';
    Xn = (R{k} ./ libSize.') * cfac;

    ws = warning("off", "all");
    [Tb, ~, ~, xyz1] = sc_splinefit(Xn, genes, true, false);
    warning(ws);
    Tb = sortrows(Tb, "genes", "ascend");

    alpha = cfac * mean(1 ./ libSize);
    u = expm1(Tb.lgu);

    % phi fitted to the spline, and phi fitted robustly to the genes
    % phi enters only the CV branch here (the dropout branch is the Poisson
    % limit, phi_z = 0), so it is fitted against the CV coordinate alone.
    [xu, iu] = unique(xyz1(:, 1));
    muS = expm1(xu);  ysp = xyz1(iu, 2);
    objSpline = @(lp) sum((log1p(sqrt(alpha ./ muS + 10^lp)) - ysp).^2);
    phiSpline = 10 ^ fminbnd(objSpline, -6, 1);

    ok = isfinite(u) & isfinite(Tb.lgcv) & u > 0;
    objGenes = @(lp) sum(abs(log1p(sqrt(alpha ./ u(ok) + 10^lp)) - Tb.lgcv(ok)));
    phiGenes = 10 ^ fminbnd(objGenes, -6, 1);

    fprintf("sample %d: alpha = %.4f   phi(spline) = %.4f   phi(genes, L1) = %.4f\n", ...
        k, alpha, phiSpline, phiGenes);

    if opts.PhiSource == "spline"
        S(k).phi = phiSpline;
    else
        S(k).phi = phiGenes;
    end
    S(k).xyz     = [Tb.lgu, Tb.lgcv, Tb.dropr];
    S(k).xyz1    = xyz1;
    S(k).libSize = libSize;
    S(k).alpha   = alpha;
end

% ---- Figure 1: 3-D ------------------------------------------------------
cSpline = [0.85 0.33 0.10];
cAnalyt = [0.00 0.45 0.74];

fig = figure("Position", [60 60 940 440], "Color", "w", "Visible", "off");
tl = tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for k = 1:2
    xyz = S(k).xyz;  xyz1 = S(k).xyz1;
    ax = nexttile(tl);  hold(ax, "on");  grid(ax, "on");  box(ax, "off");

    scatter3(ax, xyz(:,1), xyz(:,2), xyz(:,3), 8, [0.52 0.55 0.62], ...
        "filled", "MarkerFaceAlpha", 0.16, "MarkerEdgeColor", "none");

    xl = [0 max(xyz(:,1)) * 1.02];
    yl = [0 max(xyz(:,2)) * 1.05];

    % the plane the dropout rate may not cross
    patch(ax, "XData", [xl(1) xl(2) xl(2) xl(1)], ...
              "YData", [yl(1) yl(1) yl(2) yl(2)], "ZData", [0 0 0 0], ...
        "FaceColor", [0.2 0.2 0.2], "FaceAlpha", 0.06, ...
        "EdgeColor", [0.5 0.5 0.5], "LineStyle", ":");

    [cx, cy, cz] = i_curve(xyz(:,1), S(k).libSize, cfac, S(k).alpha, S(k).phi);

    hS = plot3(ax, xyz1(:,1), xyz1(:,2), xyz1(:,3), "-",  "LineWidth", 4.5, "Color", cSpline);
    hA = plot3(ax, cx, cy, cz,                      "--", "LineWidth", 3.0, "Color", cAnalyt);
    hG = scatter3(ax, nan, nan, nan, 24, [0.45 0.47 0.55], "filled");

    bad = xyz1(:,3) < 0;
    if any(bad)
        plot3(ax, xyz1(bad,1), xyz1(bad,2), xyz1(bad,3), "-", ...
            "LineWidth", 6, "Color", [0.70 0 0]);
        [~, ib] = min(xyz1(:,3));
        text(ax, xyz1(ib,1) + 0.30, xyz1(ib,2), xyz1(ib,3) + 0.12, ...
            sprintf("spline ends here, at z = %.3f", xyz1(ib,3)), ...
            "Color", [0.70 0 0], "FontSize", 10, "FontWeight", "bold", ...
            "HorizontalAlignment", "left");
    end
    plot3(ax, xyz1(end,1), xyz1(end,2), xyz1(end,3), "o", "MarkerSize", 7, ...
        "MarkerFaceColor", "w", "MarkerEdgeColor", cSpline, "LineWidth", 1.5);

    xlim(ax, xl);  ylim(ax, yl);  zlim(ax, [-0.2 1]);
    xlabel(ax, "log(1 + mean)");  ylabel(ax, "log(1 + CV)");
    zlabel(ax, "dropout rate");
    title(ax, sprintf("Sample %d   (\\alpha = %.3f, \\phi = %.3f)", ...
        k, S(k).alpha, S(k).phi), "FontWeight", "normal");
    view(ax, -37, 18);
    set(ax, "FontSize", 10);
    if k == 1
        lg = legend(ax, [hG hS hA], {'genes', 'fitted spline', 'analytic curve'}, ...
            "Box", "off", "FontSize", 10);
        lg.Position(1:2) = [0.285 0.695];
    end
end
f3d = fullfile(opts.OutDir, "analytic_vs_spline_3d.png");
exportgraphics(fig, f3d, "Resolution", opts.Resolution);
close(fig);
fprintf("wrote %s\n", f3d);

% ---- Figure 2: coordinate by coordinate ---------------------------------
fig = figure("Position", [80 80 1200 480], "Color", "w", "Visible", "off");
for k = 1:2
    xyz = S(k).xyz;  xyz1 = S(k).xyz1;
    [xu, iu] = unique(xyz1(:,1));
    [~, cy, cz] = i_curve(xu, S(k).libSize, cfac, S(k).alpha, S(k).phi, true);

    subplot(2, 2, k);
    plot(xyz(:,1), xyz(:,2), ".", "Color", [.8 .8 .85], "MarkerSize", 3);  hold on
    plot(xu, xyz1(iu,2), "-",  "LineWidth", 2.5, "Color", cSpline);
    plot(xu, cy,         "--", "LineWidth", 2.0, "Color", cAnalyt);
    xlabel("log(1+mean)");  ylabel("log(1+CV)");
    title(sprintf("Sample %d: CV", k));
    legend({'genes', 'spline', 'analytic'}, "Location", "northeast");  box off

    subplot(2, 2, k + 2);
    plot(xyz(:,1), xyz(:,3), ".", "Color", [.8 .8 .85], "MarkerSize", 3);  hold on
    plot(xu, xyz1(iu,3), "-",  "LineWidth", 2.5, "Color", cSpline);
    plot(xu, cz,         "--", "LineWidth", 2.0, "Color", cAnalyt);
    yline(0, ":k");
    xlabel("log(1+mean)");  ylabel("dropout rate");
    title(sprintf("Sample %d: dropout (0 fitted params)", k));  box off
end
f2d = fullfile(opts.OutDir, "analytic_vs_spline.png");
exportgraphics(fig, f2d, "Resolution", 150);
close(fig);
fprintf("wrote %s\n", f2d);
end

% =========================================================================
function [cx, cy, cz] = i_curve(lgu, libSize, cfac, alpha, phi, atX)
% Analytic curve. Sampled over the full range of gene means, or, when ATX is
% true, evaluated exactly at the supplied log(1+mean) values.
if nargin < 6, atX = false; end
if atX
    mu = expm1(lgu(:));
else
    lo = max(expm1(min(lgu)), 1e-4);
    mu = logspace(log10(lo), log10(expm1(max(lgu))), 800).';
end
cx = log1p(mu);
cy = log1p(sqrt(alpha ./ mu + phi));
cz = i_dropout(mu, libSize, cfac, 0);   % dropout uses the Poisson limit
end

% =========================================================================
function z = i_dropout(mu, libSize, cfac, phi)
mu  = mu(:).';
z   = zeros(numel(mu), 1);
blk = max(1, floor(2e7 / numel(libSize)));
for i = 1:blk:numel(mu)
    j = i:min(i + blk - 1, numel(mu));
    t = libSize(:) * (mu(j) / cfac);
    if phi <= 1e-10
        z(j) = mean(exp(-t), 1).';
    else
        z(j) = mean((1 + phi * t) .^ (-1/phi), 1).';
    end
end
end
