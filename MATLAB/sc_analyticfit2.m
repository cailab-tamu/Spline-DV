function [T, sx, sy, sz, d, params] = sc_analyticfit2(X, Y, genelistx, genelisty, sortid, opts)
%SC_ANALYTICFIT2  Compare two conditions using the analytic reference curve.
%
%   The analytic counterpart of SC_SPLINEFIT2: it intersects the two gene
%   lists, fits each condition's reference curve, joins the two per-gene
%   tables, and scores every gene by dd = d2 - d1, the change in how far it
%   sits from its own condition's curve.
%
%   [T, sx, sy, sz, d] = SC_ANALYTICFIT2(X, Y, genelistx, genelisty, sortid)
%
%   INPUT CONVENTION DIFFERS FROM SC_SPLINEFIT2. That function takes
%   library-size normalized matrices; this one takes RAW counts, because the
%   analytic curve needs the library sizes to pin down alpha and the dropout
%   law, and it normalizes internally. To pass data that is already
%   normalized, set IsNormalized=true and supply LibSizeX and LibSizeY.
%
%   OUTPUTS
%     T       joined table: genes, logu1, logcv1, dropr1, d1, pval1, fdr1,
%             nearidx1, then the same for condition 2, then dd = d2 - d1.
%             Sorted by dd, descending, unless sortid is false.
%     sx, sy  foot point on each condition's curve for each of its genes.
%             As in SC_SPLINEFIT2, sx(T.nearidx1,:) is the foot point of the
%             genes in T, and likewise sy(T.nearidx2,:).
%     sz, d   Procrustes comparison of the two curves, as in SC_SPLINEFIT2:
%             d is the dissimilarity and sz the transformed second curve.
%             Empty and NaN if the two conditions kept different numbers of
%             genes, which Procrustes cannot handle.
%     params  struct with alpha1, phi1, alpha2, phi2 and the two curve
%             handles. Because the curves are analytic, this supports an
%             exact comparison that Procrustes only approximates -- see the
%             note below.
%
%   NAME-VALUE OPTIONS
%     LibSizeX, LibSizeY  raw library sizes, required when IsNormalized
%     ScaleFactor         normalization scale, default 1e4
%     IsNormalized        default false
%     Dispersion          [] to fit phi per condition (default), or a scalar
%                         to force the same phi on both, which removes the
%                         overdispersion difference from the reference
%     OneSided            default true, matching SC_SPLINEFIT
%
%   ON COMPARING THE TWO CURVES
%   Procrustes treats each curve as an unordered cloud and finds the best
%   rigid alignment, which was the only option available for two splines. The
%   analytic curves can be compared exactly instead, because they live in a
%   two-parameter family: they differ only through alpha, fixed by sequencing
%   depth, and phi, the overdispersion. params carries both, so
%
%       mu = logspace(-3, 2, 500)';
%       gap = max(abs(params.curve1.y(mu) - params.curve2.y(mu)));
%
%   gives the separation directly, in the units dd is measured in.
%
%   Prefer that to d. The Procrustes value is NOT comparable with the one
%   SC_SPLINEFIT2 reports: on the bundled data it is 0.0140 here against
%   0.0014 there, not because these curves are ten times further apart but
%   because a spline stops at the last observed gene while the analytic curve
%   carries on to the end of the mean range, so the two point clouds being
%   aligned have different extents. Both numbers are internally consistent;
%   neither can be read across.
%
%   See also SC_ANALYTICFIT, SC_SPLINEFIT2.

arguments
    X {mustBeNumeric, mustBeNonempty}
    Y {mustBeNumeric, mustBeNonempty}
    genelistx (:,1) string
    genelisty (:,1) string = genelistx
    sortid (1,1) logical = true
    opts.LibSizeX (:,1) double = []
    opts.LibSizeY (:,1) double = []
    opts.ScaleFactor (1,1) double {mustBePositive} = 1e4
    opts.IsNormalized (1,1) logical = false
    opts.Dispersion double {mustBeScalarOrEmpty, mustBeNonnegative} = []
    opts.OneSided (1,1) logical = true
end

[genelist, i, j] = intersect(genelistx, genelisty, 'stable');
if isempty(genelist)
    error("sc_analyticfit2:NoCommonGenes", "The two gene lists do not overlap.");
end
X = X(i, :);
Y = Y(j, :);

shared = {'ScaleFactor', opts.ScaleFactor, 'IsNormalized', opts.IsNormalized, ...
    'Dispersion', opts.Dispersion, 'OneSided', opts.OneSided, 'SortIt', false};

[T1, sx, curve1, p1] = sc_analyticfit(X, genelist, 'LibSize', opts.LibSizeX, shared{:});
T1.Properties.VariableNames = {'genes', 'logu1', 'logcv1', ...
    'dropr1', 'd1', 'pval1', 'fdr1', 'nearidx1'};

[T2, sy, curve2, p2] = sc_analyticfit(Y, genelist, 'LibSize', opts.LibSizeY, shared{:});
T2.Properties.VariableNames = {'genes', 'logu2', 'logcv2', ...
    'dropr2', 'd2', 'pval2', 'fdr2', 'nearidx2'};

T = join(T1, T2, 'Keys', 'genes');
T.dd = T.d2 - T.d1;
if sortid
    T = sortrows(T, 'dd', 'descend');
end

% Procrustes, as in SC_SPLINEFIT2. It needs the same number of points on
% each side, which fails if a gene was all-zero in one condition only.
if isequal(size(sx), size(sy))
    [d, sz] = procrustes(sx, sy);
else
    d = NaN;
    sz = [];
    warning("sc_analyticfit2:SizeMismatch", ...
        ["The two conditions kept different numbers of genes (%d and %d), " ...
        "so the Procrustes comparison was skipped. Use params instead."], ...
        size(sx, 1), size(sy, 1));
end

params = struct('alpha1', p1.alpha, 'phi1', p1.phi, ...
    'alpha2', p2.alpha, 'phi2', p2.phi, ...
    'curve1', curve1, 'curve2', curve2);
end
