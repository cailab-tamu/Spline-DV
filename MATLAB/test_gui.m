function test_gui(T, g, X1, X2, P1, P2, xyz1, xyz2)
%TEST_GUI  Interactive view of the Spline-DV result produced by TEST.
%
%   Run test.m first, then
%
%       test_gui(T, g, X1, X2, [px1 py1 pz1], [px2 py2 pz2], xyz1, xyz2)
%
%   The left panel shows the two gene clouds with their reference curves.
%   Click a gene to see its expression profile across the cells of each
%   sample in the two panels on the right.
%
%   INPUTS
%     T         result table, first column the gene names, sorted by DV score
%     g         gene names, in the same order as X1, X2 and P1, P2
%     X1, X2    normalized expression, genes by cells, one per sample
%     P1, P2    gene coordinates [log mean, log CV, dropout], one per sample
%     xyz1,xyz2 points along each sample's reference curve, for plotting
%
%   NOTE: this must stay a function, not a script. The callbacks below are
%   nested functions, and that is what lets them see T, g, X1 and the axes
%   handles. Local functions in a script get their own empty workspace, so
%   turning this back into a script would break every button silently.

cL1 = "Sample1";
cL2 = "Sample2";

lcolors = lines(2);
lcolor1 = lcolors(1, :);
lcolor2 = lcolors(2, :);

px1 = P1(:, 1);  py1 = P1(:, 2);  pz1 = P1(:, 3);
px2 = P2(:, 1);  py2 = P2(:, 2);  pz2 = P2(:, 3);

outfile = sprintf('%s_vs_%s', ...
    matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));

hFig = figure('Visible', 'off');
hFig.Position(3) = hFig.Position(3) * 1.8;

delete(findall(hFig, 'Tag', 'FigureToolBar'))
tb = uitoolbar('Parent', hFig);

gui.i_addbutton2fig(tb, 'off', {@in_HighlightSelectedGenes, 1}, 'list.gif', 'Select a gene to show expression profile');
gui.i_addbutton2fig(tb, 'off', {@in_HighlightSelectedGenes, 2}, 'list2.gif', 'Select a gene from sorted list');

gui.i_addbutton2fig(tb, 'on', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Select top n genes to perform web-based enrichment analysis...');
gui.i_addbutton2fig(tb, 'off', @i_genecards, 'fvtool_fdalinkbutton.gif', 'GeneCards...');
gui.i_addbutton2fig(tb, 'off', @ExportTable, 'export.gif', 'Export DV Table...');

gui.i_addbutton2fig(tb, 'on', @ChangeAlphaValue, 'plotpicker-rose.gif', 'Change MarkerFaceAlpha value');
gui.i_addbutton2fig(tb, 'off', @in_changeMarkerSize, 'icon-mat-text-fields-10.gif', 'Change marker size');
gui.gui_3dcamera(tb, 'DV_Results');
gui.i_addbutton2fig(tb, 'on', {@gui.i_resizewin, hFig}, 'HDF_pointx.gif', 'Resize Plot Window');


%% Left panel: the gene clouds and the two reference curves

hAx0 = subplot(2, 2, [1 3]);
h1 = scatter3(hAx0, px1, py1, pz1, 'filled', 'MarkerFaceAlpha', .1);
hold on
h2 = scatter3(hAx0, px2, py2, pz2, 'filled', 'MarkerFaceAlpha', .1);
plot3(hAx0, xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4, 'Color', lcolor1);
plot3(hAx0, xyz2(:, 1), xyz2(:, 2), xyz2(:, 3), '-', 'linewidth', 4, 'Color', lcolor2);

xlabel(hAx0, 'Mean+1, log');
ylabel(hAx0, 'CV+1, log');
zlabel(hAx0, 'Dropout rate (% of zeros)');

if ~isempty(g)
    dt = datacursormode(hFig);
    datacursormode(hFig, 'on');
    dt.UpdateFcn = {@in_myupdatefcn3, g};
end


%% Right panels: the expression profile of the top gene

idx = find(g == table2array(T(1, 1)));

hAx1 = subplot(2, 2, 2);
x1 = X1(idx, :);
sh1 = plot(hAx1, 1:length(x1), x1, 'Color', lcolor1);
xlim(hAx1, [1 size(X1, 2)]);
title(hAx1, strrep(sprintf('%s', g(idx)), '_', '\_'));
subtitle(hAx1, gui.i_getsubtitle(x1, cL1{1}));
xlabel(hAx1, 'Cell Index');
ylabel(hAx1, 'Expression Level');

hAx2 = subplot(2, 2, 4);
x2 = X2(idx, :);
sh2 = plot(hAx2, 1:length(x2), x2, 'Color', lcolor2);
xlim(hAx2, [1 size(X2, 2)]);
title(hAx2, strrep(sprintf('%s', g(idx)), '_', '\_'));
subtitle(hAx2, gui.i_getsubtitle(x2, cL2{1}));
xlabel(hAx2, 'Cell Index');
ylabel(hAx2, 'Expression Level');

% Handles of the lines drawn when a gene is picked; the callbacks below
% delete them before drawing new ones.
h3 = [];  h3a = [];  h3b = [];  h4 = [];  h5 = [];

i_matchylim();
hFig.Visible = true;


%% Callbacks

    function ExportTable(~, ~)
        [~, filesaved] = gui.i_exporttable(T, true, 'Tdvgenelist', outfile);
        fprintf('Result has been saved in %s\n', filesaved);
    end

    function txt = in_myupdatefcn3(src, event_obj, g)
        if ~isequal(get(src, 'Parent'), hAx0)
            txt = num2str(event_obj.Position(2));
            return;
        end
        subplot(hAx0);
        idx = event_obj.DataIndex;
        if idx > length(g) * 2
            txt = num2str(event_obj.Position(2));
            return;
        end
        if idx > length(g)                 % the second cloud continues the first
            idx = idx - length(g);
        end

        x_cleanfigspace(false);

        % Join the gene's two positions, so the DV score is visible as a
        % line between them.
        h3  = plot3(hAx0, [px1(idx) px2(idx)], [py1(idx) py2(idx)], ...
            [pz1(idx) pz2(idx)], 'k-', 'LineWidth', 1);
        h3a = plot3(hAx0, px1(idx), py1(idx), pz1(idx), 'b.', 'MarkerSize', 10);
        h3b = plot3(hAx0, px2(idx), py2(idx), pz2(idx), 'r.', 'MarkerSize', 10);

        txt = {g(idx)};
        i_showprofiles(idx);
    end

    function in_HighlightSelectedGenes(~, ~, typeid)
        if nargin < 3, typeid = 1; end
        x_cleanfigspace(true);

        switch typeid
            case 1
                gsorted = sort(g);                              % alphabetical
            case 2
                gsorted = T.(T.Properties.VariableNames{1});    % by DV score
        end

        [indx2, tf2] = listdlg('PromptString', 'Select a gene:', ...
            'SelectionMode', 'single', 'ListString', gsorted, ...
            'ListSize', [220, 300]);
        if tf2 ~= 1, return; end
        idx = find(g == gsorted(indx2));

        brushed = zeros(1, numel(g));   % BrushData wants numeric, not logical
        brushed(idx) = 1;
        h1.BrushData = brushed;
        h2.BrushData = brushed;
        datatip(h1, 'DataIndex', idx);
        datatip(h2, 'DataIndex', idx);

        x_cleanfigspace(false);

        % Draw each gene's deviation from its own reference curve.
        near1 = dsearchn(xyz1, [px1(idx) py1(idx) pz1(idx)]);
        h4 = plot3(hAx0, [px1(idx) xyz1(near1, 1)], [py1(idx) xyz1(near1, 2)], ...
            [pz1(idx) xyz1(near1, 3)], '-', 'LineWidth', 2, 'Color', lcolor1);

        near2 = dsearchn(xyz2, [px2(idx) py2(idx) pz2(idx)]);
        h5 = plot3(hAx0, [px2(idx) xyz2(near2, 1)], [py2(idx) xyz2(near2, 2)], ...
            [pz2(idx) xyz2(near2, 3)], '-', 'LineWidth', 2, 'Color', lcolor2);

        i_showprofiles(idx);
    end

    function EnrichrHVGs(~, ~)
        k = gui.i_inputnumk(200, 1, 2000, 'Select top n genes');
        if ~isempty(k)
            gsorted = T.(T.Properties.VariableNames{1});
            gselected = gsorted(1:k);
            fprintf('%d genes are selected.\n', length(gselected));
            gui.i_enrichtest(gselected, gsorted, k);
        end
    end

    function i_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', g(idx)), '-new');
    end

    function in_changeMarkerSize(~, ~)
        if h1.SizeData > 40
            h1.SizeData = 10;
            h2.SizeData = 10;
        else
            h1.SizeData = h1.SizeData + 2;
            h2.SizeData = h2.SizeData + 2;
        end
    end

    function ChangeAlphaValue(~, ~)
        % Step through a few useful opacities. This used to subtract 0.1 per
        % click and wrap only at 0.05 or below, so from the starting 0.1 the
        % first click landed exactly on 0 and both clouds disappeared, which
        % reads as the button being broken rather than as one step of a cycle.
        levels = [0.05 0.1 0.25 0.5 1];
        [~, k] = min(abs(levels - h1.MarkerFaceAlpha));
        a = levels(mod(k, numel(levels)) + 1);
        h1.MarkerFaceAlpha = a;
        h2.MarkerFaceAlpha = a;
    end


%% Small helpers, shared by the callbacks above

    function i_showprofiles(k)
        % Redraw the two expression profiles for gene k.
        xa = X1(k, :);
        if ~isempty(sh1) && isvalid(sh1), delete(sh1); end
        sh1 = plot(hAx1, 1:length(xa), xa, 'Color', lcolor1);
        xlim(hAx1, [1 size(X1, 2)]);
        title(hAx1, strrep(sprintf('%s', g(k)), '_', '\_'));
        subtitle(hAx1, gui.i_getsubtitle(xa, cL1{1}));
        xlabel(hAx1, 'Cell Index');
        ylabel(hAx1, 'Expression Level');

        xb = X2(k, :);
        if ~isempty(sh2) && isvalid(sh2), delete(sh2); end
        sh2 = plot(hAx2, 1:length(xb), xb, 'Color', lcolor2);
        xlim(hAx2, [1 size(X2, 2)]);
        title(hAx2, strrep(sprintf('%s', g(k)), '_', '\_'));
        subtitle(hAx2, gui.i_getsubtitle(xb, cL2{1}));
        xlabel(hAx2, 'Cell Index');
        ylabel(hAx2, 'Expression Level');

        i_matchylim();
    end

    function i_matchylim()
        % Put both profile panels on the same y scale, so the difference in
        % variability between the samples is visible.
        yl = cell2mat(get([hAx1, hAx2], 'Ylim'));
        set([hAx1, hAx2], 'Ylim', [min(yl(:, 1)), max(yl(:, 2))]);
    end

    function x_cleanfigspace(deldatatip)
        if nargin < 1, deldatatip = false; end
        if deldatatip
            delete(findobj(h1, 'Type', 'datatip'));
            delete(findobj(h2, 'Type', 'datatip'));
        end
        delete(h3(isgraphics(h3)));
        delete(h3a(isgraphics(h3a)));
        delete(h3b(isgraphics(h3b)));
        delete(h4(isgraphics(h4)));
        delete(h5(isgraphics(h5)));
    end

end
