% function test_gui
cL1 = "Sample1";
cL2 = "Sample2";

lcolors = lines(2);
lcolor1 = lcolors(1,:);
lcolor2 = lcolors(2,:);


outfile = sprintf('%s_vs_%s', ...
        matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));

% ---------
%answer = questdlg('Show gene expression profiles?','');
%if ~strcmp(answer,'Yes'), return; end

pause(1);
hFig = figure('Visible','off');
hFig.Position(3) = hFig.Position(3)*1.8;
% gui.i_movegui2parent(hFig, FigureHandle);

delete(findall(hFig, 'Tag', 'FigureToolBar'))
tb = uitoolbar('Parent', hFig);
% tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
% uipushtool(tb, 'Separator', 'off');


pkg.i_addbutton2fig(tb, 'off', {@in_HighlightSelectedGenes, 1}, 'list.gif', 'Selet a gene to show expression profile');
% pkg.i_addbutton2fig(tb, 'off', @in_HighlightSelectedGenes, 'xplotpicker-qqplot.gif', 'Highlight selected genes');
pkg.i_addbutton2fig(tb, 'off', {@in_HighlightSelectedGenes, 2}, 'list2.gif', 'Selet a gene from sorted list');

pkg.i_addbutton2fig(tb, 'on', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Select top n genes to perform web-based enrichment analysis...');
pkg.i_addbutton2fig(tb, 'off',  @i_genecards, 'fvtool_fdalinkbutton.gif', 'GeneCards...');
pkg.i_addbutton2fig(tb, 'off', @ExportTable, 'export.gif', 'Export HVG Table...');

pkg.i_addbutton2fig(tb, 'on', @ChangeAlphaValue, 'plotpicker-rose.gif', 'Change MarkerFaceAlpha value');
pkg.i_addbutton2fig(tb, 'off', @in_changeMarkerSize, 'icon-mat-text-fields-10.gif', 'ChangeFontSize');
gui.gui_3dcamera(tb, 'DV_Results');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_resizewin, hFig}, 'HDF_pointx.gif', 'Resize Plot Window');


hAx0 = subplot(2,2,[1 3]);
h1 = scatter3(hAx0, px1, py1, pz1, 'filled', 'MarkerFaceAlpha', .1);
hold on
h2 = scatter3(hAx0, px2, py2, pz2, 'filled', 'MarkerFaceAlpha', .1);
plot3(hAx0, xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4, 'Color',lcolor1);
plot3(hAx0, xyz2(:, 1), xyz2(:, 2), xyz2(:, 3), '-', 'linewidth', 4, 'Color',lcolor2);

xlabel(hAx0,'Mean+1, log');
ylabel(hAx0,'CV+1, log');
zlabel(hAx0,'Dropout rate (% of zeros)');


if ~isempty(g)
    dt = datacursormode(hFig);
    datacursormode(hFig, 'on');
    dt.UpdateFcn = {@in_myupdatefcn3, g};
end


idx = find(g==table2array(T(1,1)));


hAx1 = subplot(2,2,2);
x1 = X1(idx,:);
sh1 = plot(hAx1, 1:length(x1), x1, 'Color',lcolor1);
xlim(hAx1,[1 size(X1,2)]);

title(hAx1, strrep(sprintf('%s',g(idx)),'_','\_') );
[titxt] = gui.i_getsubtitle(x1, cL1{1});
subtitle(hAx1, titxt);
xlabel(hAx1,'Cell Index');
ylabel(hAx1,'Expression Level');

hAx2 = subplot(2,2,4);
x2 = X2(idx,:);
sh2 = plot(hAx2, 1:length(x2), x2, 'Color',lcolor2);
xlim(hAx2,[1 size(X2,2)]);

title(hAx2, strrep(sprintf('%s',g(idx)),'_','\_'));
[titxt] = gui.i_getsubtitle(x2, cL2{1});
subtitle(hAx2, titxt);
xlabel(hAx2,'Cell Index');
ylabel(hAx2,'Expression Level');
h3 = []; h4 = []; h5 = [];
h3a = []; h3b = [];
yl = cell2mat(get([hAx1, hAx2], 'Ylim'));
ylnew = [min(yl(:, 1)), max(yl(:, 2))];
set([hAx1, hAx2], 'Ylim', ylnew);

hFig.Visible=true;




   function ExportTable(~, ~)                
       % gui.i_exporttable(T, true, 'Tsplinefitg', 'SplinefitGTable');
        [~, filesaved] = gui.i_exporttable(T, true, 'Tdvgenelist', outfile);
        fprintf('Result has been saved in %s\n',filesaved);        
        % if ~isempty(filesaved)
        %    waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
        % end        
    end

    function txt = in_myupdatefcn3(src, event_obj, g)
        if isequal(get(src, 'Parent'), hAx0)
            subplot(hAx0);
            % dtp = findobj(h1, 'Type', 'datatip');
            % if ~isempty(dtp), delete(dtp); end
            % dtp = findobj(h2, 'Type', 'datatip');
            % if ~isempty(dtp), delete(dtp); end            
            idx = event_obj.DataIndex;
            if idx > length(g)*2
                txt = num2str(event_obj.Position(2)); 
                return;
            end
            if idx > length(g)
                idx = idx - length(g);
            end
            
            x_cleanfigspace(false);

            % if ~isempty(h3), delete(h3); end
            % if ~isempty(h3a), delete(h3a); end
            % if ~isempty(h3b), delete(h3b); end
            % if ~isempty(h4), delete(h4); end
            % if ~isempty(h5), delete(h5); end

                h3 = plot3(hAx0, [px1(idx) px2(idx)], ...
                    [py1(idx), py2(idx)], ...
                    [pz1(idx), pz2(idx)],'k-','LineWidth',1);

                h3a = plot3(hAx0, px1(idx), py1(idx), pz1(idx), 'b.','MarkerSize',10);
                h3b = plot3(hAx0, px2(idx), py2(idx), pz2(idx), 'r.','MarkerSize',10);
                % daspect(hAx0,'auto');
                %h3b = arrow3([px1(idx), py1(idx), pz1(idx)], ...
                %    [px2(idx), py2(idx), pz2(idx)]);
                %h3b = quiver3(hAx0, px1(idx), py1(idx), pz1(idx), px2(idx), py2(idx), pz2(idx));
            txt = {g(idx)};
            x1 = X1(idx, :);
            if ~isempty(sh1) && isvalid(sh1), delete(sh1); end
            sh1 = plot(hAx1, 1:length(x1), x1, 'Color',lcolor1);
            xlim(hAx1,[1 size(X1,2)]);
            title(hAx1, strrep(sprintf('%s',g(idx)),'_','\_') );
            [titxt] = gui.i_getsubtitle(x1, cL1{1});
            subtitle(hAx1, titxt);
            xlabel(hAx1,'Cell Index');
            ylabel(hAx1,'Expression Level');

            x2 = X2(idx, :);
            if ~isempty(sh2) && isvalid(sh2), delete(sh2); end
            sh1 = plot(hAx2, 1:length(x2), x2, 'Color',lcolor2);
            xlim(hAx2,[1 size(X2,2)]);
            title(hAx2, strrep(sprintf('%s',g(idx)),'_','\_'));
            [titxt] = gui.i_getsubtitle(x2, cL2{1});
            subtitle(hAx2, titxt);
            xlabel(hAx2,'Cell Index');
            ylabel(hAx2,'Expression Level');

            yl = cell2mat(get([hAx1, hAx2], 'Ylim'));
            ylnew = [min(yl(:, 1)), max(yl(:, 2))];
            set([hAx1, hAx2], 'Ylim', ylnew);
        else
            txt = num2str(event_obj.Position(2));
        end
    end

    function x_cleanfigspace(deldatatip)
        if nargin<1, deldatatip = false; end
        if deldatatip
           dtp = findobj(h1, 'Type', 'datatip');
           if ~isempty(dtp), delete(dtp); end
           dtp = findobj(h2, 'Type', 'datatip');
           if ~isempty(dtp), delete(dtp); end
        end
        
       if ~isempty(h3), delete(h3); end
       if ~isempty(h3a), delete(h3a); end
       if ~isempty(h3b), delete(h3b); end
       if ~isempty(h4), delete(h4); end
       if ~isempty(h5), delete(h5); end
    end

    function in_HighlightSelectedGenes(~, ~, typeid)
       if nargin < 3, typeid = 1; end

       x_cleanfigspace(true);

       % dtp = findobj(h1, 'Type', 'datatip');
       % if ~isempty(dtp), delete(dtp); end
       % dtp = findobj(h2, 'Type', 'datatip');
       % if ~isempty(dtp), delete(dtp); end
       % if ~isempty(h3), delete(h3); end
       % if ~isempty(h3a), delete(h3a); end
       % if ~isempty(h3b), delete(h3b); end
       % if ~isempty(h4), delete(h4); end
       % if ~isempty(h5), delete(h5); end
       
       switch typeid
           case 1
                gsorted = natsort(g);
           case 2
                gsorted = T.(T.Properties.VariableNames{1});
       end

        [indx2, tf2] = listdlg('PromptString', ...
            'Select a gene:', ...
            'SelectionMode', 'single', ...
            'ListString', gsorted, ...
            'ListSize', [220, 300]);
        if tf2 == 1
            idx = find(g == gsorted{indx2});
        else
            return;
        end
        h1.BrushData = idx;
        h2.BrushData = idx;
        dt1 = datatip(h1, 'DataIndex', idx);
        dt2 = datatip(h2, 'DataIndex', idx);
        
        x_cleanfigspace(false);
        %if ~isempty(h3), delete(h3); end
        %if ~isempty(h4), delete(h4); end
        %if ~isempty(h5), delete(h5); end
         
         % h3 = plot3(hAx0, [px1(idx) px2(idx)], ...
         %     [py1(idx), py2(idx)], ...
         %     [pz1(idx), pz2(idx)],'k--','LineWidth',1);   

         [nearidx] = dsearchn(xyz1, [px1(idx) py1(idx) pz1(idx)]);
         
         % assignin('base',"xyz1",xyz1);
         % h4 = arrow3(xyz1(nearidx, :), [px1(idx), py1(idx), pz1(idx)]);

        h4 = plot3(hAx0, [px1(idx) xyz1(nearidx, 1)], ...
            [py1(idx), xyz1(nearidx, 2)], ...
            [pz1(idx), xyz1(nearidx, 3)],'-','LineWidth',2,'Color',lcolor1);
         
         [nearidx] = dsearchn(xyz2, [px2(idx) py2(idx) pz2(idx)]);

         % h5 = arrow3(xyz2(nearidx, :), [px2(idx), py2(idx), pz2(idx)]);
         
         h5 = plot3(hAx0, [px2(idx) xyz2(nearidx, 1)], ...
             [py2(idx), xyz2(nearidx, 2)], ...
             [pz2(idx), xyz2(nearidx, 3)],'k-','LineWidth',2,'Color',lcolor2);   

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
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', g(idx)),'-new');
    end

    function in_changeMarkerSize(~, ~)
        % h1.Marker
        if h1.SizeData > 40
            h1.SizeData = 10;
            h2.SizeData = 10;
        else
            h1.SizeData = h1.SizeData + 2; 
            h2.SizeData = h2.SizeData + 2;
        end
    end


    function ChangeAlphaValue(~, ~)
        if h1.MarkerFaceAlpha <= 0.05
            h1.MarkerFaceAlpha = 1;
            h2.MarkerFaceAlpha = 1;
 
        else
            h1.MarkerFaceAlpha = h1.MarkerFaceAlpha - 0.1;
            h2.MarkerFaceAlpha = h2.MarkerFaceAlpha - 0.1;
 
        end
    end
% end