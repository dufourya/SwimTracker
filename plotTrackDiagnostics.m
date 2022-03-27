function plotTrackDiagnostics(fileName,hpcc)

if nargin<2 || isempty(hpcc)
    visfig = 'On';
elseif hpcc > 0
    visfig = 'Off';
end

if nargin<1 || isempty(fileName)
    track_files = dir('*.rad_swimtracker.mat');
else
    fileName = strsplit(fileName,'.');
    track_files = dir(strcat(fileName{1},'.rad_swimtracker.mat'));
end
%%
screensize = get( groot, 'Screensize' );

for k = 1:numel(track_files)
    %%
    fname = track_files(k).name;
    sname = strsplit(fname,'.');
    sname = sname{1};
    fprintf('Plotting %s...',sname);
    if exist(strcat(sname,'.diagnostics.png'),'file') % && exist(strcat(sname,'.diagnostics.svg'),'file')
        fprintf('already done.\n');
        continue;
    end
    load(fname);
    %%
    if numel(tracks)>1
        
        h = figure('color','w','Name',sname,'Position',[round((screensize(3)-900)/2)...
            round((screensize(4)-300)/2) 900 300],'visible',visfig);
        % detectobj = cellfun(@(x) size(x,1), {rad_detection.xCoord});
        % detecttime = [rad_detection.time];
        objtime = vertcat(tracks.time);
        objtime = floor(objtime/ceil((max(objtime)-min(objtime))/10))*ceil((max(objtime)-min(objtime))/10);
        t = unique(objtime);
        tedges = [t; t(end)+(t(2)-t(1))];
        % detectobj = arrayfun(@(x) mean(detectobj(detecttime>=tedges(x) & detecttime<tedges(x+1))),1:numel(tedges)-1);
        N = histcounts(objtime,tedges);
        Nt = histcounts(unique(vertcat(tracks.time)),tedges);
        t = t+(t(2)-t(1))/2;
        objtime = objtime+(t(2)-t(1))/2;
        g1 = subplot(1,3,1);
        % hold on,
        scatter(t,N./Nt,'filled');
        % plot(t,detectobj,'*');
        ax = gca;
        xlim([0 max(tedges)]);
        ylim([0 ceil(ax.YLim(2)/10^floor(log10(ax.YLim(2))))*10^floor(log10(ax.YLim(2)))]);
        axis square;
        xlabel('Time (s)');
        ylabel('Number of objects');
        box on;
        g1.Position = [0.1 0.1 0.2 0.9];
        g1.XTick = t;
        g2 = subplot(1,3,2);
        data = vertcat(tracks.objintensity);
        boxplot(data,objtime,'PlotStyle','compact','Positions',t,'Symbol','','Jitter',0.1,'Color',lines(1),'LabelOrientation','horizontal')
        ylim([min(0,quantile(data,0.025)) quantile(data,0.985)]);
        xlabel('Time (s)');
        ylabel('Object intensities (AU)');
        axis square;
        g2.Position = [0.4 0.1 0.2 0.9];
        g3 = subplot(1,3,3);
        data = vertcat(tracks.speed);
        boxplot(data,objtime,'PlotStyle','compact','Positions',t,'Symbol','','Jitter',0.1,'Color',lines(1),'LabelOrientation','horizontal')
        ylim([min(0,quantile(data,0.025)) quantile(data,0.985)]);
        xlabel('Time (s)');
        ylabel('Object speed (\mum/s)');
        axis square;
        g3.Position = [0.7 0.1 0.2 0.9];
        set(h,'PaperPositionMode','auto')
        print(h,strcat(tracks(1).metadata.name_,'.diagnostics.png'),'-dpng','-r150','-noui');
        crop(strcat(tracks(1).metadata.name_,'.diagnostics.png'));
%         print(h,strcat(tracks(1).metadata.name_,'.diagnostics.svg'),'-dsvg','-r150','-noui');
        fprintf('done.\n');
        
    end
    %%
end
end