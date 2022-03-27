function h = plotTrackStats(tracks, labels, metadata, nbin, figtitle, varargin)

% fast = 0;
prow = 3;
fast = false;

if ismember('fast',varargin)
    fast = true;
    prow = 2;
end
%
% end

if nargin<5 || ~ischar(figtitle)
    figtitle = 'Tracks statistics';
end
if nargin<4 || isempty(nbin) || nbin<1
    nbin = 20;
end
if nargin<3 || numel(metadata)~=numel(tracks)
    metadata = cell(numel(tracks),1);
    for i = 1:numel(tracks)
        metadata{i} = tracks(i).metadata.name_;
    end
end

if nargin<2
    labels = unique(metadata);
end

if ischar(labels)
    labels = {labels};
end

if isnumeric(metadata) || islogical(metadata)
    if isrow(metadata)
        metadata = cellstr(num2str(metadata'));
    else
        metadata = cellstr(num2str(metadata));
    end
end

if isnumeric(labels) || islogical(labels)
    if isrow(labels)
        labels = cellstr(num2str(labels'));
    else
        labels = cellstr(num2str(labels));
    end
end

if size(labels,2)>1
    labels = labels';
end
%%
ind = zeros(numel(metadata),1);
for i = 1:numel(labels)
    ind = ind | contains(metadata,labels(i),'IgnoreCase',true);
end
tracks = tracks(ind);
metadata = metadata(ind);

if isempty([tracks.tumblebias])
    prow = prow-1;
end
%%
visfig = 'On';
if ismember('hpcc',varargin)
    visfig = 'Off';
end
screensize = get( groot, 'Screensize' );
h = figure('Name',figtitle,'Position',[round((screensize(3)-1200)/2)...
    round((screensize(4)-prow*300)/2) 1200 prow*300],'visible',visfig);
c = lines(numel(labels));
l = cell(numel(labels),1);
%% trajectory length histogram
subplot(prow,4,1), hold on;
s = ceil(max([tracks.trajtime]))/nbin;
r = 0:s:max([tracks.trajtime])+1;
for i=1:numel(labels)
    ind = contains(metadata,labels(i),'IgnoreCase',true);
    temp = tracks(ind);
    a = zeros(numel(r),1);
    for j = 1:numel(r)-1
        ii = [temp.trajtime]<r(j+1) & [temp.trajtime]>=r(j);
        a(j) = sum([temp(ii).trajtime]);
    end
    a = a/sum(a)/s;
    stairs(r,a,'color',c(i,:),'linewidth',2);
    l{i} = num2str(round(sum([temp.trajtime])/60));
end
xlim([0 ceil(max([tracks.trajtime]))]);
xlabel('Trajectory length (s)','FontSize',14);
ylabel('Probability density','FontSize',14);
axis square;
set(gca,'XMinorTick','on');
legend(strcat(labels,{' ('},l,{'min)'}),'Interpreter','none',...
    'FontSize',12,'Location','best');
%% Diff coeff
subplot(prow,4,2), hold on,

[x,ix] = sort([tracks.diffcoeff_cve_mean]);
w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
xl = floor(log10(min(x(w>0.01)))*2)/2;
xh = ceil(log10(max(x(w<0.99)))*2)/2;
s = (xh-xl)/nbin;
r = xl-s:s:xh+s;

for i=1:numel(labels)
    ind = contains(metadata,labels(i),'IgnoreCase',true);
    temp = tracks(ind);
    a = zeros(numel(r),1);
    for j = 1:numel(r)-1
        ii = log10([temp.diffcoeff_cve_mean])<r(j+1) &...
            log10([temp.diffcoeff_cve_mean])>=r(j);
        a(j) = sum([temp(ii).trajtime]);
    end
    a = a/sum(a)/s;
    stairs(r,a,'color',c(i,:),'linewidth',2);
end
xlim([xl xh]);
ax = gca;
ax.XTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
ax.XTickLabel = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'...
    '10^3' '10^4' '10^5' '10^6' '10^7' '10^8'};
xlabel('Diffusion coefficient (\mum^2/s)','FontSize',14);
ylabel('Probability density','FontSize',14);
axis square;
%% Swimming speed
subplot(prow,4,3), hold on,

[x,ix] = sort([tracks.meanspeed]);
w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);

maxs = ceil(max(x(w<0.98)));
s = maxs/nbin;
r = 0:s:maxs+s;

for i=1:numel(labels)
    ind = contains(metadata,labels(i),'IgnoreCase',true);
    temp = tracks(ind);
    if ~isempty([temp.meanspeed])
        a = zeros(numel(r),1);
        for j = 1:numel(r)-1
            ii = [temp.meanspeed]<r(j+1) & [temp.meanspeed]>=r(j);
            a(j) = sum([temp(ii).trajtime]);
        end
        a = a/sum(a)/s;
        stairs(r,a,'color',c(i,:),'linewidth',2);
    end
end
xlim([0 maxs]);
xlabel('Swimming speed (\mum/s)','FontSize',14);
ylabel('Probability density','FontSize',14);
axis square;
set(gca,'XMinorTick','on');

%% Variance Swimming speed
subplot(prow,4,4), hold on,
% [x,ix] = sort(sqrt([tracks.varspeed])./[tracks.meanspeed]);
% w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
% xh = ceil(max(x(w<=0.99))*20)/20;
% xl = floor(min(x(w>=0.01))*20)/20;
% s = (xh-xl)/nbin;
% r = xl-s:s:xh+s;
%
% for i=1:numel(labels)
%     ind = contains(metadata,labels(i),'IgnoreCase',true);
%     temp = tracks(ind);
%     if ~isempty([temp.varspeed])
%         a = zeros(numel(r),1);
%         for j = 1:numel(r)-1
%             ii = sqrt([temp.varspeed])./[temp.meanspeed]<r(j+1) & sqrt([temp.varspeed])./[temp.meanspeed]>=r(j);
%             a(j) = sum([temp(ii).trajtime]);
%         end
%         a = a/sum(a)/s;
%         stairs(r,a,'color',c(i,:),'linewidth',2);
%     end
% end
% xlim([xl xh]);
% % ax = gca;
% % ax.XTick = [-3 -2 -1 0 1 2 3 4 5];
% % ax.XTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'...
% %     '10^3' '10^4' '10^5'};
% xlabel('CV speed','FontSize',14);
% ylabel('Probability density','FontSize',14);
% axis square;
% set(gca,'XMinorTick','on');

[x, ix] = sort(log10([tracks.diffcoeff_cve_runtime]));
w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
xh = ceil(max(x(w<0.99))*5)/5;
xl = floor(min(x(w>0.01))*5)/5;
s = (xh-xl)/nbin;
r = xl-s:s:xh+s;
for i=1:numel(labels)
    ind = contains(metadata,labels(i),'IgnoreCase',true);
    temp = tracks(ind);
    if ~isempty([temp.diffcoeff_cve_runtime])
        a = zeros(numel(r),1);
        for j = 1:numel(r)-1
            ii = log10([temp.diffcoeff_cve_runtime])<r(j+1) & log10([temp.diffcoeff_cve_runtime])>=r(j);
            a(j) = sum([temp(ii).trajtime]);
        end
        a = a/sum(a)/s;
        stairs(r,a,'color',c(i,:),'linewidth',2);
    end
end
xlim([xl xh]);
ax = gca;
ax.XTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
ax.XTickLabel = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'...
    '10^3' '10^4' '10^5'};
xlabel('Dir. persistence time (s)','FontSize',14);
ylabel('Probability density','FontSize',14);
axis square;

if ~fast
    
    %% MSD
    subplot(prow,4,4*(prow-1)+2), hold on,
    msd_all = cell(numel(labels),1);
    for i=1:numel(labels)
        ind = contains(metadata,labels(i),'IgnoreCase',true);
        temp = tracks(ind);
        ut = unique(round(cell2mat(cellfun(@(x) (x - x(1)),{temp.time}','UniformOutput',0)),5));
        meanmsd = zeros(numel(ut),1);
        w = zeros(numel(ut),1);
        
        for j = 1:numel(temp)
            for k = 2:(numel(temp(j).x))
                x_sq_disp = (temp(j).x(k:end) - temp(j).x(1:end-k+1)).^2;
                y_sq_disp = (temp(j).y(k:end) - temp(j).y(1:end-k+1)).^2;
                if sum(~isnan(x_sq_disp))>0
                    meanmsd(k) = ((meanmsd(k) * w(k)) + nansum(x_sq_disp + y_sq_disp))/(w(k)+sum(~isnan(x_sq_disp)));
                    w(k) = w(k)+sum(~isnan(x_sq_disp));
                end
            end
        end
        
        if numel(w)>1
            scatter(flipud(log10(ut)),flipud(log10(meanmsd)),[],flipud(1-(1-c(i,:)).*(0.05+(w/w(2)))/1.05),'Marker','o','MarkerFacecolor','flat',...
                'MarkeredgeColor','none');
        end
        msd_all{i} = meanmsd;
    end
    ax = gca;
    ax.XTick = [-3 -2 -1 0 1 2 3 4];
    ax.XTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2' '10^3'...
        '10^4'};
    ax.YTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
    ax.YTickLabel = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2' '10^3'...
        '10^4' '10^5' '10^6' '10^7' '10^8'};
    xlabel('Time delay (s)','FontSize',14);
    ylabel('MSD (\mum^2)','FontSize',14);
    axis square;
    yl = ylim;
    plot([-2 2], [yl(1) yl(1)+4],':','color','k');
    plot([-2 2], [yl(2)-2*4 yl(2)],':','color','k');
    xlim([-2 2]);
    ylim(yl);
    %% VACF
    
    subplot(prow,4,4*(prow-1)+1), hold on,
    for i = 1:numel(labels)
        meanmsd = msd_all{i};
        ind = contains(metadata,labels(i),'IgnoreCase',true);
        if sum(ind)>0
            temp = tracks(ind);
            t = unique(round(cell2mat(cellfun(@(x) (x - x(1)),{temp.time}','UniformOutput',0)),5));
            vacf = zeros(numel(t),1);
            n = zeros(numel(t),1);
            for k = 1:numel(temp)
                for j = 1:numel(temp(k).time)
                    d = dot([temp(k).xvec(1:(end-j+1)), ...
                        temp(k).yvec(1:(end-j+1))], ...
                        [temp(k).xvec(j:(end)), ...
                        temp(k).yvec(j:(end))],2);
                    if sum(~isnan(d))>0
                        vacf(j) = (vacf(j)*n(j) + nansum(d))/(n(j) + sum(~isnan(d)));
                        n(j) = n(j) + sum(~isnan(d));
                    end
                end
            end
            
            vacfmsd = diff(diff(meanmsd))./diff(t(1:end-1)).^2./2;
            vacf = nanmean(horzcat(vacf, [NaN; vacfmsd; NaN]),2);
            
            scatter(flipud(log10(t)),flipud((vacf/vacf(1))),[],flipud(1-(1-c(i,:)).*(0.05+(n/n(1)))/1.05),'Marker','o','MarkerFacecolor','flat',...
                'MarkeredgeColor','none');
            f = fit(t,(vacf/vacf(1)),'exp(-x/a)*cos(x*w)',...
                'Weight',n, 'Startpoint',[1 1],'Lower',[0 0],'Upper',[Inf Inf]);
            fout = feval(f,10.^(-3:0.04:2));
            plot(-3:0.04:2,fout,':','linewidth',1,'color',c(i,:));
            
        end
    end
    
    plot([-3 2], [0 0], 'k:');
    ax = gca;
    ax.XTick = [-3 -2 -1 0 1 2];
    ax.XTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'};
    xlabel('Time delay (s)','FontSize',14);
    ylabel('Dir. autocorrelation','FontSize',14);
    axis square;
    set(gca,'XMinorTick','on');
    xlim([-3 1]);
    %     yl = ylim;
    %     ylim([yl(1) 1]);
    ylim([-0.5 1]);
    
    
    %% alpha
    subplot(prow,4,4*(prow-1)+3), hold on,
    for i=1:numel(labels)
        ind = contains(metadata,labels(i),'IgnoreCase',true);
        if sum(ind)>0
            temp = tracks(ind);
            t = round(cell2mat(cellfun(@(x) (x - x(1)),{temp.time}','UniformOutput',0)),5);
            f = cell2mat(cellfun(@(x) (x - x(1)),{temp.frozen}','UniformOutput',0));
            t = t(~f);
            ut = unique(t);
            if numel(ut)>3
                alpha = vertcat(temp.alphamsd);
                n = vertcat(temp.diffcoeff_cve_n);
                meanalpha = zeros(numel(ut),1);
                w = zeros(numel(ut),1);
                for j = 1:numel(ut)
                    meanalpha(j) = nansum(n(t==ut(j)).*alpha(t==ut(j))) / nansum(n(t==ut(j)));
                    w(j) = nansum(n(t==ut(j)));
                end
                
                meanalpha(end) = NaN;
                
                %                 [~, xp, yp] = slmengine(-log10(ut),meanalpha,...
                %                     'plot','off','weights',w/sum(w),'endconditions','notaknot',...
                %                     'interiorknots','fixed','knots',-5:4/nbin:5,...
                %                     'order',1,'predictions',101);
                %                 plot(xp, yp,'-','linewidth',2,'color',c(i,:));
            end
        end
    end
    plot([-3 2],[1 1],':','color','k');
    ax = gca;
    ax.XTick = [-3 -2 -1 0 1 2];
    ax.XTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'};
    ylim([0 2]);
    xlim([-2 2]);
    xlabel('Frequency (s^{-1})','FontSize',14);
    ylabel('Anomalous diff. exp.','FontSize',14);
    axis square;
    %% object intensity
    subplot(prow,4,4*(prow-1)+4), hold on,
    objintensity = cellfun(@(x) log10(nanmean(abs(x))), {tracks.objintensity});
    [x,ix] = sort(objintensity);
    w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
    minint = floor(min(x(w>0.01))*5)/5;
    maxint = ceil(max(x(w<0.99))*5)/5;
    s = (maxint-minint)/nbin;
    r = minint-s:s:maxint+s;
    
    for i=1:numel(labels)
        ind = contains(metadata,labels(i),'IgnoreCase',true);
        temp = tracks(ind);
        
        objintensity = cellfun(@(x) log10(nanmean(abs(x))), {temp.objintensity});
        a = zeros(numel(r),1);
        for j = 1:numel(r)-1
            ii = objintensity<r(j+1) & objintensity>=r(j);
            a(j) = sum([temp(ii).trajtime]);
        end
        a = a/sum(a)/s;
        stairs(r,a,'color',c(i,:),'linewidth',2);
    end
    xlim([minint maxint]);
    xlabel('Object intensity (AU)','FontSize',14);
    ylabel('Probability density','FontSize',14);
    axis square;
    
    set(gca,'XMinorTick','on');
end

indDiff = [tracks.diffcoeff_cve_mean]>=10;
metadata=metadata(indDiff);
tracks=tracks(indDiff);

if ~isempty([tracks.tumblebias])
    
    %% Tumble bias
    subplot(prow,4,5), hold on,
    s = 1/nbin;
    r = 0:s:1;
    for i=1:numel(labels)
        ind = contains(metadata,labels(i),'IgnoreCase',true);
        temp = tracks(ind);
        a = zeros(numel(r),1);
        for j = 1:numel(r)-1
            ii = [temp.tumblebias]<r(j+1) & [temp.tumblebias]>=r(j);
            a(j) = sum([temp(ii).trajtime]);
        end
        a = a/sum(a)/s;
        stairs(r,a,'color',c(i,:),'linewidth',2);
    end
    xlim([0 1]);
    xlabel('Tumble bias','FontSize',14);
    ylabel('Probability density','FontSize',14);
    axis square;
    set(gca,'XMinorTick','on');
    
    %% Tumble angle
    subplot(prow,4,6), hold on,
    s = 180/nbin;
    r = 0:s:180;
    for i=1:numel(labels)
        ind = contains(metadata,labels(i),'IgnoreCase',true);
        temp = tracks(ind);
        tumbleangles = abs(cat(2,temp.tumbleangle));
        if ~isempty(tumbleangles)
            [a,~] = histc(tumbleangles,r);
            a = a/sum(a)/s;
            stairs(r,a,'color',c(i,:),'linewidth',2);
        end
    end
    xlim([0 180]);
    xlabel('Tumble angle','FontSize',14);
    ylabel('Probability density','FontSize',14);
    axis square;
    set(gca,'XMinorTick','on');
    
    %% Tumble time
    subplot(prow,4,7), hold on,
    s = 3/nbin;
    r = -2:s:1;
    for i=1:numel(labels)
        ind = contains(metadata,labels(i),'IgnoreCase',true);
        temp = tracks(ind);
        if ~isempty([temp.tumbletime])
            [a,~] = histc(log10([temp.tumbletime]),r);
            a = a/sum(a)/s;
            stairs(r,a,'color',c(i,:),'linewidth',2);
        end
    end
    xlim([-2 1]);
    xlabel('Tumble time (s)','FontSize',14);
    ylabel('Probability density','FontSize',14);
    ax = gca;
    ax.XTick = [-2 -1 0 1];
    ax.XTickLabel = {'10^{-2}' '10^{-1}' '10^0' '10^1'};
    axis square;
    %% Reversion Freq
    subplot(prow,4,8), hold on,
    xh = ceil(quantile([tracks.revfreq],0.98));
    if ~isnan(xh)
        s = xh/nbin;
        r = 0:s:(xh+1);
        for i=1:numel(labels)
            ind = contains(metadata,labels(i),'IgnoreCase',true);
            if sum(ind)>0
                temp = tracks(ind);
                if ~isempty([temp.revfreq])
                    a = zeros(numel(r),1);
                    for j = 1:numel(r)-1
                        ii = [temp.revfreq]<r(j+1) & [temp.revfreq]>=r(j);
                        a(j) = sum([temp(ii).trajtime]);
                    end
                    a = a/sum(a)/s;
                    stairs(r,a,'color',c(i,:),'linewidth',2);
                end
            end
        end
        xlim([0 xh]);
    end
    xlabel('Reversion freq. (s^{-1})','FontSize',14);
    ylabel('Probability density','FontSize',14);
    axis square;
    set(gca,'XMinorTick','on');
end
%%
if ismember('svg',varargin)
    set(h,'PaperPositionMode','auto');
    print(h,strcat(figtitle,sprintf('.%s','svg')),sprintf('-d%s','svg'),'-r150','-noui');
end
if ismember('png',varargin)
    set(h,'PaperPositionMode','auto');
    print(h,strcat(figtitle,sprintf('.%s','png')),sprintf('-d%s','png'),'-r150','-noui');
    crop(strcat(figtitle,sprintf('.%s','png')));
end

end