function plotSwimTracker(fileName,hpcc)

if nargin<2 || isempty(hpcc) || hpcc == 0
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
% options = statset('MaxIter',1000);
for k = 1:numel(track_files)
    %%
    fname = track_files(k).name;
    sname = strsplit(fname,'.');
    sname = sname{1};
    fprintf('Plotting %s...',sname);
    if exist(strcat(sname,'.png'),'file') % && exist(strcat(sname,'.svg'),'file')
        fprintf('already done.\n');
        continue;
    end
    load(fname);
    %%
    if numel(tracks)>1
        nbin = 20;
        sc = get( groot, 'Screensize' );
        
        if ~isempty([tracks.tumblebias])
            nrow = 3;
        else
            nrow = 2;
        end
        
        h = figure('color','w','Name',sname,'Position',...
            [round((sc(3)-1800)/2) round((sc(4)-300*nrow)/2) 1800 300*nrow],...
            'visible',visfig);
        %%
        if ~exist('movie_drift','var')
            xVec = mean(cat(1,tracks.xvec));
            yVec = mean(cat(1,tracks.yvec));
            movie_drift.total = sqrt(xVec^2+yVec^2);
        end
        %         ind = true(numel(tracks),1);
        %         ind(isnan([tracks.diffcoeff_cve_mean])) = false;
        %
        %         tracks_rejected = tracks(~ind);
        %         tracks = tracks(ind);
        
        clust = gmClusterTracks(tracks,7,1)';
        tracks(clust==0) = [];
        clust(clust==0) = [];
        uclust = unique(clust);
        diffclust = uclust;
        cumulclust = uclust;
        
        for jj = 1:numel(uclust)
            diffclust(jj) = sum(log10([tracks(clust==uclust(jj)).diffcoeff_cve_mean]).*[tracks(clust==uclust(jj)).trajtime])...
                /sum([tracks(clust==uclust(jj)).trajtime]);
            cumulclust(jj) = sum([tracks(clust==uclust(jj)).trajtime]);
        end
        
        cumulclust = cumsum(cumulclust)/sum(cumulclust);
        cumulclust = [0 cumulclust(1:end-1)];
        
        if isempty(tracks)
            subplot(nrow,6,1),
            text(0,0,'No tracks!');
            axis off;
            continue;
        end
        
        %         mSpeed = sum([tracks.meanspeed].*[tracks.trajtime])/sum([tracks.trajtime]);
        %         mDiff = sum(log10([tracks.diffcoeff_cve_mean]).*[tracks.trajtime])/sum([tracks.trajtime]);
        %% trajectories
        subplot(nrow,6,1),
        hold on,
        
        x = [tracks.trajtime]'.*rand(numel(tracks),1);
        sx = sort(x,'descend');
        ind = find(x>=sx(min(50,numel(sx))));
        c = lines(max(uclust));
        for j = ind'
            plot(tracks(j).x, tracks(j).y,'linewidth',1.5,'color',c(clust(j),:));
        end
        axis equal
        axis off
        xl = xlim;
        yl = ylim;
        xd = xl(2)-xl(1);
        yd = yl(2)-yl(1);
        text(xl(1)+xd/2-0.2*xd,yl(1)-yd/10,strcat(sprintf('Drift = %0.2f',movie_drift.total),' \mum/s'));
        title(sname,'Interpreter', 'none');
        %% histogram time
        subplot(nrow,6,2),
        hold on,
        s = ceil(max([tracks.trajtime]))/nbin;
        r = 0:s:max([tracks.trajtime])+1;
        a = zeros(numel(r),max(uclust));
        %         ar = zeros(numel(r),numel(uclust));
        for j = 1:numel(r)-1
            for jj = 1:max(uclust)
                ii = [tracks.trajtime]<r(j+1) & [tracks.trajtime]>=r(j) & clust == jj;
                a(j,jj) = sum([tracks(ii).trajtime]);
            end
        end
        %         totala = sum([a; ar]);
        an = a/sum(a(:))/s;
        %         ar = ar/totala/s;
        ha = bar(r'+s/2,an,'stacked');
        for jj = 1:max(uclust)
            ha(jj).BarWidth = 1;
            ha(jj).FaceColor = c(jj,:);
            ha(jj).EdgeColor = 'none';
            %         ha(2).FaceColor = c(2,:);
        end
        %         ha(2).EdgeColor = 'none';
        xlim([0 ceil(max([tracks.trajtime]))]);
        xlabel('Trajectory length (s)','FontSize',12);
        ylabel('Probability density','FontSize',12);
        set(gca,'XMinorTick','on');
        %         title(sprintf('Mean: %g s',round(sum(s*(r(2:end)-diff(r)/2)'.*a(1:end-1)),2,'significant')));
        axis square;
        box on;
        legend(compose('%G min',round(sum(a)/60,3,'significant')),'location','Best');
        
        %% Diff coeff
        subplot(nrow,6,3),
        [x,ix] = sort([tracks.diffcoeff_cve_mean]);
        w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
        mindiff = floor(2*log10(min(x(w>=0.005))))/2;
        maxdiff = ceil(2*log10(max(x(w<=0.995))))/2;
        
        s = (maxdiff-mindiff)/nbin;
        %         r = min(log10([tracks.diffcoeff_cve_mean]))-s:s:max(log10([tracks.diffcoeff_cve_mean]))+s;
        r = mindiff-s:s:maxdiff+s;
        a = zeros(numel(r),max(uclust));
        for j = 1:numel(r)-1
            for jj = 1:max(uclust)
                ii = log10([tracks.diffcoeff_cve_mean])<r(j+1) &...
                    log10([tracks.diffcoeff_cve_mean])>=r(j) & ...
                    clust == jj;
                a(j,jj) = sum([tracks(ii).trajtime]);
            end
        end
        a = a/sum(a(:))/s;
        hb = bar(r'+s/2,a,'stacked');
        for jj = 1:max(uclust)
            hb(jj).BarWidth = 1;
            hb(jj).FaceColor = c(jj,:);
            hb(jj).EdgeColor = 'none';
        end
        xlim([mindiff maxdiff]);
        ax = gca;
        ax.XTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
        ax.XTickLabel = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}'...
            '10^0' '10^1' '10^2' '10^3' '10^4' '10^5' '10^6' '10^7' '10^8'};
        xlabel('Diffusion coefficient (\mum^2/s)','FontSize',12);
        ylabel('Probability density','FontSize',12);
        axis square;
        %         title(strcat(sprintf('Mean: %g',round(10^mDiff,2,'significant')),' \mum^2/s'));
        
        %% Directional persistence
        %         subplot(nrow,6,4), hold on,
        %         meanrt = 10^(nansum(log10([tracks.diffcoeff_cve_runtime]).*[tracks.trajtime])./sum([tracks.trajtime]));
        %         [x, ix] = sort(log10([tracks.diffcoeff_cve_runtime]));
        %         w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
        %         xh = ceil(2*max(x(w<0.995)))/2;
        %         xl = floor(2*min(x(w>0.005)))/2;
        %         s = (xh-xl)/nbin;
        %         r = xl-s:s:xh+s;
        %
        %         a = zeros(numel(r),1);
        %         for j = 1:numel(r)-1
        %             ii = log10([tracks.diffcoeff_cve_runtime])<r(j+1) & log10([tracks.diffcoeff_cve_runtime])>=r(j);
        %             a(j) = sum([tracks(ii).trajtime]);
        %         end
        %         a = a/sum(a)/s;
        %         hb = bar(r,a,'histc');
        %         set(hb,'facecolor',c(1,:),'edgecolor','none');
        %
        %         xlim([xl xh]);
        %         ax = gca;
        %         ax.XTick = [-3 -2 -1 0 1 2 3 4 5];
        %         ax.XTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'...
        %             '10^3' '10^4' '10^5'};
        %         xlabel('Dir. persistence time (s)','FontSize',12);
        %         ylabel('Probability density','FontSize',12);
        %         axis square;
        %         box on;
        %         title(strcat(sprintf('Mean: %g',round(meanrt,2,'significant')),' s'));
        
        %% Swimming speed
        subplot(nrow,6,4),
        [x,ix] = sort([tracks.meanspeed]);
        w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
        maxs = ceil(max(x(w<=0.995)));
        s = maxs/nbin;
        r = 0:s:maxs+s;
        a = zeros(numel(r),max(uclust));
        for j = 1:numel(r)-1
            for jj = 1:max(uclust)
                ii = [tracks.meanspeed]<r(j+1) & [tracks.meanspeed]>=r(j) & clust == jj;
                a(j,jj) = sum([tracks(ii).trajtime]);
            end
        end
        a = a/sum(a(:))/s;
        hb = bar(r'+s/2,a,'stacked');
        for jj = 1:max(uclust)
            hb(jj).BarWidth = 1;
            hb(jj).FaceColor = c(jj,:);
            hb(jj).EdgeColor = 'none';
        end
        xlim([0 maxs]);
        xlabel('Mean speed (\mum/s)','FontSize',12);
        ylabel('Probability density','FontSize',12);
        axis square;
        set(gca,'XMinorTick','on');
        %         title(strcat(sprintf('Mean: %g',round(mSpeed,2,'significant')),' \mum/s'));
        
        %% Swimming acceleration
        subplot(nrow,6,5),
        %         mAcc = sum([tracks.meanacceleration].*[tracks.trajtime])/sum([tracks.trajtime]);
        [x,ix] = sort([tracks.meanacceleration]);
        w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
        mins = floor(min(x(w>=0.005))*100)/100;
        maxs = ceil(max(x(w<=0.995))*100)/100;
        s = (maxs-mins)/nbin;
        r = mins-s:s:maxs+s;
        a = zeros(numel(r),max(uclust));
        for j = 1:numel(r)-1
            for jj = 1:max(uclust)
                ii = [tracks.meanacceleration]<r(j+1) & [tracks.meanacceleration]>=r(j) & clust == jj;
                a(j,jj) = sum([tracks(ii).trajtime]);
            end
        end
        a = a/sum(a(:))/s;
        hb = bar(r'+s/2,a,'stacked');
        for jj = 1:max(uclust)
            hb(jj).BarWidth = 1;
            hb(jj).FaceColor = c(jj,:);
            hb(jj).EdgeColor = 'none';
        end
        xlim([mins maxs]);
        xlabel('Mean acceleration (\mum/s^2)','FontSize',12);
        ylabel('Probability density','FontSize',12);
        axis square;
        set(gca,'XMinorTick','on');
        %         title(strcat(sprintf('Mean: %g',round(mAcc,2,'significant')),' \mum/s^2'));
        
        %% Variance Swimming speed
        subplot(nrow,6,6),
        %         mStdS = sum(sqrt([tracks.varspeed])./[tracks.meanspeed].*[tracks.trajtime])/sum([tracks.trajtime]);
        [x,ix] = sort(sqrt([tracks.varspeed])./[tracks.meanspeed]);
        w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
        mins = floor(min(x(w>=0.005))*20)/20;
        maxs = ceil(max(x(w<=0.995))*20)/20;
        s = (maxs-mins)/nbin;
        r = mins-s:s:maxs+s;
        a = zeros(numel(r),max(uclust));
        for j = 1:numel(r)-1
            for jj = 1:max(uclust)
                
                ii = sqrt([tracks.varspeed])./[tracks.meanspeed]<r(j+1) & sqrt([tracks.varspeed])./[tracks.meanspeed]>=r(j)& clust == jj;
                a(j,jj) = sum([tracks(ii).trajtime]);
            end
        end
        a = a/sum(a(:))/s;
        hb = bar(r'+s/2,a,'stacked');
        for jj = 1:max(uclust)
            hb(jj).BarWidth = 1;
            hb(jj).FaceColor = c(jj,:);
            hb(jj).EdgeColor = 'none';
        end
        xlim([mins maxs]);
        xlabel('CV speed','FontSize',12);
        ylabel('Probability density','FontSize',12);
        axis square;
        set(gca,'XMinorTick','on');
        %         title(strcat(sprintf('Mean: %g',round(mStdS,2,'significant'))));
        
        %% Delta XY
        subplot(nrow,6,7), hold on,
        pmax = -Inf;
        for jj = 1:max(uclust(cumulclust<0.9))
            temp = tracks(clust == jj);
            
            dx = cat(1,temp.xvec)*mean(diff(temp(1).time));
            dy = cat(1,temp.yvec)*mean(diff(temp(1).time));
            dd = [dx; dy];
            %         tumble = cat(1,temp.tumble);
            %         if numel(tumble)==numel(dx)
            %             dx = dx(tumble==0);
            %             dy = dy(tumble==0);
            %         end
            smin = prctile(dd,0.005);
            smax = prctile(dd,0.995);
            smin = min(smin, -smax);
            smax = max(-smin, smax);
            r = (smax-smin)/100;
            s = smin-r:r:smax+r;
            [a,b] = histcounts(dd,s);
            b = (b(1:end-1) + b(2:end))/2;
            plot(b,a./sum(a)/r,'o','Markerfacecolor',c(jj,:),'MarkerEdgeColor','none','Markersize',4);
            
            plot(b, (pdf(makedist('normal','mu',nanmean(dd),'sigma', nanstd(dd)),b)),...
                ':','color',c(jj,:),'linewidth',1);
            %         [a,~] = histcounts(dy,s);
            %         plot(b,(a./sum(a)/r),'o','Markerfacecolor',cl(2,:),'MarkerEdgeColor','none','Markersize',4);
            %         plot(b, (pdf(makedist('normal','mu',nanmean(dy),'sigma', nanstd(dy)),b)),...
            %             ':','color',cl(2,:),'linewidth',1.5);
            pmax = max(pmax,10.^(ceil(2*log10(max((a./sum(a)/r))))/2));
        end
        
        ylim([pmax*10^-4 pmax]);
        xlabel('\DeltaXY (\mum)','FontSize',12);
        ylabel('Probability density','FontSize',12);
        axis square;
        %         legend({'X','Fit X','Y','Fit Y'},'location','south');
        set(gca,'XMinorTick','on');
        set(gca,'YScale','log');
        box on;
        xl = xlim;
        xlim([-min(abs(xl)) min(abs(xl))]);
        
        %% MSD
        subplot(nrow,6,9), hold on,
        t = round(cell2mat(cellfun(@(x) (x - x(1)),{tracks.time}','UniformOutput',0)),5);
        ut = unique(t);
        
        
        for jj = 1:max(uclust(cumulclust<0.9))
            meanmsd = zeros(numel(ut),1);
            w = zeros(numel(ut),1);
            temp = tracks(clust == jj);
            for j = 1:numel(temp)
                for i = 2:(numel(temp(j).x))
                    x_sq_disp = (temp(j).x(i:end) - temp(j).x(1:end-i+1)).^2;
                    y_sq_disp = (temp(j).y(i:end) - temp(j).y(1:end-i+1)).^2;
                    if sum(~isnan(x_sq_disp))>0
                        meanmsd(i) = ((meanmsd(i) * w(i)) + nansum(x_sq_disp + y_sq_disp))/(w(i)+sum(~isnan(x_sq_disp)));
                        w(i) = w(i)+sum(~isnan(x_sq_disp));
                    end
                end
            end
            
%             [~, xp, yp] = slmengine(log10(ut),log10(meanmsd),'minslope',0,...
%                 'plot','off','weights',w,'endconditions','notaknot',...
%                 'interiorknots','free','knots',6,'maxslope',2,...
%                 'order',4,'predictions',101);
%             plot(xp, yp,':','linewidth',1,'color',c(jj,:));
            
            scatter(flipud(log10(ut)),flipud(log10(meanmsd)),10,flipud(1-(1-c(jj,:)).*(0.05+(w/w(2))/1.05)),'Marker','o','MarkerFacecolor','flat',...
                'MarkeredgeColor','none');
        end
        
        ax = gca;
        ax.XTick = [-3 -2 -1 0 1 2 3 4 5 6 7];
        ax.XTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2' '10^3'...
            '10^4' '10^5' '10^6' '10^7'};
        ax.YTick = [-3 -2 -1 0 1 2 3 4 5 6 7];
        ax.YTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2' '10^3'...
            '10^4' '10^5' '10^6' '10^7'};
        xlim([floor(log10(ut(2))*2)/2 ceil(log10(ut(end))*2)/2]);
        xlabel('Time delay (s)','FontSize',12);
        ylabel('Mean sq. displacement (\mum^2)','FontSize',12);
        axis square;
        %         xl = xlim;
        yl = ylim;
        plot([-2 2], [yl(1) yl(1)+4],':','color','k');
        plot([-2 2], [yl(2)-2*4 yl(2)],':','color','k');
        %                 plot([xl(1) xl(2)], [yl(1) yl(1)+xl(2)-xl(1)],':','color','k');
        %         plot([xl(1) xl(2)], [yl(2)-2*(xl(2)-xl(1)) yl(2)],':','color','k');
        box on;
        %         legend({'Data','Fit'},'location','NorthWest');
        xlim([-2 2]);
        ylim(yl);
        
        %% VACF
        subplot(nrow,6,8),hold on;
        t = unique(round(cell2mat(cellfun(@(x) (x - x(1)),{tracks.time}','UniformOutput',0)),5));
        
        for jj = 1:max(uclust(cumulclust<0.9))
            temp = tracks(clust == jj);
            vacf = zeros(numel(t),1);
            n = zeros(numel(t),1);
            for i = 1:numel(temp)
                for j = 1:numel(temp(i).time)
                    d = dot([temp(i).xvec(1:(end-j+1)), ...
                        temp(i).yvec(1:(end-j+1))], ...
                        [temp(i).xvec(j:(end)), ...
                        temp(i).yvec(j:(end))],2);
                    if sum(~isnan(d))>0
                        vacf(j) = (vacf(j)*n(j) + sum(d(~isnan(d))))/(n(j) + sum(~isnan(d)));
                        n(j) = n(j) + sum(~isnan(d));
                    end
                end
            end
            %         vacfmsd = diff(diff(meanmsd))./diff(ut(1:end-1)).^2./2;
            %         vacf = nanmean(horzcat(vacf, [NaN; vacfmsd; NaN]),2);
            
            
            %         scatter(flipud(log10(t)),flipud(([NaN; vacfmsd; NaN]/vacf(1))),10,flipud(1-((1-c(2,:)).*(0.05+n/n(1)/1.05))),'Marker','o','MarkerFacecolor','flat',...
            %             'MarkeredgeColor','none');
            %         plot(ut(2:end-1),vacfmsd/vacf(1),'o','MarkerFacecolor',[0 0 0],...
            %             'MarkeredgeColor','none');
            
            f = fit(t(2:end),vacf(2:end),'a*exp(-x/b)*cos(x*c)',...
                'Weight',n(2:end), 'Startpoint',[vacf(1) 1 1],'Lower',[0 0 0],'Upper',[Inf Inf Inf]);
            
            %         xh = max(1,ceil(10*f.a));
            
            %         xlim([0 xh]);
            fout = feval(f,10.^(-2:0.04:2))/f.a;
            
            scatter(flipud(log10(t)),flipud(vacf)./f.a,10,flipud(1-((1-c(jj,:)).*(0.05+n/n(1)/1.05))),...
                'Marker','o','MarkerFacecolor','flat','MarkeredgeColor','none');
            plot(-2:0.04:2,fout,':','linewidth',1,'color',c(jj,:));
        end
        plot([-2 2], [0 0], 'k:');
        ax = gca;
        ax.XTick = [-2 -1 0 1 2];
        ax.XTickLabel = {'10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'};
        xlabel('Time delay (s)','FontSize',12);
        ylabel('Dir. autocorrelation','FontSize',12);
        axis square;
        set(gca,'XMinorTick','on');
        box on
        %         legend({'Data','Fit'},'location','NorthEast');
        %         title(strcat(sprintf('Persistence time: %g',round(f.a,2,'significant')),' s'));
        xlim([-2 2]);
        ylim([-0.4 1]);
        %% Alpha MSD
        if isfield(tracks, 'alphamsd')
            subplot(nrow,6,10),
            
            for jj = 1:max(uclust(cumulclust<0.9))
                
                t = round(cell2mat(cellfun(@(x) (x - x(1)),{tracks(clust == jj).time}','UniformOutput',0)),5);
                f = cell2mat(cellfun(@(x) (x - x(1)),{tracks(clust == jj).frozen}','UniformOutput',0));
                t = t(~f);
                alpha = vertcat(tracks(clust == jj).alphamsd);
                n = vertcat(tracks(clust == jj).diffcoeff_cve_n);
                ut = unique(t);
                meanalpha = zeros(numel(ut),1);
                w = zeros(numel(ut),1);
                
                for j = 1:numel(ut)
                    w(j) = nansum(n(t==ut(j)));
                    meanalpha(j) = nansum(n(t==ut(j)).*alpha(t==ut(j))) / w(j);
                end
                
%                 meanalpha(end) = NaN;
%                 if sum(~isnan(meanalpha))>1
%                     [~, xp, yp] = slmengine(-log10(ut),meanalpha,...
%                         'plot','off','weights',w/sum(w),'endconditions','notaknot',...
%                         'interiorknots','fixed','knots',-5:2/nbin:5,...
%                         'order',1,'predictions',101);
%                     plot(xp, yp,'-','linewidth',1.5,'color',c(jj,:));
                    
%                     hold on;
%                 end
            end
            plot([-3 2],[1 1],':','color','k');
            ax = gca;
            ax.XTick = [-3 -2 -1 0 1 2];
            ax.XTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'};
            ylim([0 2]);
            xlim([-2 2])
            %             xlim([floor(-log10(ut(end-1))*2)/2 ceil(-log10(ut(2))*2)/2]);
            xlabel('Frequency (s^{-1})','FontSize',12);
            ylabel('Anomalous diff. exp.','FontSize',12);
            axis square;
            %             legend({'Fit'},'location','NorthWest');
        end
        
        %% object intensity
        subplot(nrow,6,11),
        objintensity = [tracks.meanobjintensity];
        %         meanint = nansum(objintensity.*[tracks.trajtime])./sum([tracks.trajtime]);
        [x,ix] = sort(round(objintensity,2,'significant'));
        w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
        minint = floor(min(x(w>=0.005)));
        maxint = ceil(max(x(w<=0.995)));
        
        if minint==maxint
            minint = min(objintensity);
            maxint = max(objintensity);
        end
        
        s = (maxint-minint)/nbin;
        r = minint-s:s:maxint+s;
        a = zeros(numel(r),max(uclust));
        for j = 1:numel(r)-1
            for jj = 1:max(uclust)
                ii = objintensity<r(j+1) & objintensity>=r(j) & clust == jj;
                a(j,jj) = sum([tracks(ii).trajtime]);
            end
        end
        a = a/sum(a(:))/s;
        hb = bar(r'+s/2,a,'stacked');
        for jj = 1:max(uclust)
            hb(jj).BarWidth = 1;
            hb(jj).FaceColor = c(jj,:);
            hb(jj).EdgeColor = 'none';
        end
        xlim([minint maxint]);
        xlabel('Object intensity (AU)','FontSize',12);
        ylabel('Probability density','FontSize',12);
        axis square;
        %         title(strcat(sprintf('Mean: %g',round(meanint,2,'significant')),' AU'));
        
        %% CV object intensity
        subplot(nrow,6,12),
        stdobjintensity = sqrt([tracks.varobjintensity])./[tracks.meanobjintensity];
        meanint = sum(stdobjintensity.*[tracks.trajtime])./sum([tracks.trajtime]);
        [x,ix] = sort(round(stdobjintensity,2,'significant'));
        w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);
        minint = floor(min(x(w>=0.005))*20)/20;
        maxint = ceil(max(x(w<=0.995))*20)/20;
        
        s = (maxint-minint)/nbin;
        r = minint-s:s:maxint+s;
        a = zeros(numel(r),max(clust));
        for j = 1:numel(r)-1
            for jj = 1:max(clust)
                ii = stdobjintensity<r(j+1) & stdobjintensity>=r(j)& clust == jj;
                a(j,jj) = sum([tracks(ii).trajtime]);
            end
        end
        a = a/sum(a(:))/s;
        hb = bar(r'+s/2,a,'stacked');
        for jj = 1:max(uclust)
            hb(jj).BarWidth = 1;
            hb(jj).FaceColor = c(jj,:);
            hb(jj).EdgeColor = 'none';
        end
        xlim([minint maxint]);
        xlabel('CV object intensity','FontSize',12);
        ylabel('Probability density','FontSize',12);
        axis square;
        %         title(strcat(sprintf('Mean: %g',round(meanint,2,'significant'))));
        %%
        if nrow == 3

        %% Tumble bias
        if ~isempty([tracks.tumblebias]) && sum(diffclust>1)>0
            subplot(nrow,6,13),
            %             mTB = nansum([tracks.tumblebias].*[tracks.trajtime])/sum([tracks.trajtime]);
            s = 1/nbin;
            r = 0:s:1;
            a = zeros(numel(r),max(clust));
            for j = 1:numel(r)-1
                for jj = 1:max(clust)
                    if diffclust(jj) > 1
                        ii = [tracks.tumblebias]<r(j+1) & [tracks.tumblebias]>=r(j)& clust == jj;
                        a(j,jj) = sum([tracks(ii).trajtime]);
                    end
                end
            end
            a = a/sum(a(:))/s;
            hb = bar(r'+s/2,a,'stacked');
            for jj = 1:max(uclust)
                hb(jj).BarWidth = 1;
                hb(jj).FaceColor = c(jj,:);
                hb(jj).EdgeColor = 'none';
            end
            xlim([0 1]);
            ax = gca;
            ax.XTick = [0 0.2 0.4 0.6 0.8 1];
            ax.XTickLabel = {'0' '0.2' '0.4' '0.6' '0.8' '1'};
            xlabel('Tumble bias','FontSize',12);
            ylabel('Probability density','FontSize',12);
            axis square;
            set(gca,'XMinorTick','on');
            %             title(sprintf('Mean: %g',round(mTB,2,'significant')));
        end
        %% Tumble angle
        tumbleangles = abs(cat(2,tracks.tumbleangle));
        if numel(tumbleangles)>2 && sum(diffclust>1)>0

            subplot(nrow,6,14),
            s = 180/nbin;
            r = 0:s:180;
            a = zeros(numel(r),max(clust));
            for j = 1:numel(r)-1
                for jj = 1:max(clust)
                    if diffclust(jj) > 1
                        ii = [tracks(clust == jj).tumbleangle]<r(j+1) & [tracks(clust == jj).tumbleangle]>=r(j);
                        a(j,jj) = sum(ii);
                    end
                end
            end
            a = a/sum(a(:))/s;
            hb = bar(r'+s/2,a,'stacked');
            for jj = 1:max(uclust)
                hb(jj).BarWidth = 1;
                hb(jj).FaceColor = c(jj,:);
                hb(jj).EdgeColor = 'none';
            end
            xlim([0 180]);
            xlabel('Tumble angle (deg)','FontSize',12);
            ylabel('Probability density','FontSize',12);
            %             title(sprintf('Mean: %g deg.',round(mean(tumbleangles),2,'significant')));
            axis square;
        end
        %% Tumble time

        if numel([tracks.tumbletime])>2 && sum(diffclust>1)>0
            subplot(nrow,6,15),
            s = 3/nbin;
            r = -2:s:1;
            a = zeros(numel(r),max(clust));
            for j = 1:numel(r)-1
                for jj = 1:max(clust)
                    if diffclust(jj) > 1
                        ii = log10([tracks(clust == jj).tumbletime])<r(j+1) & log10([tracks(clust == jj).tumbletime])>=r(j);
                        a(j,jj) = sum(ii);
                    end
                end
            end
            a = a/sum(a(:))/s;
            hb = bar(r'+s/2,a,'stacked');
            for jj = 1:max(uclust)
                hb(jj).BarWidth = 1;
                hb(jj).FaceColor = c(jj,:);
                hb(jj).EdgeColor = 'none';
            end
            xlim([-2 1]);
            xlabel('Tumble time (s)','FontSize',12);
            ylabel('Probability density','FontSize',12);
            ax = gca;
            ax.XTick = [-2 -1 0 1];
            ax.XTickLabel = {'10^{-2}' '10^{-1}' '10^0' '10^1'};
            %             title(sprintf('Mean: %g s',round(mean([tracks.tumbletime]),2,'significant')));
            axis square;
        end
        %% Rev freq

        if numel([tracks.revfreq])>2 && sum(diffclust>1)>0

            subplot(nrow,6,16),
            %             ind = ~isnan([tracks.revfreq]);
            %             meanrev = sum([tracks(ind).revfreq].*[tracks(ind).trajtime])./sum([tracks(ind).trajtime]);
            
            [x,ix] = sort(round([tracks.revfreq],2,'significant'));
            w = cumsum([tracks(ix).trajtime])/sum([tracks(ix).trajtime]);

            minint = floor(min(x(w>=0.005))*10)/10;
            maxint = ceil(max(x(w<=0.995))*10)/10;

            if maxint<=minint
                minint = min(x);
                maxint = max(x);
            end
            s = (maxint-minint)/nbin;
            r = minint-s:s:maxint+s;
            a = zeros(numel(r),max(uclust));
            
            for j = 1:numel(r)-1
                for jj = 1:max(uclust)
                    if diffclust(jj) > 1
                        ii = [tracks.revfreq]<r(j+1) & [tracks.revfreq]>=r(j) & clust ==jj;
                        a(j,jj) = sum([tracks(ii).trajtime]);
                    end
                end
            end
            a = a/sum(a(:))/s;
            hb = bar(r'+s/2,a,'stacked');
            for jj = 1:max(uclust)
                hb(jj).BarWidth = 1;
                hb(jj).FaceColor = c(jj,:);
                hb(jj).EdgeColor = 'none';
            end
            xlim([minint maxint]);
            xlabel('Reversion frequency (s^{-1})','FontSize',12);
            ylabel('Probability density','FontSize',12);
            axis square;
            %             title(strcat(sprintf('Mean: %g',round(meanrev,2,'significant')),' s^{-1}'));
        end
        %% Diff coeff th
        if ~isempty([tracks.tumblebias])
            subplot(nrow,6,18),
            %scatter([tracks.tumblebias],[tracks.revfreq],'.')
            ind = ~isnan([tracks.diffcoeff]);
            if sum(ind)>2
            %             diffth = (0.5.*[tracks(ind).meanspeed].^2.*[tracks(ind).meanruntime]);
            scatter(log10([tracks(ind).diffcoeff_cve_mean]),log10([tracks(ind).diffcoeff]),[],c(clust(ind),:),'Marker','.');
            hold on,
            plot([-5 8],[-5 8],':','color','k');
            xlim([floor(2*min(mindiff,min([tracks(ind).diffcoeff])))/2 ceil(2*max(maxdiff,max(min(mindiff,min([tracks(ind).diffcoeff]))))/2)]);
            ylim([floor(2*min(mindiff,min([tracks(ind).diffcoeff])))/2 ceil(2*max(maxdiff,max(min(mindiff,min([tracks(ind).diffcoeff]))))/2)]);
            ax = gca;
            ax.XTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
            ax.XTickLabel = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}'...
                '10^0' '10^1' '10^2' '10^3' '10^4' '10^5' '10^6' '10^7' '10^8'};
            ax.YTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
            ax.YTickLabel = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}'...
                '10^0' '10^1' '10^2' '10^3' '10^4' '10^5' '10^6' '10^7' '10^8'};
            xlabel('Eff. Diff. coeff. (\mum^2/s)','FontSize',12);
            ylabel('Th. Diff. coeff. (\mum^2/s)','FontSize',12);
            axis square;
            end
        end
        
        %% TB rev freq
        if ~isempty([tracks.tumblebias])
            yh = ceil(quantile([tracks.revfreq],0.995));
            subplot(nrow,6,17),
            scatter([tracks.tumblebias],[tracks.revfreq],[],c(clust,:),'Marker','.');
            xlim([0 1]);
            ylim([0 yh]);
            xlabel('Tumble bias','FontSize',12);
            ylabel('Reversion freq. (s^{-1})','FontSize',12);
            axis square;
        end
        end
        %%
        set(h,'PaperPositionMode','auto');
        %         set(h,'PaperOrientation','landscape');
        print(h,strcat(tracks(1).metadata.name_,'.png'),'-dpng','-r150','-noui');
        crop(strcat(tracks(1).metadata.name_,'.png'));
%         print(h,strcat(tracks(1).metadata.name_,'.svg'),'-dsvg','-r150','-noui');
        fprintf('done.\n');
    else
        fprintf('No track to plot!\n');
    end
end
end
