function tracks = filterAnomalousTracks(tracks, window, msd_rep, fdr,diag_plot)
% This filter calculates the square displacement of cells 'sdisp' over a time 
% window 'window' (1 second default) that is converted to number of frames 'w'.
% The square displacement 'sdisp' is normalized by the mean square 
% displacement 'msdisp' and smoothed with a moving average over a window 'w2'. 
% 'w2' is calculated such that the average is done with 'msd_rep' number of 
% data points and that 'w2' is an odd number to keep the symmetry in the data 
% vector.
% The log of the normalized displacement is fitted with a gaussian and the 
% false discovery rate 'fdr' is used to determine a cutoff for anomalous 
% diffusion.

% window = 1;
% msd_rep = 10;
% fdr = 0.05;
% diag_plot = 1;
if numel(tracks)<2
    return
end
%recommened setting window = 1 sec, msd_rep = 10, fdr =0.05
w = round(window/(tracks(1).time(2)-tracks(1).time(1)));
w2 = 2*floor((w+msd_rep)/2)+1;

ind1 =[];  % good tracks
ind2 = []; % filtered tracks
ind3 = []; % short tracks

msd = cell(numel(tracks),1);

for i = 1:numel(tracks)
    if numel(tracks(i).x) > w2
        %         dx = (tracks(i).x(w:end)-tracks(i).x(1:end-w+1))' ./ tracks(i).meanspeed;
        %         dy = (tracks(i).y(w:end)-tracks(i).y(1:end-w+1))' ./ tracks(i).meanspeed;
        
        dx = (tracks(i).x(w:end)-tracks(i).x(1:end-w+1))';
        dy = (tracks(i).y(w:end)-tracks(i).y(1:end-w+1))';
        
        dx = [ones(1,w-1)*NaN dx];
        dy = [ones(1,w-1)*NaN dy];
        
        dx(logical(tracks(i).frozen)) = NaN;
        dy(logical(tracks(i).frozen)) = NaN;
        
        sdisp = sqrt(dx.^2+dy.^2);
        msdisp = mean(sdisp(~isnan(sdisp)));
        
        %         msdtmp = tsmovavg((0.5*(dx.^2+dy.^2)./(tracks(i).time(w)-tracks(i).time(1))'), 's',w2);
%         msdtmp = tsmovavg(sdisp, 's',w2);
%         msd{i} = [msdtmp(ceil(w2/2):end) ones(1,(w2-1)/2)*NaN];
msd{i} = std(sdisp(~isnan(sdisp)))/msdisp;
    end
end
%

allmsd = ([msd{:}])';
allmsd = allmsd(~isnan(allmsd));
try
    opts = statset('MaxIter',500);
    obj2 = gmdistribution.fit(allmsd,2,'Options',opts);
%     [~, b] = max(obj2.ComponentProportion);
%     obj = gmdistribution(obj2.mu(b),obj2.Sigma(1,1,b),1);
%     idx = cluster(obj2,allmsd);
catch
    return
end
%
% anomalous = cell(numel(tracks),1);
anomalous = zeros(numel(tracks),1);

for i = 1:numel(tracks);
    if ~isempty(msd{i})
        idx = cluster(obj2,log10((msd{i}')));
        anomalous(i) = sum(idx==1)/numel(idx);
    end
end

hist(anomalous,100);

tracks = tracks(ind1);
% for i = 1:numel(tracks);
%     if ~isempty(msd{i})
%         x = cdf(obj,log10(sqrt(msd{i}')));
%         y  = x * sum(~isnan(x));
%         z = zeros(numel(x),1);
%         
%         for j = 1:numel(x)
%             z(j) = y(j)/sum(x<=x(j));
%         end
%         
%         anomalous{i} = z <= fdr;
%         
%         if sum(z<=fdr)==0
%             ind1 = [ind1 i];
%         else
%             ind2= [ind2 i];
%         end
%     else
%         ind3 = [ind3 i];
%     end
% end
% 
% tracks = tracks(ind1);

% %%
% if diag_plot
%     
%     figure,
%     subplot(3,2,1),
%     hist(log10([msd{:}]),100);
%     %     xlim([0 2]);
%     
%     subplot(3,2,2),
%     qqplot(log10([msd{:}]),random('norm',0,1,100,1));
%     %     xlim([0 2]);
%     
%     subplot(3,2,3),
%     hist(log10([msd{ind1}]),100);
%     %     xlim([0 2]);
%     
%     subplot(3,2,4),
%     qqplot(log10([msd{ind1}]),random('norm',0,1,100,1));
%     %     xlim([0 2]);
%     
%     subplot(3,2,5),
%     hist(log10([msd{ind2}]),100);
%     %     xlim([0 2]);
%     
%     figure, hold on,
%     for i = ind1
%         plot(tracks(i).x,tracks(i).y,'color',random('unif',0,1,1,3));
%     end
%     axis tight
%     axis equal
%     
%     figure, hold on,
%     for i = ind2
%         plot(tracks(i).x,tracks(i).y,'color',random('unif',0,1,1,3));
%     end
%     axis tight
%     axis equal
%     
%     figure, hold on,
%     for i = ind3
%         plot(tracks(i).x,tracks(i).y,'color',random('unif',0,1,1,3));
%     end
%     axis tight
%     axis equal
%     
%     if 0
%         figure,
%         for i = ind2
%             
%             h1 = subplot(2,1,1);
%             plot(tracks(i).x,tracks(i).y,'k')
%             hold on
%             scatter(tracks(i).x(anomalous{i}),tracks(i).y(anomalous{i}),'r.')
%             scatter(tracks(i).x(logical(tracks(i).frozen)),tracks(i).y(logical(tracks(i).frozen)),'b')
%             hold off
%             axis equal
%             axis tight
%             p = get(h1, 'pos');
%             p(1) = p(1) - 0.05;
%             p(2) = p(2) - 0.05;
%             p(3) = p(3) + 0.1;
%             p(4) = p(4) + 0.1;
%             set(h1, 'pos', p);
%             
%             h2 = subplot(2,1,2);
%             plot(tracks(i).time,tracks(i).speed,'k')
%             hold on
%             scatter(tracks(i).time(anomalous{i}),tracks(i).speed(anomalous{i}),'r.')
%             scatter(tracks(i).time(logical(tracks(i).frozen)),tracks(i).speed(logical(tracks(i).frozen)),'b')
%             hold off
%             p = get(h2, 'pos');
%             p(1) = p(1) - 0.05;
%             p(2) = p(2) - 0.05;
%             p(3) = p(3) + 0.1;
%             p(4) = p(4) + 0.1;
%             set(h2, 'pos', p);
%             pause();
%         end
%     end
% end
% %% rescue anomalous tracks
% ii = 1;
% new_tracks =cell(1,1);
% for i = ind2
%     pos = find(anomalous{i});
%     pos = [1-w2; pos; numel(tracks(i).x)+w2];
%     l = find(pos(2:end) - pos(1:end-1)> 4*w2);
%     
%     for j = 1:numel(l)
%         new_track = chop_track(tracks(i), pos(l(j))+w2, pos(l(j)+1)-w2);
%         
%         if 0
%             h1 = subplot(2,1,1);
%             plot(tracks(i).x,tracks(i).y,'k')
%             hold on
%             plot(new_track.x,new_track.y,'g')
%             scatter(tracks(i).x(anomalous{i}),tracks(i).y(anomalous{i}),'r.')
%             hold off
%             axis equal
%             axis tight
%             p = get(h1, 'pos');
%             p(1) = p(1) - 0.05;
%             p(2) = p(2) - 0.05;
%             p(3) = p(3) + 0.1;
%             p(4) = p(4) + 0.1;
%             set(h1, 'pos', p);
%             
%             h2 = subplot(2,1,2);
%             plot(tracks(i).time,tracks(i).speed,'k')
%             hold on
%             plot(new_track.time,new_track.speed,'g')
%             scatter(tracks(i).time(anomalous{i}),tracks(i).speed(anomalous{i}),'r.')
%             scatter(tracks(i).time(logical(tracks(i).frozen)),tracks(i).speed(logical(tracks(i).frozen)),'b')
%             hold off
%             p = get(h2, 'pos');
%             p(1) = p(1) - 0.05;
%             p(2) = p(2) - 0.05;
%             p(3) = p(3) + 0.1;
%             p(4) = p(4) + 0.1;
%             set(h2, 'pos', p);
%             pause();
%         end
%         new_tracks{ii} = new_track;
%         ii = ii +1;
%     end
% end
% new_tracks(cellfun('isempty',new_tracks))=[];
% new_tracks = horzcat(new_tracks{:});
% %%
% tracks = horzcat(tracks(ind1), new_tracks);
end
