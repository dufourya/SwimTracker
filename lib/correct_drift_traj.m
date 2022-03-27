function [trajectories, movie_drift] = correct_drift_traj(trajectories,correct)
%%
diag_plot = 0;

% temp = trajectories(cellfun(@(x) size(x,1), trajectories)>3);
% dx = cellfun(@(x) x(end,1)-x(1,1), temp);
% dy = cellfun(@(x) x(end,2)-x(1,2), temp);
% dt = cellfun(@(x) x(end,11)-x(1,11), temp);

% xvec = sum(dx)/sum(dt);
% yvec = sum(dy)/sum(dt);

% %%
xvec = cell2mat(cat(1,cellfun(@(x) squeeze(x(:,7)),trajectories,'unif',0)));
yvec = cell2mat(cat(1,cellfun(@(x) squeeze(x(:,8)),trajectories,'unif',0)));
frozen = cell2mat(cat(1,cellfun(@(x) squeeze(x(:,12)),trajectories,'unif',0)));
% nframes = numel(unique(cell2mat(cat(1,cellfun(@(x) squeeze(x(:,11)),trajectories,'unif',0)))));
% 
xvec1 = (xvec(~isnan(xvec) & ~isnan(yvec) & ~frozen));
yvec1 = (yvec(~isnan(xvec) & ~isnan(yvec) & ~frozen));
% 
%factor = 1/(1+1/((numel(xvec1)/nframes)/10)^2);
factor = 1;
% 
% [f, xi] = ksdensity([xvec1,yvec1]);
% [~, fi] = max(f);
% dx = max(diff(xi(:,1)));
% dy = max(diff(xi(:,2)));
% [X, Y] = meshgrid(xi(fi,1)-dx:dx/100:xi(fi,1)+dx, xi(fi,2)-dy:dy/100:xi(fi,2)+dy);
% f = ksdensity([xvec1,yvec1],[X(:), Y(:)]);
% [~, fi] = max(f);
% 
xvec = mean(xvec1)* factor;
yvec = mean(yvec1)* factor;
% 
if diag_plot
    c = lines(10);
    figure,
    hold on,
    for i=1:10
        plot(trajectories{i}(:,1),trajectories{i}(:,2),'color',c(i,:));
    end
end
% p = atan2(yvec, xvec);
% v = sqrt(xvec.^2 + yvec.^2);
% pr = round(p,2,'decimals');
%
% a = unique(pr(~isnan(pr)));
% va = zeros(numel(a),1);
%
% for i = 1:numel(a)
%     va(i) = mean(v(pr == a(i)));
% end
%
% xvec = sum(cos(a).*va)/629;
% yvec = sum(sin(a).*va)/629;

%%
movie_drift = struct('x',xvec,'y',yvec,'total',sqrt(xvec^2+yvec^2));

if correct
    fprintf('Correcting drift...');
    for i = 1:numel(trajectories)
        trajectories{i}(:,1) = trajectories{i}(:,1) - trajectories{i}(:,11)*xvec;
        trajectories{i}(:,2) = trajectories{i}(:,2) - trajectories{i}(:,11)*yvec;
        
        % delta t for each pair of tracked frames.
        timediff = [trajectories{i}(2:end,11) - trajectories{i}(1:end-1,11); ...
            NaN];
        
        % xvec: the x component of the velocity vector.
        trajectories{i}(:,7) = [trajectories{i}(2:end,1)-trajectories{i}(1:end-1,1); ...
            NaN]./timediff;
        
        % yvec: the y component of the velocity vector.
        trajectories{i}(:,8) = [trajectories{i}(2:end,2)-trajectories{i}(1:end-1,2); ...
            NaN]./timediff;
        
        % speed: the magnitude of the velocity vector.
        trajectories{i}(:,9) = sqrt(trajectories{i}(:,7).^2 + trajectories{i}(:,8).^2);
        
        % angle
        angle = atan2(trajectories{i}(1:end-1,7).*trajectories{i}(2:end,8)- ...
            trajectories{i}(1:end-1,8).*trajectories{i}(2:end,7), ...
            trajectories{i}(1:end-1,7).*trajectories{i}(2:end,7)+ ...
            trajectories{i}(1:end-1,8).*trajectories{i}(2:end,8));
    
        trajectories{i}(:,10) = [NaN; angle];
        % divide angle by time to get angular velocity
        trajectories{i}(:,10) = trajectories{i}(:,10)./timediff;
    end
end

if diag_plot
    for i=1:10
        plot(trajectories{i}(:,1),trajectories{i}(:,2),':','color',c(i,:));
    end
end

end
