% smooth cell trajectories
function s_traj = smooth_traj(traj)
%     warning off
s_traj = traj;

for i = 1:numel(s_traj)

    
%     cs = spaps(s_traj{i}(:,11),[s_traj{i}(:,1) s_traj{i}(:,2)]',numel(s_traj{i}(:,11))*tol,2);
%     s = fnval(cs, s_traj{i}(:,11))';
    
    X = [s_traj{i}(:,11), s_traj{i}(:,1), s_traj{i}(:,2)];
    X = fillmissing(X,'makima');
    
    s_traj{i}(:,1) = X(:,2);
    s_traj{i}(:,2) = X(:,3);
    
    % calculate xy vectors
    timediff = [s_traj{i}(2:end,11)-s_traj{i}(1:(end-1),11); NaN];
    
    xvec = s_traj{i}(2:end,1)-s_traj{i}(1:(end-1),1);
    yvec = s_traj{i}(2:end,2)-s_traj{i}(1:(end-1),2);
    
    s_traj{i}(:,7) = [xvec; NaN]./timediff;
    s_traj{i}(:,8) = [yvec; NaN]./timediff;
    
    % calculate speed
    s_traj{i}(:,9) = sqrt(s_traj{i}(:,7).^2 + s_traj{i}(:,8).^2);
    
    angle = atan2(s_traj{i}(1:end-1,7).*s_traj{i}(2:end,8)- ...
            s_traj{i}(1:end-1,8).*s_traj{i}(2:end,7), ...
            s_traj{i}(1:end-1,7).*s_traj{i}(2:end,7)+ ...
            s_traj{i}(1:end-1,8).*s_traj{i}(2:end,8));
    
    s_traj{i}(:,10) = [NaN; angle];

    % divide angle by time to get angular velocity
    s_traj{i}(:,10) = s_traj{i}(:,10)./timediff;
   
end
