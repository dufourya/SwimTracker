% smooth cell trajectories
function s_traj = smooth_traj(traj)
    warning off
    s_traj = traj;
    
    parfor i = 1:numel(s_traj)
        % smooth trajectory
        warning('off','all');
        cs = csaps(s_traj{i}(:,11),[s_traj{i}(:,1) s_traj{i}(:,2)]');
        warning('on','all');
        s = fnval(cs, s_traj{i}(:,11))';
        
        s_traj{i}(:,1) = s(:,1);
        s_traj{i}(:,2) = s(:,2);
        
        % calculate xy vectors
        times = [s_traj{i}(2:end,11)-s_traj{i}(1:(end-1),11); 1];
        xvec = s_traj{i}(2:end,1)-s_traj{i}(1:(end-1),1);
        yvec = s_traj{i}(2:end,2)-s_traj{i}(1:(end-1),2);
        
        s_traj{i}(:,7) = [xvec; 0]./times;
        s_traj{i}(:,8) = [yvec; 0]./times;
        
        s_traj{i}(end,7) = s_traj{i}(end-1,7);
        s_traj{i}(end,8) = s_traj{i}(end-1,8);
        
        % calculate speed
        s_traj{i}(:,9) = sqrt(s_traj{i}(:,7).^2 + s_traj{i}(:,8).^2);
        
        % calculate change in angle
        angles=zeros(length(s_traj{i}(:,7)),1);
        for j=2:length(angles),
            angles(j) = atan2(s_traj{i}(j-1,7)*s_traj{i}(j,8)-s_traj{i}(j-1,8)*s_traj{i}(j,7) , s_traj{i}(j-1,7)*s_traj{i}(j,7)+s_traj{i}(j-1,8)*s_traj{i}(j,8));
        end
        s_traj{i}(:,10) = angles;
        s_traj{i}(1,10) = angles(2);
    end
