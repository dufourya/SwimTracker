function tracks = calculateTumbleBasedStatistics(tracks)
    nTracks = numel(tracks);
    for i = 1:nTracks
        times = [tracks(i).time(2:end)-tracks(i).time(1:(end-1)); ...
                 tracks(i).time(2)-tracks(i).time(1)];
        
        tracks(i).tumblebias = sum(tracks(i).tumble==1) / ...
                               sum(~isnan(tracks(i).tumble));
        
        tracks(i).tumbletime  = NaN;
        tracks(i).runtime     = NaN;
        tracks(i).runlength   = NaN;
        tracks(i).tumbleangle = NaN;
        
        tumbleind = find(tracks(i).tumble==1);
        tumbletime = 0;
        
        if length(tumbleind)>1
            k=1;
            for j=2:length(tumbleind)
                tumbletime = tumbletime + (tracks(i).time(tumbleind(j-1)) - ...
                                           tracks(i).time(tumbleind(j-1)-1));
                if (tumbleind(j)-tumbleind(j-1))>1
                    tracks(i).tumbletime = [tracks(i).tumbletime tumbletime];
                    tumbletime = 0;
                    runtime = tracks(i).time(tumbleind(j)) - ...
                              tracks(i).time(tumbleind(j-1)+1);
                    tracks(i).runtime = [tracks(i).runtime runtime];
                    runlength = ...
                        sum(tracks(i).speed(tumbleind(j-1):tumbleind(j)-1) .* ...
                        times(tumbleind(j-1):tumbleind(j)-1));
                    tracks(i).runlength = [tracks(i).runlength, runlength];
                    tracks(i).tumbleangle = [tracks(i).tumbleangle, ...
                        180/pi*atan2(tracks(i).xvec(tumbleind(k)-1)* ...
                        tracks(i).yvec(tumbleind(j-1)+1)- ...
                        tracks(i).yvec(tumbleind(k)-1)* ...
                        tracks(i).xvec(tumbleind(j-1)+1), ...
                        tracks(i).xvec(tumbleind(k)-1)* ...
                        tracks(i).xvec(tumbleind(j-1)+1)+ ...
                        tracks(i).yvec(tumbleind(k)-1)* ...
                        tracks(i).yvec(tumbleind(j-1)+1))];
                    k=j;
                end
            end
        elseif length(tumbleind)==1
            if tumbleind == numel(tracks(i).time)
                tracks(i).tumbletime = NaN;
                tracks(i).tumbleangle = NaN;
            else
                tracks(i).tumbletime = (tracks(i).time(tumbleind(1)+1)- ...
                    tracks(i).time(tumbleind(1)));
                tracks(i).tumbleangle = [tracks(i).tumbleangle, ...
                    180/pi*atan2(tracks(i).xvec(tumbleind(1)-1)* ...
                    tracks(i).yvec(tumbleind(1)+1)- ...
                    tracks(i).yvec(tumbleind(1)-1)* ...
                    tracks(i).xvec(tumbleind(1)+1), ...
                    tracks(i).xvec(tumbleind(1)-1)* ...
                    tracks(i).xvec(tumbleind(1)+1)+ ...
                    tracks(i).yvec(tumbleind(1)-1)* ...
                    tracks(i).yvec(tumbleind(1)+1))];
            end
        end
        
        tumble = tracks(i).tumble;
        tumble(isnan(tumble)) = [];
        tracks(i).meanrunspeed = nanmean(tracks(i).speed(tumble == 0));
        tracks(i).runspeed     = tracks(i).runlength./tracks(i).runtime;
        tracks(i).meanruntime  = nanmean(tracks(i).runtime);
        tracks(i).revfreq      = 2*(length(tracks(i).runtime)+1)/tracks(i).trajtime;
        tracks(i).diffcoeff    = 0.5.*tracks(i).meanruntime.*tracks(i).meanspeed.^2;
    end
end