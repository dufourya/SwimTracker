function traj0 = uTrack2mat(tracksFinal, movieInfo, pixel_size,n_before_frozen)
    %
    % Converts the output of u-track to something sane.
    % The following is copied from the u-track documentation
    % (TrackingProcess.doc):
    %
    %    tracksCoordAmpCG: The positions and amplitudes of the tracked 
    %       particles, after gap closing.
    %
    %    Number of rows: Number of tracks merging with each other and 
    %       splitting from each other (i.e., involved in compound track).
    %       Number of columns  : 8 × number of frames the compound track spans.
    %       For every frame, the matrix stores the particle’s x-coordinate,
    %       y-coordinate, z-coordinate (0 if 2D), amplitude, x-coordinate 
    %       standard deviation, y-coordinate standard deviation, z-coordinate
    %       standard deviation (0 if 2D) and amplitude standard deviation.
    %       NaNs indicate frames where a track does not exist, either because 
    %       those frames are before the track starts or after it ends, or 
    %       because of temporary particle disappearance.
    %
    %    tracksFeatIndxCG: Connectivity matrix of particles between frames,
    %       after gap closing.  Number of rows = Number of tracks merging with 
    %       each other and splitting from each other (i.e., involved in 
    %       compound track).  Number of columns = Number of frames the 
    %       compound track spans.  Zeros indicate frames where a track does not
    %       exist, either because those frames are before the track starts or 
    %       after it ends, or because of temporary particle disappearance.
    %
    %    seqOfEvents: Matrix storing the sequence of events in a compound track
    %        (i.e. track start, track end, track splitting and track merging).
    %        Number of rows = number of events in a compound track.  Number of
    %        columns = 4.  In every row, the columns mean the following:
    %        1st column indicates frame index where event happens.
    %        2nd column indicates whether the event is the start of end of a 
    %        track.  1 = start, 2 = end.  3rd column indicates the index of
    %        the track that starts or ends (The index is “local�?, within the
    %        compound track. It corresponds to the track’s row number in
    %        tracksfeatIndxCG and tracksCoordAmpCG).  4th column indicates 
    %        whether a start is a true initiation or a split, and whether an 
    %        end is a true termination of a merge. In particular, if the 4th
    %        column is NaN, then a start is an initiation and an end is a 
    %        termination. If the 4th column is a number, then the start is a
    %        split and the end is a merge, where the track of interest splits
    %        from / merges with the track indicated by the number in the 4th
    %        column.
     
    traj0 = cell(numel(tracksFinal),1);
    time = sort([movieInfo.time]);
    time_frozen = time(n_before_frozen);
	
    for i = 1:numel(traj0)
        nrow = numel(tracksFinal(i).tracksCoordAmpCG)/8;
        traj0{i} = zeros(nrow, 12);
        
        % x position
        traj0{i}(:,1) = pixel_size*tracksFinal(i).tracksCoordAmpCG(1:8:8*nrow);
        
        % y position
        traj0{i}(:,2) = pixel_size*tracksFinal(i).tracksCoordAmpCG(2:8:8*nrow);

        % intensity
        traj0{i}(:,4) = tracksFinal(i).tracksCoordAmpCG(4:8:8*nrow);
        
        % tracks feature index
        traj0{i}(:,6) = tracksFinal(i).tracksFeatIndxCG;
        
        % vector of times particle was tracked.
        start = tracksFinal(i).seqOfEvents(1,1);
        stop = tracksFinal(i).seqOfEvents(end,1);
        traj0{i}(:,11) = time(start:stop);
        
        % delta t for each pair of tracked frames.
        timediff = [traj0{i}(2:end,11) - traj0{i}(1:end-1,11); ...
                    traj0{i}(end,11) - traj0{i}(end-1,11)];
        
        % xvec: the x component of the velocity vector.
        traj0{i}(:,7) = [traj0{i}(2:end,1)-traj0{i}(1:end-1,1); ...
                         traj0{i}(end,1)-traj0{i}(end-1,1)]./timediff;
       
        % yvec: the y component of the velocity vector.
        traj0{i}(:,8) = [traj0{i}(2:end,2)-traj0{i}(1:end-1,2); ...
                         traj0{i}(end,2)-traj0{i}(end-1,2)]./timediff;
       
        % speed: the magnitude of the velocity vector.
        traj0{i}(:,9) = sqrt(traj0{i}(:,7).^2 + traj0{i}(:,8).^2);
        
        % angle
        for j = 2:nrow
            traj0{i}(j,10) = atan2(traj0{i}(j-1,7)*traj0{i}(j,8)- ...
                                   traj0{i}(j-1,8)*traj0{i}(j,7), ...
                                   traj0{i}(j-1,7)*traj0{i}(j,7)+ ...
                                   traj0{i}(j-1,8)*traj0{i}(j,8));
        end
        
        % duplicate the second angle into the first position, which has an
        % undefined angle.
        traj0{i}(1,10) = traj0{i}(2,10);
        
        % add frozen frame logical to ignore frames in downstream analysis
        traj0{i}(:,12) = traj0{i}(:,11) > time_frozen;
    end
end
