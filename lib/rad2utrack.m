function rad_detection = rad2utrack(rad_objs, meta, file_info)
    % column key
    X = 1;
    Y = 2;
    AMP = 3;
    PARTICLE_ID = 4;
    FRAME = 5;
    TRACK_ID = 6;
    SIGMA = 7;
    
    interval_sec = meta.getImageIntervalMs() ./ 1000;
    
    is_cell = false;
    if iscell(rad_objs)
        is_cell = true;
    elseif ismatrix(rad_objs)
        is_cell = false;
    else
        error('rad_objs is incorrect type.')
    end
    
    n_frames = 0;
    if is_cell
        n_frames = size(rad_objs,1);
    else
        n_frames = max(rad_objs(FRAME,:));
    end

    rad_detection(n_frames).xCoord = [];
    rad_detection(n_frames).yCoord = [];
    rad_detection(n_frames).amp    = [];

    for i = 1:n_frames
        if is_cell
            frame_data = rad_objs{i};
        else
            frame_data = rad_objs(:, ismember(rad_objs(FRAME,:), i));
        end

        %frame_data = rad_objs{i};
        n_objs = size(frame_data,2);
        rad_detection(i).xCoord = [frame_data(X,:)', zeros(n_objs,1)];
        rad_detection(i).yCoord = [frame_data(Y,:)', zeros(n_objs,1)];
        rad_detection(i).amp    = [frame_data(AMP,:)', zeros(n_objs,1)];
        rad_detection(i).time   = interval_sec .* (i-1);
    end
end