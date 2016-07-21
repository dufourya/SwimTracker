function runTrackSingle(movie_name, params)
    % 1. Detects cells using particle_center
    %    (http://physics.uoregon.edu/~raghu/particle_tracking.html) or
    %    u-track and tracks cells using u-track
    %    (http://lccb.hms.harvard.edu/software.html).
    % 2. Converts the output of u-track to a more intuitive format, which
    %    is used by runTumbleDetection to analyze runs and tumbles.
    % 3. Runs runTumbleDetectionSingle with default parameters.
    %
    % INPUTS:
    % movie_file   : the name of the movie to analyze.
    % params       : a structure defined in parameters.m
        
    % rename variables for compactness
    mem_usage        = params.runTrackSingle.mem_usage;
    max_object_count = params.runTrackSingle.maximum_object_count;

    if mem_usage > 1
        error('Can''t request more than the available amount of memory.')
    end

    close all
     
    startMatlabPool();

    saveResults.dir = '.';
        
    file_info = RadFileInfo(movie_name);
    
    fprintf('* %s\n',file_info.name);
    
    if ~exist(file_info.metaFile,'file');
        fprintf('Missing metadata file, skipping.\n');
        return
    else
        metadata = Metadata(file_info.metaFile);
        metadata.read();
        pixel_size = metadata.getPixelSize();
        n_images = metadata.getNumberOfImages();
        n_before_frozen = metadata.getNImagesBeforeFreeze();
        if isempty(n_before_frozen)
            n_before_frozen = n_images;
        end
    end
    
    try
        if exist(file_info.detectionFile, 'file') && ...
           exist(file_info.trackingFile, 'file') && ...
           exist(file_info.swimtrackerFile, 'file')
            swimtracker = load(file_info.swimtrackerFile);
            if isempty([swimtracker.tracks.tumblebias])
                runTumbleDetectionSingle(movie_name, params);
                return;
            else
                fprintf('nothing to do\n');
                return;
            end
        elseif ~exist(file_info.detectionFile,'file');
            
            detection_res = radDetect(file_info, metadata, params);

            % check the number of detected objects.
            if tooManyObjects(detection_res, max_object_count)
                return;
            end

            [tracksFinal, params] = track(detection_res, file_info, metadata, params);           
        elseif ~exist(file_info.trackingFile,'file');
            dat = load(file_info.detectionFile, file_info.detectionDat);
            detection_res = dat.(file_info.detectionDat);
            if tooManyObjects(detection_res, max_object_count)
                return;
            end
            [tracksFinal, params] = track(detection_res, file_info, metadata, params);            
        elseif ~exist(file_info.swimtrackerFile,'file');
            detection_dat = load(file_info.detectionFile, ...
                                 file_info.detectionDat);
            tracking_dat = load(file_info.trackingFile, file_info.tracksFinal);
            detection_res = detection_dat.(file_info.detectionDat);
            tracksFinal = tracking_dat.(file_info.tracksFinal);
        end
                                  
        % convert u-track datastructure to new datastructure
        trajectories = uTrack2mat(tracksFinal, detection_res, ...
                                   pixel_size, n_before_frozen);
        
        if params.runTrackSingle.smooth_trajectories
            trajectories = smooth_traj(trajectories);
        end
        
        tracks = cell(length(trajectories),1);

        % transform_data_structure also does many calculations.
        parfor i=1:length(trajectories)
            tracks{i} = transformDataStructure(trajectories{i}, file_info);
        end
        
        tracks(cellfun('isempty',tracks))=[];
        tracks = horzcat(tracks{:});
        tracks([tracks.trajtime] == 0) = [];

        save(file_info.swimtrackerFile, 'tracks', 'params');

        runTumbleDetectionSingle(movie_name, params);

        fprintf('done\n');
    catch e
        fprintf('%s\n', getReport(e, 'extended'));
    end
end

function skip = tooManyObjects(detect_res, max_obj)
    mean_n_obj = mean(cell2mat(arrayfun(@(x) size(x.xCoord(:,1),1), ...
                                        [detect_res], 'uni', 0)));
    if mean_n_obj > max_obj
        fprintf('mean number of objects over threshold\n');
        skip = 1;
    else
        skip = 0;
    end
end
