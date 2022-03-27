function runTrackSingle(movie_name,hpcc)

fprintf('* %s\n',movie_name);

warning('off','MATLAB:DELETE:FileNotFound');

if nargin < 2 || isempty(hpcc)
    hpcc = 1;
end

params           = parameters(movie_name,hpcc);
mem_usage        = params.runTrackSingle.mem_usage;
detection_method = params.runTrackSingle.detection_method;
max_object_count = params.runTrackSingle.maximum_object_count;
psfSigma         = params.utrack.psfSigma;
alphaLocMax      = params.utrack.alphaLocMax;
doMMF            = params.utrack.doMMF;
testAlpha        = params.utrack.testAlpha;

if mem_usage > 1
    error('Can''t request more than the available amount of memory.')
end

detectionParam = [];
is_utrack = true;

if (strcmp(detection_method, 'gaussian'))
    detectionParam.psfSigma = psfSigma;
    detectionParam.alphaLocMax = alphaLocMax;
    detectionParam.doMMF = doMMF;
    detectionParam.testAlpha = testAlpha;
elseif (strcmp(detection_method, 'radcenter'))
    is_utrack = false;
else
    error_string = ['Don''t recognize detection method ''' ...
        detection_method '''.\nOptions are ''gaussian'' or ''radcenter''.\n'];
    error(error_string);
end

saveResults.dir = '.';
if is_utrack
    file_info = UtrackFileInfo(movie_name);
else
    file_info = RadFileInfo(movie_name);
end

if ~exist(file_info.metaFile,'file')
    fprintf('Missing metadata file, skipping.\n');
    return
else
    metadata = Metadata(file_info.metaFile);
    metadata.read();
    pixel_size = metadata.getPixelSize();
    n_images = metadata.getNumberOfImages();
    n_before_frozen = metadata.getNImagesBeforeFreeze();
    if isempty(n_before_frozen) || n_before_frozen > n_images
        n_before_frozen = n_images;
    end
end

%% Detection file
if ~exist(file_info.detectionFile,'file') && exist(file_info.binFile,'file')
    delete(file_info.trackingFile, file_info.swimtrackerFile, ...
        file_info.tableFile, file_info.pngFile, file_info.svgFile, ...
        file_info.diagpngFile, file_info.diagsvgFile);
    if is_utrack
        movieParam.imageDir = '.';
        movieParam.imageSizeY = metadata.getImageHeightPixel();
        movieParam.imageSizeX = metadata.getImageWidthPixel();
        movieParam.imageName = file_info.name;
        movieParam.imageFile = strcat(movieParam.imageDir, filesep, ...
            file_info.binFile);
        movieParam.metaFile = strcat(movieParam.imageDir, filesep, ...
            file_info.metaFile);
        saveResults.filename = file_info.detectionFile;
        [detection_res,~,~,~,~] = ...
            detectSubResFeatures2D_StandAlone_Yann(movieParam, ...
            detectionParam, ...
            saveResults);
    else
        detection_res = radDetect(file_info, metadata, params);
    end
elseif exist(file_info.detectionFile,'file') && ~exist(file_info.trackingFile,'file')
    fprintf('Loading detection file.\n');
    dat = load(file_info.detectionFile, file_info.detectionDat);
    detection_res = dat.(file_info.detectionDat);
end

%% Tracking file
if exist('detection_res','var')
    delete(file_info.swimtrackerFile, ...
        file_info.tableFile, file_info.pngFile, file_info.svgFile, ...
        file_info.diagpngFile, file_info.diagsvgFile);
    if tooManyObjects(detection_res, max_object_count)
        return;
    end
    fprintf('Tracking objects...\n');
    [tracksFinal, params] = track(detection_res, file_info, metadata, params);
elseif exist(file_info.trackingFile,'file') && ~exist(file_info.swimtrackerFile,'file')
    fprintf('Loading tracking file.\n');
    %     dat = load(file_info.detectionFile, file_info.detectionDat);
    %     detection_res = dat.(file_info.detectionDat);
    tracking_dat = load(file_info.trackingFile, file_info.tracksFinal);
    tracksFinal = tracking_dat.(file_info.tracksFinal);
end

%% Swimtracker file
if exist('tracksFinal','var')
    delete(file_info.tableFile, file_info.pngFile, file_info.svgFile, ...
        file_info.diagpngFile, file_info.diagsvgFile);
    %     tracksFinal(cellfun(@numel,{tracksFinal.tracksFeatIndxCG})<4)=[];
    if isempty(tracksFinal)
        fprintf('0 track to analyze\n');
        return;
    end
    t = metadata.getImageIntervalMs() ./ 1000 .* (0:(n_images-1));
    trajectories = uTrack2Yann(tracksFinal, t, ...
        pixel_size, n_before_frozen);
    %     fprintf('Calculating movie drift...\n');
    [trajectories, movie_drift] = correct_drift_traj(trajectories,params.runTrackSingle.correct_drift);
    fprintf('%f um/s\n',movie_drift.total);
    if params.runTrackSingle.smooth_trajectories
        trajectories = smooth_traj(trajectories);
    end
    tracks = cell(length(trajectories),1);
    fprintf('Calculating trajectory statistics...');
    if hpcc
        for i = 1:length(trajectories)
            tracks{i} = transformDataStructure(trajectories{i}, file_info, i);
        end
    else
        parfor i = 1:length(trajectories)
            tracks{i} = transformDataStructure(trajectories{i}, file_info, i);
        end
    end
    
    %     tracks(cellfun('isempty',tracks))=[];
    tracks = horzcat(tracks{:});
    tracks([tracks.trajtime] == 0) = [];
    if params.runTrackSingle.filter_tracks
        window    = params.filterAnomalousTracks.window;
        msd_rep   = params.filterAnomalousTracks.msd_rep;
        fdr       = params.filterAnomalousTracks.fdr;
        diag_plot = params.filterAnomalousTracks.diag_plot;
        tracks = filterAnomalousTracks(tracks, window, msd_rep, ...
            fdr, diag_plot);
    end
    fprintf('done.\n');
    
    if ~isempty(params.tumble_model)
        tracks = runTumbleDetectionSingle_hmm(tracks, params, 0);
    end
    fprintf('Saving file...');
    save(file_info.swimtrackerFile, 'tracks', 'params', 'movie_drift');
    fprintf('done.\n');
    
elseif exist(file_info.swimtrackerFile,'file') && ~exist(file_info.tableFile,'file')
    fprintf('Loading swimtracker file.\n');
    tracking_dat = load(file_info.swimtrackerFile, file_info.tracksDat);
    tracks = tracking_dat.tracks;
end

if ~exist(file_info.pngFile,'file')
    plotSwimTracker(file_info.name,1);
end

if ~exist(file_info.diagpngFile,'file')
    plotTrackDiagnostics(file_info.name,1);
end

if exist(file_info.picturesFile,'file') && ~exist(file_info.FASTfile,'file')
    FASTanalysis(file_info.name,'hpcc');
end

% Table file
if exist('tracks','var')
    fprintf('Writing table file.\n');
    tracksTable = reduceTrackStruct(tracks);
    writetable(tracksTable,strcat(tracksTable.metadata{1},'.table.txt'));
end

warning('on','all');
end

function skip = tooManyObjects(detect_res, max_obj)
mean_n_obj = mean(cell2mat(arrayfun(@(x) size(x.xCoord(:,1),1), ...
    detect_res, 'uni', 0)));
if isnan(mean_n_obj)
    mean_n_obj = 0;
end

if mean_n_obj > max_obj || mean_n_obj < 1
    fprintf('%d detected objects on average\nToo many or too few objects to track! (%d cells max)\n',[round(mean_n_obj),max_obj]);
    skip = 1;
else
    fprintf('%d detected objects on average\n',round(mean_n_obj));
    skip = 0;
end
end
