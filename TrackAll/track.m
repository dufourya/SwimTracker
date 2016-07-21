function [utracks, params] = track(movieInfo, file_info, metadata, params)
    mean_cell_speed = params.utrack.mean_cell_speed;
    
    % minimum allowed search radius. The search radius is calculated on the spot in
    % the code given a feature's motion parameters. If it happens to be smaller
    % than this minimum, it will be increased to the minimum.
    params.utrack.costMatrices(1).parameters.minSearchRadius = mean_cell_speed * metadata.getImageIntervalMs()/1000 / metadata.getPixelSize();
    params.utrack.costMatrices(2).parameters.minSearchRadius = params.utrack.costMatrices(1).parameters.minSearchRadius;
    
    % maximum allowed search radius. Again, if a feature's calculated search
    % radius is larger than this maximum, it will be reduced to this maximum.
    params.utrack.costMatrices(1).parameters.maxSearchRadius = 2 * mean_cell_speed * metadata.getImageIntervalMs()/1000 / metadata.getPixelSize();
    params.utrack.costMatrices(2).parameters.maxSearchRadius = params.utrack.costMatrices(1).parameters.maxSearchRadius;
    
    % minimum track segment length to classify it as linear or random.
    params.utrack.costMatrices(2).parameters.lenForClassify = round(0.5 / metadata.getImageIntervalMs());
    
    costMatrices    = params.utrack.costMatrices;
    gapCloseParam   = params.utrack.gapCloseParam;
    kalmanFunctions = params.utrack.kalmanFunctions;
    probDim         = params.utrack.probDim;
    verbose         = params.utrack.verbose;
    
    saveResults.filename = file_info.trackingFile;
    
    [utracks,kalmanInfoLink,errFlag] = ...
        trackCloseGapsKalmanSparse(movieInfo, costMatrices, ...
                                   gapCloseParam, kalmanFunctions, ...
                                   probDim, saveResults, verbose);
end