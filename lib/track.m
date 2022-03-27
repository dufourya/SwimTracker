function [utracks, params] = track(movieInfo, file_info, metadata, params)

mean_cell_speed = params.radcenter.mean_cell_speed;

cell_speed_error = metadata.getResolution + 2*sqrt(params.radcenter.cell_area/pi) + sqrt(10*metadata.getImageIntervalMs()/1000);

params.utrack.gapCloseParam.timeWindow = ...
    ceil(1000 * params.utrack.gapCloseParam.timeWindow / metadata.getImageIntervalMs);
params.utrack.costMatrices(1).parameters.nnWindow = ...
    ceil(1000 * params.utrack.costMatrices(1).parameters.nnWindow / metadata.getImageIntervalMs);

params.utrack.costMatrices(2).parameters.timeReachConfB = params.utrack.gapCloseParam.timeWindow;
params.utrack.costMatrices(2).parameters.timeReachConfL = params.utrack.gapCloseParam.timeWindow;


params.utrack.costMatrices(2).parameters.resLimit = ...
    metadata.getResolution / metadata.getPixelSize();

% minimum allowed search radius. The search radius is calculated on the spot in
% the code given a feature's motion parameters. If it happens to be smaller
% than this minimum, it will be increased to the minimum.
params.utrack.costMatrices(1).parameters.minSearchRadius = ...
    1/2 * mean_cell_speed * metadata.getImageIntervalMs() / 1000 / metadata.getPixelSize();%mean_cell_speed /3 * metadata.getImageIntervalMs()/1000 / metadata.getPixelSize();

params.utrack.costMatrices(2).parameters.minSearchRadius = ...
    params.utrack.costMatrices(1).parameters.minSearchRadius;

% maximum allowed search radius. Again, if a feature's calculated search
% radius is larger than this maximum, it will be reduced to this maximum.
params.utrack.costMatrices(1).parameters.maxSearchRadius = ...
    2 * ((mean_cell_speed * metadata.getImageIntervalMs() / 1000) + cell_speed_error) / metadata.getPixelSize();

params.utrack.costMatrices(2).parameters.maxSearchRadius = ...
    params.utrack.costMatrices(1).parameters.maxSearchRadius;

% minimum track segment length to classify it as linear or random.
params.utrack.costMatrices(2).parameters.lenForClassify = ceil(1 / metadata.getImageIntervalMs() * 1000);

params.utrack.costMatrices(2).parameters.brownStdMult = ...
    params.utrack.costMatrices(1).parameters.brownStdMult*ones(params.utrack.gapCloseParam.timeWindow,1);

params.utrack.costMatrices(2).parameters.linStdMult = ...
    params.utrack.costMatrices(2).parameters.linStdMult*ones(params.utrack.gapCloseParam.timeWindow,1);


% params.utrack.costMatrices(1).parameters.kalmanInitParam.convergePoint = [];
% params.utrack.costMatrices(1).parameters.kalmanInitParam.searchRadiusFirstIteration = ceil(mean_cell_speed * metadata.getImageIntervalMs() / 1000 / metadata.getPixelSize());
% params.utrack.costMatrices(1).parameters.kalmanInitParam.initVelocity = 0;

costMatrices    = params.utrack.costMatrices;
gapCloseParam   = params.utrack.gapCloseParam;
kalmanFunctions = params.utrack.kalmanFunctions;
probDim         = params.utrack.probDim;
verbose         = params.utrack.verbose;

saveResults.filename = file_info.trackingFile;

[utracks,kalmanInfoLink,errFlag] = ...
    trackCloseGapsKalmanSparse(movieInfo, costMatrices, ...
    gapCloseParam, kalmanFunctions, probDim, saveResults, verbose);
end