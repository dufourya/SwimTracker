function params = parameters()
    %% Filter used by any of the tracking components
    params.filter = [];

    %% runTrackSingle
    
    % the fraction of the computer's available memory
    % to use for processing images. This determines
    % the number of images to process at one time.
    params.runTrackSingle.mem_usage = 0.7;
    
    % the amount of time in seconds to use for
    % background subtraction. The movie will be
    % broken into segments corresponding to
    % 'background_sample_time' seconds. This set
    % of images are then averaged to produce the
    % background.
    params.runTrackSingle.background_sample_time = 6;
    
    % the method to use for detection.  Can be
    % 'gaussian' to use 2D gaussian fitting from the
    % utrack package or 'radcenter' to use the radial
    % symmetry method of Parthasarathy (2012).
    % 'radcenter' is ~30x faster.
    params.runTrackSingle.detection_method = 'radcenter';
    
    params.runTrackSingle.smooth_trajectories = 1;
        
    % Sets the maximum number of objects that will be considered for
    % tracking.
    params.runTrackSingle.maximum_object_count = 5e3;
    
    %% radcenter detection

    % If true, opens the gui associated with theradcenter code. Useful for
    % determining radcenter parameter values.

    params.radcenter.run_gui = 0;
    
    % Should the background be subtracted?
    params.radcenter.subtract_background = 1;
    
    % the number of standard deviations above
    % background to be considered a particle.
    params.radcenter.detection_threshold = -6;
    
    % in um^2, 1um * 0.5um * Pi for a typical E. coli. Phase
    % contrast makes this larger (~4 um^2).
    params.radcenter.cell_area = 4; % um^2
    
    params.radcenter.bpfiltsize = [];
    params.radcenter.nsize = [];
    
    
    %% utrack tracking

    % The mean cell speed.
    params.utrack.mean_cell_speed = 30; % um/sec

    % the following explanations are paraphrased from trackCloseGapsKalmanSparse.m
    
    % Names of Kalman filter functions for self-adaptive tracking.
    % For non-self-adaptive tracking, enter [].
    
    % Reserves memory for kalmanFilterInfo.
    params.utrack.kalmanFunctions.reserveMem = 'kalmanResMemLM';
    
    % Initializes the Kalman filter for an appearing feature.
    params.utrack.kalmanFunctions.initialize = 'kalmanInitLinearMotion';
    
    % Calculates the Kalman gain after linking.
    params.utrack.kalmanFunctions.calcGain = 'kalmanGainLinearMotion';
    
    % Reverses time (and associated variables) in kalmanInfoLink between the
    % different frame-to-frame linking steps.
    params.utrack.kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
    
    % maximum allowed time gap (in frames) between a track segment end and a track
    % segment start that allows linking them.
    params.utrack.gapCloseParam.timeWindow = 6;
    
    % 1 if merging and splitting are to be considered, 2 if only merging is to be
    % considered, 3 if only splitting is to be considered, 0 if no merging or
    % splitting are to be considered.
    params.utrack.gapCloseParam.mergeSplit = 0;
    
    % minimum length of track segments from linking to be used in gap closing.
    params.utrack.gapCloseParam.minTrackLen = 2;
    
    % 1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
    params.utrack.gapCloseParam.diagnostics = 0;
    
    % u-track 2.1
    %costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';
    %costMatrices(1).parameters.linearMotion = 2;
    
    % current u-track
    params.utrack.costMatrices(1).funcName = 'costMatLinearMotionLink2';
    
    % use linear motion Kalman filter.
    params.utrack.costMatrices(1).parameters.linearMotion = 1;
    
    
    % multiplication factor to calculate search radius from standard deviation.
    params.utrack.costMatrices(1).parameters.brownStdMult = 3;
    
    % 1 if you want to expand the search radius of isolated features in the linking
    % (initial tracking) step.
    params.utrack.costMatrices(1).parameters.useLocalDensity = 0;
    
    % number of frames before the current one where you want to look to see a
    % feature's nearest neighbor in order to decide how isolated it is (in the
    % initial linking step).
    params.utrack.costMatrices(1).parameters.nnWindow = params.utrack.gapCloseParam.timeWindow;
    
    % if you want to plot the histogram of linking distances up to certain frames,
    % indicate their numbers; 0 or empty otherwise. Does not work for the first or
    % last frame of a movie.
    params.utrack.costMatrices(1).parameters.diagnostics = [];
    
    % u-track 2.1
    %costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';
    
    % current u-track
    params.utrack.costMatrices(2).funcName = 'costMatLinearMotionCloseGaps2';
    params.utrack.costMatrices(2).parameters = params.utrack.costMatrices(1).parameters;
    params.utrack.costMatrices(2).parameters.brownStdMult = ...
    params.utrack.costMatrices(1).parameters.brownStdMult*ones(params.utrack.gapCloseParam.timeWindow,1);
    
    % power for scaling the Brownian search radius with time, before and after
    % timeReachConfB (next parameter).
    params.utrack.costMatrices(2).parameters.brownScaling = [0.5 0.01];
    
    % before timeReachConfB, the search radius grows with time with the power in
    % brownScaling(1); after timeReachConfB it grows with the power in
    % brownScaling(2).
    params.utrack.costMatrices(2).parameters.timeReachConfB = 4;
    
    % for merging and splitting. Minimum and maximum ratios between the intensity
    % of a feature after merging/before splitting and the sum of the intensities of
    % the 2 features that merge/split.
    params.utrack.costMatrices(2).parameters.ampRatioLimit = [0.7 4];
    
    % multiplication factor to calculate linear search radius from standard
    % deviation.
    params.utrack.costMatrices(2).parameters.linStdMult = 3*ones(params.utrack.gapCloseParam.timeWindow,1);
    
    % power for scaling the linear search radius with time (similar to
    % brownScaling).
    params.utrack.costMatrices(2).parameters.linScaling = [0.5 0.01];
    
    % similar to timeReachConfB, but for the linear part of the motion.
    params.utrack.costMatrices(2).parameters.timeReachConfL = params.utrack.gapCloseParam.timeWindow;
    
    % maximum angle between the directions of motion of two tracks that allows
    % linking them (and thus closing a gap). Think of it as the equivalent of a
    % searchRadius but for angles.
    % costMatrices(2).parameters.maxAngleVV = 30;
    params.utrack.costMatrices(2).parameters.maxAngleVV = 90;
    
    % penalty for increasing temporary disappearance time (disappearing for n
    % frames gets a penalty of gapPenalty^n).
    params.utrack.costMatrices(2).parameters.gapPenalty = 1.5;
    
    % resolution limit, which is generally equal to 3 * point spread function
    % sigma.
    %costMatrices(2).parameters.resLimit = 3 * detectionParam.psfSigma;
    params.utrack.costMatrices(2).parameters.resLimit = 3.4;
    
    % print output during tracking
    params.utrack.verbose = 1;
    
    % Problem dimensionality. 2 (for 2D) or 3 (for 3D).
    % Optional. If not input, dimensionality will be derived from movieInfo.
    params.utrack.probDim = 2;
    
    %% runTumbleDetectionSingle
    
    % If 1, model parameters are generated at each iteration, which slows
    % down computation. If 0, uses model defined below based on RP437.
    % May be necessary if cells are very different from RP437.
    params.runTumbleDetectionSingle.learn_model = 0;
    
    % the minimum acceptable percent difference between the runspeed of the
    % previous iteration and the current one.
    params.runTumbleDetectionSingle.tol = 0.01;
    
    % the maximum number of iterations to perform.
    params.runTumbleDetectionSingle.max_iter = 20;
    
    % show some diagnostic plots of the whole dataset.
    params.runTumbleDetectionSingle.show_diagnostics = 0;
    
    % Component means
    %    speed               delta theta            acceleration
    params.tumble_model.mu = ...
        [0.991878498881794,  -0.000174698179524619, -0.00770215381116974; % run
         0.704904005311772,  -0.000795217497216972, -0.0309194682008243;  % intermediate
         0.324364764701188,   0.00161684957539444,   0.0625190586216278]; % tumble
    
    % Component covariances
    params.tumble_model.sigma = zeros(3,3,3);
    
    params.tumble_model.sigma(:,:,1) = ...
        [ 0.0218    0.0000   -0.0044 % run
          0.0000    0.0094   -0.0000 % tumble
         -0.0044   -0.0000    0.0090 % intermediate
        ];
    
    params.tumble_model.sigma(:,:,2) = ...
        [ 0.0554   -0.0001   -0.0192
         -0.0001    0.0920    0.0003
         -0.0192    0.0003    0.0517
        ];
    
    params.tumble_model.sigma(:,:,3) = ...
        [ 0.0313    0.0002   -0.0146
          0.0002    2.2841   -0.0006
         -0.0146   -0.0006    0.0544
        ];
    
    % The mixing proportion of each component.
    %    run                intermediate       tumble
    params.tumble_model.component_proportion = ...
        [0.334785512263232, 0.417120796327256, 0.248093691409512];
    
    % Generate the model.
    params.tumble_model.model = gmdistribution(params.tumble_model.mu, ...
                                               params.tumble_model.sigma, ...
                                               params.tumble_model.component_proportion);
    
    %% runTumbleDetectionLearn
    params.runTumbleDetectionLearn.show_diagnostics = 0;
end
