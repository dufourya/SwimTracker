function runTumbleDetectionSingle(data, params)
    % Categorizes each segment of a track into one of three states: 'run',
    % 'tumble', or 'transition' (accelerating from a run to a tumble, or
    % vice versa), based on the variables speed, acceleration, and change
    % in angle (theta).
    %
    % First, the speed and acceleration of each track are normalized to the
    % mean run speed (MRS) of that track.  Then, the data are fit using a model
    % with three (one for each state) tri-variate (one for each variable)
    % gaussians.  Finally, clustering is used to assign a state to each
    % frame of each track.
    %
    % The values for speed, acceleration, and theta are re-calculated and
    % used in the next iteration of fitting/clustering until 'tol' falls
    % below a specified value (see below) or until 'max_iter' iterations
    % have been performed.
    %
    % Convergence:  'tol' is calculated as 100*(MRS_prev - MRS_cur)/MRS_prev
    %
    % Initially, the entire track is used to calculate MRS.  After
    % one iteration, parts of the track will be assigned to other (slower)
    % states, so the MRS should be positive.  In subsequent iterations, it may
    % be positive or negative, but should always be decreasing.

    close all
    RELEASE = version('-release');
    startMatlabPool();

    filter = params.filter;
    tol              = params.runTumbleDetectionSingle.tol;
    max_iter         = params.runTumbleDetectionSingle.max_iter;
    show_diagnostics = params.runTumbleDetectionSingle.show_diagnostics;
    learn_model      = params.runTumbleDetectionSingle.learn_model;

    tumble_model = params.tumble_model.model;

    if isstruct(data)
        tracks = data;
        file_info = RadFileInfo(tracks(1).metadata.name);
    else
        file_info = RadFileInfo(data);
        fprintf('input file: %s\n', data);
        fprintf('loading data ''%s''\n', file_info.swimtrackerFile);
        data = load(file_info.swimtrackerFile,'tracks');
        if ~isempty(filter)
            data = filterData(data, filter);
        end
        tracks = data.tracks;
        fprintf('done.\n');
    end

    fprintf('Detecting tumbles...\n');
    
    iter = 1;
    nTracks = numel(tracks);
    parfor z = 1:nTracks
        tracks(z).meanrunspeed = prctile(tracks(z).speed,95);
    end

    if ~isfield(tracks,'ppState')
        tracks(1).ppState = [];
    end
    
    if ~isfield(tracks, 'diffcoeff')
        tracks(1).diffcoeff = [];
    end

    if show_diagnostics
        figure
    end

    n_states = 3;
    sum_of_states = (n_states*(n_states+1))/2;
    RUN = 0;
    TUMBLE = 1;
    INTER = 2;

    while iter <= max_iter

        tracks2 = tracks;

        norm_speed_all = [];
        theta_all = [];
        norm_accel_all = [];

        if learn_model
            % Calculate statistics for all of the tracks.
            fprintf('Calculating track statistics...\n');
            parfor i = 1:nTracks
                if ~isempty(tracks2(i).deltaangle(2:end-2)) && ...
                        ~isempty(tracks2(i).acceleration(2:end-2))
                    norm_speed_all = ...
                        [norm_speed_all; ...
                        tracks2(i).speed(2:end-2-sum(tracks2(i).frozen)) ./ ...
                        tracks2(i).meanrunspeed];
                    theta_all = [theta_all; tracks2(i).deltaangle(2:end-2)];
                    norm_accel_all = [norm_accel_all; ...
                        tracks2(i).acceleration(2:end-2) ./ ...
                        tracks2(i).meanrunspeed];
                end
            end

            X_all = [norm_speed_all theta_all norm_accel_all];

            options = statset('Display', 'off', 'MaxIter', 1000);

            % define start struct with guesses at initial conditions for fit:
            % row 1 = run, row 2 = intermediate, row 3 = tumble.
            % col 1 = speed, col 2 = theta, col 3 = accel.
            % mu = mean
            % Sigma = covariance of each compoent.
            % PComponents = frequency of each state.
            ns   = std(norm_speed_all);
            nacc = std(norm_accel_all);
            thet = std(theta_all);
                       
            s = struct( ...
                'mu', [1.0 0 0;
                       0.5 0 0; 
                       0.2 0 0], ...
                'Sigma', cat(3, ...
                             [0.1*ns  0        0;
                              0       0.1*thet 0;
                              0       0        0.1*nacc], ...
                             [0.5*ns  0        0;
                              0       0.1*thet 0;
                              0       0        nacc], ...
                             [0.1*ns  0        0;
                              0       thet     0;
                              0       0        0.1*nacc]), ...
                'PComponents',[0.65 0.10 0.25]); % run, intermediate, tumble

            fprintf('Generating model...\n');
            tumble_model = gmdistribution.fit(X_all,3,'Options',options, ...
                                             'Replicates',1,'Start',s, ...
                                             'Regularize',0.000001);
        end

        % In general, the rows of the final state assignments in the
        % matrix may not correspond to the initial row assignments, so we
        % deduce which row is which state from the data.

        %  The row with the lowest run speed ia the tumble state.
        [~,TUMBLE_ROW] = min(tumble_model.mu(:,1));

        % The row with the highest run speed corresponds to runs.
        [~,RUN_ROW] = max(tumble_model.mu(:,1));

        % The other row is the intermediate state.
        INTER_ROW = sum_of_states - TUMBLE_ROW - RUN_ROW;

        % Determine state of each segment in each track.
        fprintf('Assigning states...\n');
        parfor i=1:nTracks
            norm_speed = ...
                tracks2(i).speed(2:end-2-sum(tracks2(i).frozen)) ./ ...
                tracks2(i).meanrunspeed;
            theta      = tracks2(i).deltaangle(2:end-2);
            norm_accel = tracks2(i).acceleration(2:end-2) ./ ...
                         tracks2(i).meanrunspeed;

            X_track = [norm_speed theta norm_accel];

            if size(X_track,2)==3

                [z,~,P] = cluster(tumble_model, X_track);
                is_tumble = z == TUMBLE_ROW;
                is_inter = z == INTER_ROW;

                tracks2(i).tumbleModel = tumble_model;

                % Assign consistent values to states.  First and last two
                % entries of each trajectory are NaNined and assigned NaN.
                % Tumble state is TUMBLE based on truth of is_tumble.
                % Intermediate state is assigned INTER.  "The rest" belong to run
                % state and are left with RUN.

                tracks2(i).tumble = [NaN;is_tumble;NaN;NaN; ...
                                     nan(sum(tracks2(i).frozen),1)];
                tracks2(i).tumble([NaN;is_inter;NaN;NaN; ...
                                   nan(sum(tracks2(i).frozen),1)]>0) = INTER;
                tracks2(i).meanrunspeed = mean([tracks(i).meanrunspeed, ...
                                                mean(tracks2(i).speed(tracks2(i).tumble == RUN))]);
                                            
                % col 1 = run, col 2 = intermediate, col 3 = tumble
                tracks2(i).ppState = [[NaN NaN NaN];P;[NaN NaN NaN];[NaN NaN NaN]; ...
                                      nan(sum(tracks2(i).frozen),3)];

            else
                tracks2(i).tumbleModel  = tumble_model;
                tracks2(i).tumble       = nan(numel(tracks2(i).x),1);
                tracks2(i).meanrunspeed = mean(tracks2(i).speed);
                tracks2(i).ppState      = nan(numel(tracks2(i).x),3);
            end
        end
        
        % initial guess for mean run speed.
        for z = find(isnan([tracks2.meanrunspeed])|[tracks2.meanrunspeed]<0.1)
            tracks2(z).meanrunspeed = prctile(tracks2(z).speed,95);
        end

        if show_diagnostics
            hold on
            subplot(ceil(max_iter/2),floor(max_iter/2), iter);
            scatter([tracks.meanrunspeed], [tracks2.meanrunspeed]);
            title(iter);
            xlabel('previous meanrunspeed');
            ylabel('current meanrunspeed');
            hold off;
        end
        
        tol_calc = mean(abs([tracks.meanrunspeed] - [tracks2.meanrunspeed]) ./ ...
                        [tracks.meanrunspeed]) * 100;
        fprintf('iteration %d, %% change: %.3f\n', iter, tol_calc);
        
        % update data.
        tracks = tracks2;
        params.tumble_model.model = tumble_model;
        params.tumble_model.mu    = tumble_model.mu;
        params.tumble_model.sigma = tumble_model.Sigma;
        if strcmp(RELEASE, '2013b') || strcmp(RELEASE, '2014a')
            params.tumble_model.component_proportion = tumble_model.PComponents;
        elseif strcmp(RELEASE, '2014b') || strcmp(RELEASE, '2015a')
            params.tumble_model.component_proportion = tumble_model.ComponentProportion;
        end
        
        if tol_calc <= tol
            break;
        end
        iter = iter + 1;
    end % end while
    
    tracks = calculateTumbleBasedStatistics(tracks);
    fprintf('done.\n');

    fprintf('saving data...\n');
    file_name = [];
    save_version = '';
    if params.runTumbleDetectionSingle.learn_model == 1
        file_name = 'all_tracks';
        save_version = '-v7.3';
    else
        file_name = file_info.swimtrackerFile;
    end
    save(file_name, 'tracks', 'params', save_version);
    fprintf('done with tumble detection\n');
end
