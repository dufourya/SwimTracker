function runTumbleDetectionLearn(params)
    % Concatenates all tracks below current folder and uses the combined 
    % data to learn the tumble bias model.
    
    close all
    startMatlabPool();
    
    filter = params.filter;
    
    filelisting = rdir(fullfile('.', '**', strcat('*', RadFileInfo.swimtrackerExt)));
    
    all_tracks = [];
    fprintf('Concatenating all tracks\n');
    parfor i = 1:numel(filelisting)
        fprintf('   loading file: ''%s''\n', filelisting(i).name);
        data = load(filelisting(i).name, 'tracks');
        if ~isempty(filter)
            data = filterData(data, filter);
        end
        if ~isempty(data)
            if ~isfield(data.tracks(1), 'diffcoeff')
                data.tracks(1).diffcoeff = [];
            end
            all_tracks = [all_tracks data.tracks];
        end
    end
    save('all_tracks', 'all_tracks', 'params', '-v7.3');
    
    if isempty(all_tracks)
        error('No swimtracker files found!');
    end
    
    if params.runTumbleDetectionLearn.show_diagnostics
        figure
        hold on
        s(1) = subplot(2,2,1);
        hist(s(1), log10([all_tracks.fit_runtime]),100);
        xlabel(s(1), 'log10[fit\_runtime]');
        
        s(2) = subplot(2,2,2);
        hist(s(2), log10([all_tracks.fit_diffcoeff_nocurve]),100);
        xlabel(s(2), 'log10[fit\_diffcoeff\_nocurve]');
        
        s(3) = subplot(2,2,3);
        hist(s(3), ([all_tracks.meanspeed]),100);
        xlabel(s(3), 'meanspeed');
        
        s(4) = subplot(2,2,4);
        scatter(s(4), log10([all_tracks.fit_angle]), ...
                [all_tracks.meanspeed],'.');
        xlabel(s(4), 'log10[fit\_runtime]');
        ylabel(s(4), 'meanspeed');
        
        mtit('cumulative statistics');
        
        drawnow;
    end
    
    clear filelisting i idx j tracks tracks0 tracks1;
    
    params.runTumbleDetectionSingle.learn_model = 1;
    runTumbleDetectionSingle(all_tracks, params);
end
