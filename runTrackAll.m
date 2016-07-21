function runTrackAll(params)
    % Runs runTrackSingle on all appropriate files in the current folder.
    % These can be movies (*.bin), detection files (*_detection.mat),
    % tracking files (*_tracking.mat), or the final dataset results
    % (*_swimtracker.mat).

    startMatlabPool();
    
    extensions = sprintf('%s|%s|%s|%s', FileInfo.binExt, ...
                                        RadFileInfo.detectionExt, ...
                                        RadFileInfo.trackingExt, ...
                                        RadFileInfo.swimtrackerExt);
    files = dir();
    files = files(find(~cellfun(@isempty, regexp({files(:).name}, extensions))));
    names = cellfun(@(x) regexp(x, '^(.*?)\.', 'tokens'), {files(:).name}, 'uni', 0);
    unique_names = unique(cellfun(@(x) x{1,1}, names));
    
    if numel(files) == 0
        fprintf('No appropriate files found in this folder\n')
    end
    
    saveResults.dir = '.';
    for i = 1:numel(unique_names)
        try
            runTrackSingle(unique_names{i}, params)
        catch e
            fprintf('%s\n', getReport(e, 'extended'));
        end
     end
