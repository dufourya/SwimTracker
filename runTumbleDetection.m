function runTumbleDetection(params)
    % runs runTumbleDetectionSingle on all files in the current folder.

    close all

    startMatlabPool();

    filelisting = dir(fullfile(strcat('*', RadFileInfo.swimtrackerExt)));
    
    for j = 1:numel(filelisting)
        file_name = filelisting(j).name;
        if regexp(file_name, '^\.');
            continue
        end
        fprintf('   file: ''%s''\n', file_name);
        runTumbleDetectionSingle(file_name, params);
    end

    clear all
    fprintf('all done\n');
end
