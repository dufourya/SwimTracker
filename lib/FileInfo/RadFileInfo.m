classdef RadFileInfo < FileInfo
    properties (Constant = true)
        detectionExt      = '.rad_detection.mat';
        trackingExt       = '.rad_tracking.mat';
        swimtrackerExt    = '.rad_swimtracker.mat';
        detectionDat      = 'rad_detection';
    end
    methods
        function f = RadFileInfo(file)
            f = f@FileInfo(file);
        end
    end
end
