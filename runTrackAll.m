function runTrackAll()
%     startMatlabPool();

extensions = {FileInfo.metaExt,...
    FileInfo.binExt, ...
    RadFileInfo.detectionExt, ...
    RadFileInfo.trackingExt, ...
    RadFileInfo.swimtrackerExt};

files = dir('*.*');
files = files(~cellfun(@isempty, regexp({files.name}','^\w')));
files = files(~contains({files(:).name},'.pictures.meta'));
files = files(contains({files(:).name}, extensions));
names = replace({files(:).name}',extensions,'');
unique_names = unique(names);

if numel(files) == 0
    fprintf('No appropriate files found in this folder\n')
end

%     saveResults.dir = '.';
for i = 1:numel(unique_names)
    try
        runTrackSingle(unique_names{i}, 0)
    catch e
        fprintf('%s\n', getReport(e, 'extended'));
    end
end
