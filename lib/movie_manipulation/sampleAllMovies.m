function sampleAllMovies(folder, start_time, end_time)
    if nargin <1 || isempty(folder) || strcmp(folder, '.')
        folder = pwd();
    end
    movies = dir(fullfile(folder, strcat('*', FileInfo.binExt)));
    for i=1:numel(movies)
        sampleMovie(movies(i).name, start_time, end_time);
    end
    fprintf('done\n');
end