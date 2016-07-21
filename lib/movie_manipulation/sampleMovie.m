function sampleMovie(movie, start_time, end_time)
    out_folder = 'sampled';
    if ~exist(out_folder, 'dir')
        mkdir(out_folder);
    end
    
    movie_info = FileInfo(movie);
    metadata = Metadata(movie_info.metaFile);
    metadata.read();
    height_px = metadata.getImageHeightPixel();
    width_px = metadata.getImageWidthPixel();
    pixel_type = metadata.getPixelType();

    frames = 1:metadata.getNumberOfImages();
    interval_sec = metadata.getImageIntervalMs() / 1000;
    times = (0:(numel(frames)-1))*interval_sec;
        
    sample_frames = frames(times>= start_time & times<=end_time);
    n_frames = numel(sample_frames);
    read_intervals = diff([0 sample_frames]);
    
    out_name = fullfile(out_folder, movie_info.name);
    metadata.rename(strcat(out_name, FileInfo.metaExt));
    metadata.set('SequenceAcquisition', 'NumberImages', n_frames);

    fprintf('Sampling %g to %g sec of ''%s''\n', start_time, end_time, ...
            movie);
    fin = fopen(movie);
    fout = fopen(strcat(out_name, FileInfo.binExt), 'w');
    for i=1:numel(sample_frames)
        frame = fread(fin, width_px*height_px*read_intervals(i), pixel_type);
        frame = flipud(rot90(reshape(frame, [width_px height_px])));
        fwrite(fout, frame, pixel_type);

    end
    metadata.write();
    fclose('all');
end