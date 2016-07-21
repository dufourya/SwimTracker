function viewMovie(image_file, data_path, movie_mag, subtract_background, sample_time, ...
                   show_background, show_tracking, show_trails, save_movie, ...
                   track_method, filter)
    % View a '.bin' movie with or without tracks.
    %
    % image_file         : the file to view (e.g., 'my_movie.bin').
    % data_path          : the path to the tracking data (if desired).
    % movie_mag          : the scaling size to view the movie, in percent.
    %                      (e.g., 100 for 100%). [50].
    % subtract_background: should the background be subtracted? [0].
    % sample_time        : the sample time used to generate each background
    %                      image, in seconds. [25]
    % show_background    : displays background in a separate window, regardless
    %                      of 'subtract_background' option. [0].
    % show_tracking      : If tracking data is present, plots circles
    %                      around each tracked object. [0].
    % show_trails        : If tracking data is present, draws a 'trail'
    %                      behind each tracked object.
    %                      This is currently very slow. [0].
    % save_movie         : Save the movie as an avi file, without
    %                      displaying it. [0].
    % track_method       : The method used to track the cells.  Can be
    %                      'radcenter'.
    % filter             : a structure of the type used by filterData.m.
        
    if nargin < 11 || isempty(filter)
        filter = [];
    end
    if nargin < 10 || isempty(track_method)
        track_method = 'radcenter';
    end
    if nargin < 9 || isempty(save_movie)
        save_movie = 0;
    end
    if nargin < 8 || isempty(show_trails)
        show_trails = 0;
    end
    if nargin < 7 || isempty(show_tracking)
        show_tracking = 0;
    end
    if nargin < 6 || isempty(show_background)
        show_background = 0;
    end
    if nargin < 5 || isempty(show_background)
        sample_time = 25;
    end
    if nargin < 4 || isempty(subtract_background)
        subtract_background = 1;
    end
    if nargin < 3 || isempty(movie_mag)
        movie_mag = 50;
    end
    
    [path, name, ext] = fileparts(image_file);
    if isempty(path)
        path = pwd();
    end
    
    image_info = [];
    if strcmp(track_method, 'radcenter')
        image_info = RadFileInfo(image_file);
    else
        error(['do not recognize tracking method ''' track_method '''.'])
    end
    
    if ~exist(image_file, 'file')
        error(['could not find ''' image_file ''' on path.'])
    end
    
    metadata = Metadata(fullfile(data_path,image_info.metaFile));
    metadata.read();
    n_frames = metadata.getNumberOfImages();
    interval = metadata.getImageIntervalMs() / 1000;
    pixel_size = metadata.getPixelSize();
    height = metadata.getImageHeightPixel();
    width = metadata.getImageWidthPixel();
    pixel_type = metadata.getPixelType();
    
    if save_movie
        vid = VideoWriter(image_info.aviFile, 'Uncompressed AVI');
        open(vid);
    end
    
    colors = [];
    tracks = [];
    data = [];
    if show_tracking
        try
            data = load(fullfile(data_path, strcat(image_info.name, image_info.swimtrackerExt)));
            if ~isempty(filter)
                filtered_data = filterData(data, filter);
                if isempty(filtered_data)
                    return
                end
                tracks = filtered_data.tracks;
            else
                tracks = data.tracks;
            end
            colors = lines(length(tracks));
        catch e
            fprintf('%s\n', getReport(e, 'extended'));
            fprintf('WARNING: No tracking data to display\n');
        end

    end
    
    fid = fopen(image_file);
    
    sample_time = min(sample_time, interval*n_frames);
    block_size = floor(sample_time/interval);
    
    if subtract_background
        fprintf('Calculating background...');
        backgrounds = calculateBackground(image_file, metadata, block_size);
        fprintf('done.\n');
    end
    
    if show_background
        if save_movie
            background_fig = figure('Visible', 'off');
        else
            background_fig = figure;
        end
        
        bg_img = imshow(mean(backgrounds,3), [], 'InitialMagnification', ...
            movie_mag);
        
        if save_movie
            imwrite(uint16(mean(backgrounds,3)), image_info.backgroundFile);
        end
    end
    
    
    movie_fig = [];
    if save_movie
        movie_fig = figure('Visible','off');
    else
        movie_fig = figure();
    end
    
    iptsetpref('ImshowBorder','tight');

    image_name = strrep(name, '_', '\_');

    find_tol = 0.0001;
    times = zeros([1 n_frames]);
    scaling = [];
    frewind(fid);
    for i=1:n_frames
        %fprintf('frame %d\n', i);
        
        if save_movie
            fprintf('saving frame %d\n', i);
        end
        times(i) = (i-1)*interval;
        
        frame = fread(fid,width*height,pixel_type);
        if ~isempty(frame)
            frame = flipud(rot90(reshape(frame,[width height])));
            d_frame = double(frame);
            img = d_frame./mean(d_frame(:));
            if subtract_background
                scaling = [-0.2, 0.05];
                img = img - backgrounds(:,:,ceil(i/block_size));
            end
            set(0,'CurrentFigure',movie_fig);
            
            imshow(img, scaling, 'InitialMagnification', movie_mag);

            text(700, 150, sprintf('%s\nframe: %d\ntime: %.2f s\n', ...
                                  image_name, i, times(i)), ...
                                  'FontSize', 14);
            
            % draw trails
            if show_tracking
                hold on
                for j=1:length(tracks)
                    track_col = colors(j,:);
                    
                    % draw circles around cells
                    row = find(abs(tracks(j).time-times(i)) < find_tol);
                    if ~isempty(row)
                        x = tracks(j).x(row)./pixel_size;
                        y = tracks(j).y(row)./pixel_size;
                        plot(x, y, 'ro', 'color', track_col, 'MarkerSize', ...
                             10, 'LineWidth', 3);
                    end
                    
                    if show_trails
                        rows = zeros([1 i]);
                        % loop through times up to this frame
                        for k = 1:i
                            row = find(abs(tracks(j).time-times(k)) < ...
                                       find_tol);
                            if ~isempty(row)
                                rows(k) = row;
                            end
                        end
                        rows(rows==0) = [];
                        if ~isempty(rows)
                            x = tracks(j).x(rows)./pixel_size;
                            y = tracks(j).y(rows)./pixel_size;
                            plot(x, y, 'r-', 'color', track_col, ...
                                 'MarkerSize', 10, 'LineWidth', 3);
                        end
                    end
                end
            end
            drawnow();
            if save_movie
                writeVideo(vid, getframe(movie_fig));
            end
            hold off
        end
    end
    if save_movie
        close(vid);
    end
    fclose(fid);
end