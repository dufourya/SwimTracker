function rad_detection = radDetect(file_info, metadata, params)
    background_sample_time = params.runTrackSingle.background_sample_time;
    frac_mem_usage = params.runTrackSingle.mem_usage;
    thresh = params.radcenter.detection_threshold;
    cell_area = params.radcenter.cell_area;
    bpfiltsize = params.radcenter.bpfiltsize;
    nsize = params.radcenter.nsize;
    
    is_fluorescent = metadata.isFluorescent();
    
    subtract_background = params.radcenter.subtract_background;
    run_gui = params.radcenter.run_gui;
    f_name = getFunctionName();
    
    width = metadata.getImageWidthPixel();
    height = metadata.getImageHeightPixel();
    
    pixel_type = metadata.getPixelType();
    interval = metadata.getImageIntervalMs() / 1000;
    n_frames = metadata.getNumberOfImages();

    pixel_size = metadata.getPixelSize();
    bytes_per_pixel = metadata.getBytesPerPixel();
    image_size = width*height*bytes_per_pixel/1e6;
    
    if thresh > 0
        thresh = -thresh;
    end
    
    % fo4_rp options
    if isempty(bpfiltsize)
        bpfiltsize = round(cell_area * 3 / pixel_size); % cell area ~1.5 um^2 * 10 /pixelarea
    end
    if isempty(nsize)
        nsize = bpfiltsize;
    end
    fitstr = 'radial';
            
    backgrounds = [];
    background_block_size = [];
    n_backgrounds = [];
    if subtract_background
        background_sample_time = min(background_sample_time, interval*n_frames);
        background_block_size = floor(background_sample_time/interval);
        
        fprintf('Calculating background...');
        backgrounds = calculateBackground(file_info.binFile, metadata, ...
            background_block_size);
        fprintf('done.\n');
        n_backgrounds = size(backgrounds,3);
    end
    
    detected_objs = [];    
    fprintf('[%s] Detecting cells in frames...\n', f_name);
    
    freemem = determineMemory();
    % conversion to double (64-bit) increases memory requirement by a
    % factor of 4.
    mem_to_use_mb = freemem * frac_mem_usage / 3;
    block_size = floor(mem_to_use_mb / image_size);
    n_blocks = ceil(n_frames / block_size);
    %fprintf('using %d MB of memory (%d blocks)\n', mem_to_use_mb, n_blocks);
    
    tic
    fid = fopen(file_info.binFile);
    for i = 1:n_blocks
        fprintf('reading block %d of %d...', i, n_blocks);
        frames = fread(fid,width*height*block_size,pixel_type);
        frames = reshape(frames, height, width, []);
        fprintf('done.\n');
       
        fprintf('processing images...');
        if run_gui
            for j = 1:size(frames,3)
                frame_number = (i-1)*block_size+j;
                frame = double(flipud(rot90(frames(:,:,j))));
                frame = prepImage(frame, backgrounds, frame_number, ...
                                  background_block_size, ...
                                  subtract_background, is_fluorescent);

                TrackingGUI_rp(frame);
                pause;
            end
        else
            parfor j = 1:size(frames,3)
                frame_number = (i-1)*block_size+j;
                frame = double(flipud(rot90(frames(:,:,j))));
                frame = prepImage(frame, backgrounds, frame_number, ...
                                  background_block_size, ...
                                  subtract_background, is_fluorescent);
                
                tmpobj = fo4_rp(frame, [bpfiltsize nsize], thresh, fitstr);
                % set the frame (not done automatically by fo4_rp).
                tmpobj(5,:) = frame_number;
                detected_objs = [detected_objs tmpobj];
            end
        end
        fprintf('done.\n');
        fprintf('[%s] %d of %d\n', f_name, min(i*block_size, n_frames), ...
                n_frames);
    end
    toc
    fprintf('[%s] done.\n', f_name);
    
    % convert the output to the form required by utrack.
    rad_detection = rad2utrack(detected_objs, metadata);
    
    save(file_info.detectionFile, 'detected_objs', ...
         'rad_detection' , 'params');
end

function images = normalize(images)
    min_all = min(images(:));
    max_all = max(images(:));
    images = (images - min_all) ./ (max_all-min_all);
end

function playMovie(frames)
    for i = 1:size(frames,3)
        imshow(frames(:,:,i), 'InitialMagnification', 67);
        drawnow
    end
end

function frame = prepImage(frame, backgrounds, frame_number, ...
                           background_block_size, subtract_background, ...
                           is_fluorescent)
    % Normalize the frame to its mean to be on the same scale as
    % the background image.
    frame = frame/mean(frame(:));

    background = 0;
    if subtract_background
        % subtract the background and invert if it is not a fluorescent
        % image.
        background = backgrounds(:,:,ceil(frame_number/background_block_size));
    end

    if is_fluorescent
        frame = frame - background;
    else
        frame = background - frame;
    end

    gauss_filter = fspecial('gaussian',[3 3], 0.5);
    frame = imfilter(frame, gauss_filter, 'replicate');
end


