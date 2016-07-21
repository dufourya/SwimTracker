function metadata = sequenceAcquisition(mmc, file, runtime, config, fps, bit8, crop)
    % Record a movie, saved as raw binary output.
    %
    % INPUTS
    %   mmc     : the micro-manager object representing the microscope.
    %   file    : the output file name.
    %   runtime : the length of the movie to be recorded, in seconds.
    %   config  : the name of the objective configuration used to acquire the
    %               movie
    % OUTPUT
    %  A raw binary file in little-endian byte order.
    %  metadata = a Metadata object describing the acquisition.
    
    if nargin < 5 || isempty(fps) || fps < 0
        fps = 8;
    end
    if nargin < 6 || isempty(bit8) || bit8 ~=0
        bit8 = 1;
    end
    if nargin < 7 || isempty(crop) || crop ~= 1
        crop = 0;
    end
    
    %% microscope settings
    import mmcorej.mmc.*;
    f_name = getFunctionName();
    computer_name = getenv('COMPUTERNAME');
    camera = mmc.getCameraDevice();
    
    mmc.setConfig('System','Startup');
    mmc.setConfig('Channel',config);
    mmc.waitForSystem();
    checkTemperature(mmc);
    
    if mmc.getBytesPerPixel == 2
        pixelType = 'uint16';
        realbit = 16;
    else
        pixelType = 'uint8';
        realbit = 8;
    end
    
    emoscope_name = 'EMONETSCOPE';
    zenscope_name = 'T3610-EMONET';
    if strcmp(computer_name, emoscope_name)
        circ_buff_size_mb = 256;
    elseif strcmp(computer_name, zenscope_name)
        circ_buff_size_mb = 256;    
    else
        error(['[' f_name ']' ' Do not know memory requirements of computer ' ...
            '''' computer_name, '''.']);
    end
    
    if strcmp('HamamatsuHam_DCAM',camera)
        interval = str2double(mmc.getProperty(camera,'Exposure'));
    end
    
    %% open files to save data
    imgfile = strcat(file, FileInfo.binExt);
    metadata = Metadata(strcat(file, FileInfo.metaExt));
    metadata.create(mmc);
    
    fid = fopen(imgfile,'W');
    metadata.add('Experiment', 'Date', date);
    metadata.add('Experiment', 'Time', datestr(rem(now,1)));
    metadata.add('Experiment', 'DateVector', num2str(clock));
    metadata.add('Experiment', 'Computer', computer_name);
    metadata.add('Experiment', 'User', getenv('USERNAME'));
    
    %% start sequence acquisition
    mmc.setCircularBufferMemoryFootprint(circ_buff_size_mb);
    mmc.initializeCircularBuffer();
    
    mmc.prepareSequenceAcquisition(camera)
        
    if strcmp('HamamatsuHam_DCAM',char(camera))
        mmc.setProperty(camera,'TRIGGER SOURCE','INTERNAL');
    end
    
    w=mmc.getImageWidth();
    h=mmc.getImageHeight();
    
    mmc.waitForImageSynchro();
    mmc.startContinuousSequenceAcquisition(0);
    
    % HAVE to set interval AFTER starting sequence acquisition, since the camera
    % resets it when the method is called.
    if strcmp('Andor',char(camera))
        interval = str2double(char(mmc.getProperty('Andor','ActualInterval-ms')));
        realbit = 14;
    end
        
    nImages = ceil(runtime*1000/interval);    
    skip_frame = 1;
    
    if fps ~=0
        skip_frame = floor(1000/(fps * interval));
        if skip_frame < 1
            mmc.stopSequenceAcquisition();
            mmc.setShutterOpen(0);
            mmc.waitForSystem();
            fclose(fid);
            error('The target aquisition frame rate cannot be achieved with this configuration! Use fps = 0 if you want to record at the maximum possible frame rate.');
        end
    end
    
    if bit8 == 1
        fprintf('Converting images to 8 bits\n');
    else
        fprintf('Recording images in 16 bits\n');
    end
    if crop == 1
        fprintf('Frame is cropped by 1/2\n');
    else
        fprintf('Recording full frame\n');
    end    
    fprintf('Skipping %d frame(s)\n', skip_frame-1);
	fprintf('Image interval = %.2f ms\n', skip_frame*interval);
    fprintf('Recording at %.2f frames/sec\n', 1000/(skip_frame*interval));
    
    i=0;
    tic1=tic;
    nImages_actual = 0;
    
    while i < nImages
        if mmc.getRemainingImageCount()>0
            img = mmc.popNextImage();
            if mod(i, skip_frame) == 0
                img = typecast(img, pixelType);
                if crop == 1
                    img = reshape(img,w,h);
                    img = img((w/2-w/4+1):(w/2+w/4), (h/2-h/4+1):(h/2+h/4));
                    img = img(:);
                end
                if bit8 == 1 && realbit ~= 8
                    img = uint8(double(img) * (2^8 - 1) / (2^realbit - 1));
                    fwrite(fid, img, 'uint8');
                else
                    fwrite(fid, img, pixelType);
                end
                nImages_actual = nImages_actual + 1;
            end
            i = i+1;
            if mod(i,ceil(60000/interval)) == 0
                toc(tic1);
            end
        end
    end
    toc(tic1);
    
    mmc.stopSequenceAcquisition();
    mmc.waitForImageSynchro();
    mmc.setShutterOpen(0);
    
    fprintf('Images recorded: %d\n', nImages_actual);
    
    %% write more metadata
    metadata.add('SequenceAcquisition', 'ActualIntervalBurst-ms', ...
                 sprintf('%f', skip_frame * interval));
    metadata.add('SequenceAcquisition', 'NumberImages', ...
                 sprintf('%d', nImages_actual));
    
    if bit8 == 1
        metadata.add('SequenceAcquisition', 'FilePixelType', '*uint8');
    else
        metadata.add('SequenceAcquisition', 'FilePixelType', ...
                     strcat('*', pixelType));
    end
    
    if crop == 1
        metadata.set('Camera', 'ImageWidth', ...
                     sprintf('%f', mmc.getImageWidth()/2));
        metadata.set('Camera', 'ImageHeight', ...
                     sprintf('%f', mmc.getImageHeight()/2));
    end
    metadata.write();
end
