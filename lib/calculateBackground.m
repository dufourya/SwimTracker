function backgrounds = calculateBackground(binFile, metadata, block_size)
    %  Extracts the background from a set of images.
    %
    %  The set of images is split into ceil(n_frames/block_size) blocks.
    %  Each of these image blocks is averaged to generate a background
    %  image for that block. If block_size does not split the set of images
    %  evenly, the final block is calculated by moving back block_size
    %  images from the end of the set.
    %
    % INPUTS:
    %   binFile          :  the location of the images.
    %   metadata         :  a metadata object.
    %   block_size       :  each background image will be generated from a set
    %                       of 'block_size' images.
    % OUTPUTS:
    %   backgrounds      : an array of background images
    
    width = metadata.getImageWidthPixel();
    height = metadata.getImageHeightPixel();
    
    pixel_type = metadata.getPixelType();
    bytes_per_pixel = metadata.getBytesPerPixel();
    n_frames = metadata.getNumberOfImages();
    n_before_frozen = metadata.getNImagesBeforeFreeze();
    has_frozen = 0;
    
    n_blocks = ceil(n_frames/block_size);
    backgrounds = [];
    if isempty(n_before_frozen)
        backgrounds = zeros(width,height,n_blocks);
        n_before_frozen = n_frames;
    else
        has_frozen = 1;
        backgrounds = zeros(width,height,n_blocks+1);
    end
    
    fid = fopen(binFile);
    %tic
    for i = 1:n_blocks
        if i==n_blocks
            % get last *non-frozen* block of images.
            fseek(fid, -width*height*(block_size+n_frames-n_before_frozen)* ...
                  bytes_per_pixel, 'eof');
        end
        backgrounds(:,:,i) = getBackground(fid, width, height, block_size, ...
                                           pixel_type);
    end
    
    if has_frozen
        backgrounds(:,:,i+1) = backgrounds(:,:,i);
    end
    %toc
    fclose(fid);
end

function bg = getBackground(fid,width,height,block_size,pixel_type)
    bg = fread(fid,width*height*block_size, pixel_type);
    bg = double(flipud(rot90(mean(reshape(bg,height,width,block_size),3))));
    bg = bg/mean(bg(:));
end
