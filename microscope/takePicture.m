function img = takePicture(mmc, config, crop, save, img_path)
    if nargin < 5
        img_path = [];
    end
    
    if save && isempty(img_path)
        error('Must specify save path for image.')
    end
    if ~save
        img_path = 'dummy';
    end
        
    mmc.setConfig('Channel',config);
    mmc.waitForSystem();
    mmc.snapImage();
    imgtmp = mmc.getImage();
    
    pixel_type = [];
    if mmc.getBytesPerPixel == 2
        pixel_type = 'uint16';
    else
        pixel_type = 'uint8';
    end
    img = rotateFrame(typecast(imgtmp, pixel_type), mmc);

    w=mmc.getImageWidth();
    h=mmc.getImageHeight();

    [path, img_name, ext] = fileparts(img_path);
    metadata = Metadata(fullfile(path, strcat(img_name, FileInfo.metaExt)));
    metadata.create(mmc);
    
    if crop
        img = img((w/2-w/4+1):(w/2+w/4), (h/2-h/4+1):(h/2+h/4));
        metadata.set('Camera', 'ImageWidth', sprintf('%f', w/2));
        metadata.set('Camera', 'ImageHeight', sprintf('%f',h/2));
    end
    
    if save
        imwrite(img, img_path, 'TIFF');
        metadata.write();
    end
end
