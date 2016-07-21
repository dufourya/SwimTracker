function convertMetadataFiles(camera, image_mag)
    % converts metafiles in current directory from old format to new format.
    % camera     : 'ixon' or 'flash'
    % image_mag  : the objective magnification used to take the movie. 
    %              e.g., '10'
    
    pix_info = PixelInfo(camera);
    files = dir(strcat('*', FileInfo.metaExt));
    for file = files'
        file_info = FileInfo(file.name);
        metadata = Metadata(file_info.metaFile);
        fprintf('converting %s...\n', file.name);
        pixel_size = metadata.getPixelSize();
        width = metadata.getImageWidthPixel();
        height = metadata.getImageHeightPixel();
        
        if isempty(pixel_size)
            metadata.appendPixelSize(pix_info.getUmPerPixel(image_mag));
            added('pixel size');
        else
            already('pixel size');
        end
        
        if isempty(width)
            metadata.appendImageWidthPixel(pix_info.getImageWidthPixel());
            added('image width');
        else
            already('image width');
        end
        
        if isempty(height)
            metadata.appendImageHeightPixel(pix_info.getImageHeightPixel());
            added('image height');
        else
            already('image height');
        end
        
        fprintf('done\n');
    end
end

function already(val)
    fprintf('    ''%s'' already set\n', val);
end
function added(val)
    fprintf('    added entry for ''%s''\n', val);
end
