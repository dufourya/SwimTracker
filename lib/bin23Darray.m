function frame_array = bin23Darray(fid, metadata)

    n_frames = metadata.getNumberOfImages();
    width = metadata.getImageWidthPixel();
    height = metadata.getImageHeightPixel();
    pixel_type = metadata.getPixelType();

    frewind(fid);
    frame_array = zeros(width, height, n_frames);
    for i=1:n_frames
        frame_array(:,:,i) = reshape(fread(fid,width*height, pixel_type), ...
                                     [width height]);
        
    end
    fclose(fid);
end