function rad_detection = radDetect(file_info, metadata, params)

diagnostic_plot = 0;

frac_mem_usage = params.runTrackSingle.mem_usage;
% thresh = params.radcenter.detection_threshold;
cell_area = params.radcenter.cell_area;
f_name = getFunctionName();

resolution = 2*metadata.getResolution;

width = metadata.getImageWidthPixel();
height = metadata.getImageHeightPixel();

mask = true(height,width);
if exist(file_info.maskFile,'file')
    fprintf('Reading mask image...');
    mask = flipud(imread(file_info.maskFile));
    mask = mask == max(mask(:));
    fprintf('done.\n');
end

pixel_type = metadata.getPixelType();
n_frames = metadata.getNumberOfImages();
bg_sample_time = params.runTrackSingle.background_sample_time;
imageInterval = metadata.getImageIntervalMs()/1000;
bg_sample_n = min(ceil(bg_sample_time/imageInterval),floor((n_frames-1)/2));

pixel_size = metadata.getPixelSize();
bytes_per_pixel = metadata.getBytesPerPixel();
image_size = width*height*bytes_per_pixel/1e6;

isFluo = metadata.isFluorescent();

removeHalo = 0;
% if contains(metadata.getVal('TINosePiece','Label'),{'40x','100x'})
if strcmp(params.objects,'microspheres')
   removeHalo = 1;
end

% fo4_rp options
object_radius = (resolution + sqrt(cell_area/pi)) / pixel_size;
scalewavelet = max(2,1.5*object_radius);
ste = strel('disk',ceil(scalewavelet),0);  % for dilation

% if ~isFluo
%     object_radius = 1.5*object_radius;
% end

detected_objs = cell(n_frames,1);
fprintf('[%s] Detecting cells in frames...\n', f_name);

% freemem = determineMemory();
% conversion to double (64-bit) increases memory requirement by a
% factor of 4.
% mem_to_use_mb = freemem * frac_mem_usage;
% n_blocks = ceil(2*image_size*n_frames / mem_to_use_mb);
n_blocks = 1;
block_size = ceil(n_frames/n_blocks);
% block_size = 200;
tic
fid = fopen(file_info.binFile,'r');
for i = 1:n_blocks
    
    fprintf('Reading block %d of %d...', i, n_blocks);
    frames = fread(fid,width*height*block_size,pixel_type);
    frames = reshape(frames, width, height, []);
    fprintf('done.\n');
    
    fprintf('Calculating threshold...');
    bkgsample1 = cell(10,1);
    ii =1;
    
   a = 1:ceil(n_frames/12):n_frames;

    for j = a(2:11)
        ind1 = j-bg_sample_n:j+bg_sample_n;
        ind1(ind1<1) = [];
        ind1(ind1>n_frames) = [];
        ind = setdiff(ind1,j);
        frame = double(rot90(frames(:,:,j)));
        if bg_sample_n>0
            background = rot90(mean(double(frames(:,:,ind)),3));
            frame = frame - background;
        end
        if ~isFluo
            filtimg = cwtft2(frame,'wavelet','mexh','scales',scalewavelet);
        else
            filtimg = cwtft2(-frame,'wavelet','mexh','scales',scalewavelet);
        end
        immax = imregionalmax(filtimg.cfs);
        immax(~mask) = 0;
        %         filtimg = filtimg.cfs(randperm(numel(filtimg.cfs),10000));
        bkgsample1{ii} = filtimg.cfs(immax);
        ii = ii +1;
        %         bkgsample2 = [bkgsample2; filtimg(:)];
    end
    
    %%
    nbg = min(cellfun(@numel,bkgsample1));
    bkgsample2 = zeros(nbg,10);
    for j = 1:numel(bkgsample1)
        idx = randperm(numel(bkgsample1{j}));
        bkgsample2(:,j) = bkgsample1{j}(idx(1:nbg));
    end
    bkgsample1 = imquantnorm(bkgsample2);
    bkgsample1 = bkgsample1(:);
    %%
    options = statset('MaxIter',1000000, 'MaxFunEvals',1000000);
    ev = fitdist(-bkgsample1(:),'ExtremeValue');
    mix_pdf = @(x,mu1,sig1,mu2,sig2,p) p*evpdf(x,mu1,sig1) + (1-p)*evpdf(x,mu2,sig2) + 10^-37;
    phat = mle(-bkgsample1(:),'pdf',mix_pdf,'start',[ev.mu/2 ev.sigma/2, 2*ev.mu 2*ev.sigma, 0.99],...
        'lower',[-Inf 0 -Inf 0 0], 'upper',[Inf Inf Inf Inf 1],'options',options);
    %%
    
    mix_cdf_params = @(x) phat(5).*evcdf(x,phat(1),phat(2)) + (1-phat(5)).*evcdf(x,phat(3),phat(4));
    sig_cdf_params = @(x) (1-phat(5)).*evcdf(x,phat(3),phat(4));
    bg_cdf_params  = @(x) phat(5).*evcdf(x,phat(1),phat(2));
    
    trs = [fzero(@(x) sig_cdf_params(x)/mix_cdf_params(x)-(1-params.runTrackSingle.fdr_detection),-100),...
        fzero(@(x) bg_cdf_params(x)-params.runTrackSingle.pval_detection,-100)];

%        fzero(@(x) mix_cdf_params(x)-params.runTrackSingle.max_object_proportion,-100)];
%     -trs
    thresh = min(sum(bkgsample1(:)>=-min(trs))/numel(bkgsample1(:)),params.runTrackSingle.max_object_proportion);
    params.runTrackSingle.calculated_object_proportion = thresh;
    fprintf('%f...done.\n', thresh);
    
    %%
    if diagnostic_plot
        mix_pdf_params = @(x) phat(5)*evpdf(x,phat(1),phat(2)) + (1-phat(5))*evpdf(x,phat(3),phat(4));
        sig_pdf_params = @(x) (1-phat(5))*evpdf(x,phat(3),phat(4));
        bg_pdf_params  = @(x) phat(5)*evpdf(x,phat(1),phat(2));
        h = figure('color','w','Name','Diagnostics');
        subplot(1,3,1),
        histogram(-bkgsample1(:),'normalization','pdf');
        hold on, fplot(@(x) mix_pdf_params(x));
        fplot(@(x) sig_pdf_params(x));
        fplot(@(x) bg_pdf_params(x));
        plot([-min(trs) -min(trs)],[10^-10 1]);
        set(gca,'YScale','log');
%         ylim([10^-7 10^-2]);
%         xlim([-50000 10000]);
        xlabel('Background signal fit');
        axis square
        hold off;
    end
    
    %%
    
    fprintf('Processing images...');
    for j = 1:size(frames,3)
        frame_number = (i-1)*block_size+j;
        %         fprintf('processing frame %d\n', frame_number);
        
        % Normalize the frame to its mean to be on the same scale as
        % the background image.
        frame = rot90(frames(:,:,j));
        
        if diagnostic_plot
            figure(h),
            subplot(1,3,2),
            imshow(imadjust(frame,stretchlim(frame,0.001)),[]);
            xlabel('Frame');
        end
        
        frame = double(frame);
        ind1 = j-bg_sample_n:j+bg_sample_n;
        if ind1(1) < 1
            ind1 = ind1-ind1(1)+1;
        end
        if ind1(end)>size(frames,3)
            ind1 = ind1-ind1(end)+size(frames,3);
        end
        
        ind = setdiff(ind1,j);
        
        if bg_sample_n>0
            background = rot90(mean(double(frames(:,:,ind)),3));
            framecorr = frame - background;
        else
            framecorr = frame;
        end
        
        if ~isFluo
            framecorr = cwtft2(framecorr,'wavelet','mexh','scales',scalewavelet);
            frame = -frame;
        else
            framecorr = cwtft2(-framecorr,'wavelet','mexh','scales',scalewavelet);
        end
        
        framecorr = framecorr.cfs;
        %         framecorr = framecorr/sigmabackground;
        tmpobj = fo4_rp_ysd(framecorr, frame, mask, thresh, ste, removeHalo);
        %         size(tmpobj)
        
        if diagnostic_plot
            figure(h),
            subplot(1,3,3),
            imshow(framecorr,[]);
            hold on;
            plot(tmpobj(1,:), tmpobj(2,:),'g.')
            hold off;
            xlabel('Background - Frame');
            drawnow;
        end
        
        if ~isempty(tmpobj)
            % set the frame (not done automatically by fo4_rp).
            tmpobj(5,:) = frame_number;
            %             tmpobj(3,:) = tmpobj(3,:);
            detected_objs{frame_number} = tmpobj;
        end
    end
    fprintf('done.\n');
    fprintf('[%s] %d of %d\n', f_name, min(i*block_size, n_frames), ...
        n_frames);
end
detected_objs = cell2mat(detected_objs');

% %%
% histogram((detected_objs(3,:))./detected_objs(7,:))
% %%
% options = statset('MaxIter',1000);
% gm_objs = fitgmdist([log(detected_objs(3,:));detected_objs(7,:)]',2,'Options',options,...
%     'Replicates',5);
% bic = gm_objs.BIC;
% diff_bic = Inf;
% k=1;
% while diff_bic > 0
%     best_gm = gm_objs;
%     k= k+1;
%     gm_objs = fitgmdist([log(detected_objs(3,:));detected_objs(7,:)]',k,'Options',options,'CovarianceType','full',...
%         'RegularizationValue',10^-10,'Replicates',5);
%     diff_bic = bic - gm_objs.BIC;
%     bic = gm_objs.BIC;
% end
% best_gm
% %%
% clusters = best_gm.cluster([log(detected_objs(3,:));detected_objs(7,:)]');
% scatter(log(detected_objs(3,:)),detected_objs(7,:),[],clusters')
%%
if ~isempty(detected_objs)
    % convert the output to the form required by utrack.
    rad_detection = rad2utrack(detected_objs, metadata);
else
    rad_detection = [];
end

save(file_info.detectionFile, 'detected_objs', ...
    'rad_detection' , 'params');

fprintf('[%s] done.\n', f_name);
toc
end
