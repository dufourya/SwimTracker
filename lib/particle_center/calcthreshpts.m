% calcthreshpts.m
% 
% Function to determine the points that "pass" the threshold criterion in
% an image.  
% Usual application: particle tracking.  Typically first filter the image (e.g. w/ 
% bpass.m) before inputting to this function.  This function
% thresholds (various options) and looks for local maxima around which to 
% localize particles by analyzing the original (non-filtered) image.
%
% Three thresholding options, that determine the "thresh" variable.
%
% Inputs:
%    A : 2D image.
%    threshopt: determined by fo4_rp.m
%      (1) keeps all pixels with intensity above the thresh*100 percentile
%      (2) keeps pixels with intensity > thresh*std.dev. above
%          the median intensity
%      (3) keeps brightest "thresh" number of particles by finding the 
%          "thresh" brightest regions and allowing only one maximum in
%          each.
%          Note that the number of particles can be less than thresh, if
%          fewer than thresh local maxima exist.
%    thresh : intensity threshold
%          For option 2, sign should be positive (flipped by fo4_rp.m) --
%          the function makes sure of this.
%    nsize : necessary only for threshopt==3. object size used for dilating,
%        determining neighborhood size.  Roughly, the particle diameter.
%        Referred to as nsize in > June 27, 2012 fo4_rp.m,
%        since filtering object size is separated from neighborhood object
%        size; formerly called objsize.
%    try1pernhood : (optional; default false).  If true, dilate local maxima to
%        try to have only one local max per neighborhood.  Dilate by radius 
%        nsize/4 and take centroid of dilated max region.  
%        (This is always done for threshold option 3, regardless of
%        this input variable.)
%          
% Outputs
%    [y, x] : positions of local post-threshold maxima
%
% Raghuveer Parthasarathy
% extracted from fo4_rp.m, May 2, 2012
% June 28, 2012: write "try1pernhood" option
% last modified July 3, 2012


function [y, x] = calcthreshpts(A, threshopt, thresh, nsize, try1pernhood)

if ~exist('try1pernhood', 'var') || isempty(try1pernhood)
    try1pernhood = false;
end

if threshopt==2
    thresh = abs(thresh);  % make sure it's been flipped
end

showplots=false;  % for debugging -- plot things.

ste = strel('disk', floor(nsize/2),0);  % for dilation
dimg = imdilate(A, ste);
if showplots
    figure(3)
    imshow(dimg,[]); title('3 dilated')
end
% finding regional maxima (imregionalmax(A)) is much faster than A==dilation.
% However, it is inaccurate, as it does not get rid of local maxima that
% are close (within nsize/2) to brighter true particle maxima.
% imgmax = imregionalmax(A);

% Local maxima in the filtered image
BW = (A == dimg);
if showplots
    figure(4)
    imshow(BW); title('4 local maxima')
end

qradius = max(floor(nsize/4), 1);  % For 1 per region; force qradius >= 1
switch threshopt
    case 1
        % intensity thresholding (percentile)
        [hs, bins] = hist(A(:),1000);
        ch = cumsum(hs);
        ch = ch/max(ch);
        noiseind = find(ch > thresh); %
        noiseind = noiseind(1); % The index value below which "thresh" fraction
        % of the pixels lie.  (Originally noisind(2), but this can lead to errors
        % if there is only one bin above threshold.)
        % find local maxima that are also above-threshold
        if try1pernhood
            stats = dilatepeaks(A, BW, qradius); % move to a nested function
            stats([stats.MaxIntensity] <= bins(noiseind))=[]; % remove dim spots
            [x y] = getxycentroid(stats, length(stats)); % move to a nested function
        else
            % look at each local maximum
            [y,x] = find(BW & (A > bins(noiseind))); % "the magic happens here!" (Andrew DeMond)
        end
        if showplots
            figure(5)
            imshow(zeros(size(BW))); hold on; 
            plot(x,y,'wo', 'markersize', 2, 'markerfacecolor', 'w'); 
            title('5 local maxima > threshold')
        end
    case 2
        % intensity thresholding (std. dev. above background)
        medint = median(A(:));
        stdint = std(A(:));
        % find local maxima that are also above-threshold
        if try1pernhood
            stats = dilatepeaks(A, BW, qradius); % move to a nested function
            stats([stats.MaxIntensity] <= (medint + thresh*stdint))=[]; % remove dim spots
            [x y] = getxycentroid(stats, length(stats)); % move to a nested function
        else
            isbright = (A > (medint + thresh*stdint));
            % look at each local maximum
            [y,x] = find(BW & isbright);
        end
    case 3
        % thresh >= 1 : keep brightest "thresh#" of particles
        % Can't just look for the thresh# of brightest spots -- too sensitive
        % to bright noise.  Need to consider "regions"
        thresh = round(thresh);  % make sure it's an integer
        stats = dilatepeaks(A, BW, qradius); % move to a nested function
        % find the thresh# brightest points (region properties) 
        [filtmaxint, fmix] = sort([stats(:).MeanIntensity],'descend');
        % find thresh number of particles, or number of local maxima, whichever is smaller
        numtofind = min([thresh length(fmix)]);  
        x = zeros(numtofind,1); y = zeros(numtofind,1);
        for k=1:numtofind
            xy =  stats(fmix(k)).Centroid;
            x(k) = xy(1);
            y(k) = xy(2);
        end
        
end

        function BWstats = dilatepeaks(A, BW, qradius)
            % Dilate, so that local maxima in the same nsize region are one
            % "object."  Dilate by a half-sized structuring element.  (Closing
            % by ste doesn't work -- disjoint points remain points
            BWdil = imdilate(BW, strel('disk', qradius));
            BWstats = regionprops(BWdil, A, 'Centroid', 'MeanIntensity', 'MaxIntensity');
        end

        function [xc yc] = getxycentroid(stats, N)
            xy = [stats.Centroid];
            xc = xy(1:2:(2*N-1));
            yc = xy(2:2:(2*N));
        end
end
