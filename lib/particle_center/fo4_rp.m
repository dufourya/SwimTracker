% fo4_rp.m
% fo4_rp : finds objects in images, first filtering and determining local
% maxima (calls calcthreshpts.m),
% and then refining positions with one of several possible methods
%
% Inputs:
% img = image to locate objects within
% objsize : size in pixels of objects to find. Can be a two-element array
%    (1) bpfiltsize, used to determine the low-pass filter size sent to 
%        bpass.m.  "0" indicates no filtering.
%    (2) nsize = size of the "neighborhood" around each particle, used to
%        make the structuring element for dilation (local max finding) and 
%        to set the array size for single-particle localization.
%    If one element, use the same value for both
% thresh : intensity threshold 
%    *Various options,* inferred from the form of the input
%    (1) if a number in [0, 1), sets the local max intensity threshold, 
%        keeping all points with intensity above the thresh*100 percentile
%    (2) if a number <0, keeps pixels with intensity > thresh*std.dev. above
%          the median intensity
%    (3) if a number >= 1, keeps brightest "thresh" number of particles
%          by finding the "thresh" brightest regions and allowing only one
%          maximum in each
% fitstr : (Optional) string that selects the fitting option 
%          -- [default] 'radial'.  Radial-symmetry based fit.  Fast, accurate -- see
%             notes July-August 2011, and RP's Nature Methods paper, June
%             2012.
%             If the number of points in the image is <10, call
%             radialcenter.m .  Otherwise, call radialcenter_stk.m (stack
%             input, avoiding redundant grid calculations for speed).
%          --  'gaussmle' , radially symmetric 2D Gaussian, fit by maximum
%             likelihood estimation.  Slow, most accurate.
%          --  'nonlineargauss' , radially symmetric 2D Gaussian
%             nonlinear fit.  Slow, accurate.
%             Fit to intensity = offset  A*exp(-(x^2+y^2)/2/sigma^2);
%            The "intensity" (row 3 of the object matrix) is pi*A*sigma_x*sigma_y
%          -- 'lineargauss' , 2D Gaussian, linear fit (i.e. parabolic 
%             fit to log(intensity).  Fast, moderate accuracy
%          -- 'centroid'.  Centroid fit.  Least accurate (heavily biased); 
%             may be necessary for particles with saturated image
%             intensities.
%          -- 'weightedlineargauss' , Linearized Gaussian fit, weighted 
%             for noise (Stephen M. Anthony, Steve Granick -- see Langmuir paper)
%             DON'T USE THIS
% lsqoptions : [optional] options structure for nonlinear least-squares
%              fitting, from previously running 
%              "lsqoptions = optimset('lsqnonlin');"
%              Inputting this speeds up the function by avoiding redundant 
%              calls to optimset.
%              This variable is only used for non-linear Gaussian fitting.
%
% Output
%   objs : array of object properties; one column per object.  Rows:
%     [x;
%      y;
%      mass;  (i.e. brightness)
%      particleid;
%      frame;
%      trackid;
%      sigma]
% sigma is the 2D Gaussian width (std)('nonlineargauss'), the avg. of the x and
%    y Gauss fit ('lineargauss') or is zero ('centroid')
% frame and trackid are set to zero for all objects by fo4() and must be
% dealt with by the functions calling fo4().
%
% Raghuveer Parthasarathy, April 2007
% Modifications -- 
%   March 24, 2011 (option for nonlinear 2D Gaussian fit)
%   March 31, 2011 (returns the Gaussian width, sigma)
%   June 8, 2011: Allow threshold to set a max no. of particles
%   August 8, 2011: Allowing thresholding to find centers of only 
%      the 'n' brightest objects
%   August 31, 2011: Allow radial symmetry based fitting
%   October 24, 2011: std. dev. based thresholding
%   Feb. 11, 2012: add max. likelihood Gaussian fit
%   May 2, 2012 : move thresholding to a function: calcthreshimg.m
%   June 27, 2012: allow two-element object size array 
%   July 3, 2012: Calls radialcenter_stk.m if there are >10 pts for radial
%     symmetry fitting, for greater speed.
% Last modified June 27, 2012

function [objs] = fo4_rp(img, objsize, thresh, fitstr, lsqoptions)

% Fitting option; default is nonlinear Gaussian
if ~exist('fitstr', 'var') || isempty(fitstr)
    fitstr = 'radial';
end
if ~exist('lsqoptions', 'var') || isempty(lsqoptions)
    % This variable is only used for non-linear Gaussian fitting.
    if strcmpi(fitstr, 'nonlineargauss')
        lsqoptions = optimset('lsqnonlin');
    else
        lsqoptions = [];
    end
end

% Size parameters
if length(objsize)==1
    bpfiltsize = objsize;
    nsize = objsize;
else
    bpfiltsize = objsize(1);
    nsize = objsize(2);
end

showplots=false;  % for debugging -- plot things.

% Determine thresholding option -- see header comments 
if thresh >= 1.0
    threshopt = 3;
elseif thresh >= 0.0
    threshopt = 1;
else
    threshopt = 2;
    thresh = -thresh;  % note that negative thresh is the indicator of this option;
                       % flip so it's positive.
end

% make img processable
img = double(img);

% Bandpass filter (if nsize > 0) -- use Grier et al bpass.m
if bpfiltsize>0
    noisesize = 1;  % size of noise, px
    filtimg = bpass(img,noisesize,bpfiltsize);
else
    filtimg = img;
end

if showplots
    figure(1)
    imshow(img,[]); title('1 original image')
    figure(2)
    imshow(filtimg,[]); title('2 bandpass filtered')
end

% Three options for thresholding -- see above
% For the chosen option, determine the points that "pass" the threshold.
% Move to a separate function, so it can be called by the GUI
[y, x] = calcthreshpts(filtimg, threshopt, thresh, nsize);

% Get rid of maxima too close to the edge
lenx = 2*floor(nsize/2) + 1;  % 'floor' isn't really necessary, but this 
    % is the size of "nhood = getnhood(ste);" for a disk structuring
    % element of that size
leny = lenx;  % in principle, could make different
edgeind = ((x < lenx/2) | (x > (size(img,2) - lenx/2))) | ...
    ((y < leny/2) | (y > (size(img,1) - leny/2)));
x(edgeind) = [];
y(edgeind) = [];
   
% Compute "masses"
savemass = zeros(1, length(x));
rect = zeros(length(x), 4);
% Compute the first neighborhood of a maximum, to know the image size
% Can skip if there are no objects to find
if ~isempty(x)
    rect(1,:) = [round(x(1) - lenx/2) round(y(1) - leny/2) (lenx-1) (leny-1)];
    cropimg1 = imcrop(img, rect(1,:));
    % all the other neighborhoods
    cropimg = repmat(cropimg1, [1 1 length(x)]);
    for k = 2:length(x)
        rect(k,:) = [round(x(k) - lenx/2) round(y(k) - leny/2) (lenx-1) (leny-1)];
        cropimg(:,:,k) = imcrop(img, rect(k,:));
    end
end

% Calculate "mass" (intensity)
nhood = getnhood(strel('disk', floor(nsize/2),0));  % somewhat silly
for k = 1:length(x)
    tempreg = cropimg(:,:,k);
    if(size(cropimg1) == size(nhood))
        cropreg = tempreg(nhood);
    else
        cropreg = tempreg;
    end
    savemass(k) = sum(cropreg(:));
end

% Do refinement (find center) around each local max
% If there are many points, use radialcenter_stk.m rather than radialcenter.m
% for radial symmetry method -- even faster!
xcent = zeros(1,length(x));
ycent = zeros(1,length(x));
sigma = zeros(1,length(x));
lsumx = 1:size(cropimg,2);
lsumy = 1:size(cropimg,1);
Lx = lsumx(end);
Ly = lsumy(end);
% for j = 1:length(x)
    switch lower(fitstr)
        case {'radial'}
            % Radial-symmetry based fit -- fast, accurate
            % If <10 points, use radialcenter_stk.m for extra speed (avoid
            % redundant grid calculations); else radialcenter.m
            if length(x) < 10
                for j = 1:length(x)
                    [xcent(j), ycent(j), sigma(j)] = radialcenter(cropimg(:,:,j));
                end
            else
                [xcent ycent sigma] = radialcenter_stk(cropimg) ;
            end
            % Is the center within reasonable bounds?
            % If not, replace with centroid
            % frequency of bad cases ~ 1/100,000 !  (Extremely rare)
            % See notes Oct. 26, 2011
            % This conditional statement can slow things; Delete?
            badcase = find(abs(xcent - Lx/2)>1.5*Lx | abs(ycent - Ly/2)>1.5*Ly);
            for j = badcase
                ci = cropimg(:,:,j);
                xcent(j) = sum(sum(ci) .* lsumx) / sum(sum(ci));
                ycent(j) = sum(sum(ci,2) .* lsumy') / sum(sum(ci));
            end
        case {'gaussmle'}
            % Gaussian fit via maximum likelihood estimmation -- most accurate
            for j = 1:length(x)
                [A, xcent(j), ycent(j), sigma(j), offset] = gaussfit2DMLE(cropimg(:,:,j));
            end
        case {'nonlineargauss'}
            % Gaussian fit via nonlinear least squares
            for j = 1:length(x)
                [A, xcent(j), ycent(j), sigma(j), offset] = gaussfit2Dnonlin(cropimg(:,:,j), [], [], [], [], lsqoptions);
                % savemass(j) = A*2*pi*sigma(j)*sigma(j);  % Area under a 2D gaussian
                %        Don't use this, since gives large values for weak Gaussians
            end
        case {'lineargauss'}
            % Linear gaussian fit
            % Linear regression -- fit of cropimg intensity to a 2D Gaussian,
            % via polynomial fit of log(intensity),
            % using gaussfit2D.m with a 0.2 threshold
            for j = 1:length(x)
                [A, x0, sigma_x, y0, sigma_y] = gaussfit2D(lsumx, lsumy, cropimg(:,:,j), 0.2);
                if imag(x0)>0.0
                    xcent(j) = 0.0;  % return zero -- nonsense
                else
                    xcent(j) = x0;
                end
                if imag(y0)>0.0
                    ycent(j) = 0.0;  % return zero -- nonsense
                else
                    ycent(j) = y0;
                end
                sigma(j) = 0.5*(sigma_x+sigma_y);  % mean Gaussian width
            end
        case {'weightedlineargauss'}
            % Linearized Gaussian fit, weighted for noise (Stephen M.
            % Anthony, Steve Granick -- see Langmuir paper)
            % Need "noiselevel"  std dev of background noise
            % Function gauss2dcirc.m written by Stephen M.
            % Anthony.
            noiselevel = 220;
            disp('hardwiring noise level!!!  -- re-write this')
            for j = 1:length(x)
                [xcent(j),ycent(j),A,sigma(j)] = ...
                    gauss2dcirc(cropimg(:,:,j),repmat(lsumx,size(cropimg,1),1),...
                    repmat(lsumy',1,size(cropimg,2)),noiselevel);
            end
        case {'centroid'}
            % centroid (center of mass) fit
            % don't subtract background
            % consider all points at once (not looping)
            sumcropimg = sum(sum(cropimg,1),2); % length = length(x)
            xcent = sum(sum(cropimg,1).*repmat(lsumx, [1,1,length(x)]),2) ./ sumcropimg;
            xcent = reshape(xcent, [1 length(x)]);
            ycent = sum(sum(cropimg,2).*repmat(lsumy', [1,1,length(x)]),1) ./ sumcropimg;
            ycent = reshape(ycent, [1 length(x)]);
        otherwise
            disp('Unknown method! [fo4_rp.m]');
    end
%end
% center position relative to image boundary
xn = xcent + rect(:,1)' - 1; % -1 is to correct for matlab indexing
yn = ycent + rect(:,2)' - 1;

if showplots
    figure(6)
    imshow(zeros(size(img))); title('6 particle centers')
    for j=1:length(xn)
        rectangle('Position', [xn(j)-1 yn(j)-1 2 2], 'Curvature', [1 1], ...
            'Linewidth', 2.0, 'EdgeColor', [1.0 1.0 0.3]);
    end
end

objs = zeros(6, length(xn));
objs(1,:) = xn;
objs(2,:) = yn;
objs(3,:) = savemass;
objs(4,:) = 1:length(x);
objs(7,:) = sigma;
