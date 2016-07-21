function objs = im2obj_rp(im, objsize, thresh, fitstr)
%im2obj : finds objects in image stack
%
% im      : 3-d matrix of images
% objsize : size in pixels of objects to find. Can be a two-element array
%    (1) bpfiltsize, used to determine the low-pass filter size sent to 
%        bpass.m.  "0" indicates no filtering.
%    (2) nsize = size of the "neighborhood" around each particle, used to
%        make the structuring element for dilation (local max finding) and 
%        to set the array size for single-particle localization.
%    Can also be one element; fo4_rp.m will use the same value for both.
%    im2obj_rp.m only uses objsize for close-pair determination; if objsize
%        has two elements, use nsize (element 2)
% thresh  : number in [0 1] to threshold objects with
% fitstr : (Optional) string that selects the fitting option (e.g.
%    radial symmetry based, Gaussian fitting, centroid, ...
%    See fo4_rp.m for the full list of options. Default is radial symmetry
%   ('radial'). 
%
% Output:
% objs    : object matrix with following form:
%
%             objs = [x;
%                     y; 
%                     mass; (brightness)
%                     particleid; 
%                     frame; 
%                     trackid]
%           frame field is set by im2obj().  
%
% based on Andrew Demond's im2obj.m
% Modified by Raghuveer Parthasarathy, 16 April, 2007:
%   Calls fo4_rp -- RP's modified fo4.m -- no t-test, just std test
%   (deleted iput p : standard deviation value for maxima test in fo4_rp()
%           -- recommended: not > 1 )
% June 3, 2011: use distance.m for close pair calculation -- very small
% speedup
% Oct. 24, 2011: modify close pair criterion to be objsize, instead of
% 2*objsize
% July 10, 2012: avoid distance.m function call; keep the method!
% Last modified: July 10, 2012 (allow 2-component objsize)

% Fitting option
if ~exist('fitstr', 'var') || isempty(fitstr)
    disp('Default fitting option: radial-symmetry-based');
    fitstr = 'radial';
end
if strcmpi(fitstr, 'centroid')
    disp('Center of mass (centroid) fit -- AVOID unless necessary (saturated images)');
end

% Get nonlinear fitting options, to avoid repeated calls
% These are only used for non-linear Gaussian fitting, but it doesn't hurt
% to define them and pass them on to fo4_rp.m 
lsqoptions = optimset('lsqnonlin');

objs = [];
nf = size(im,3);
progtitle = sprintf('im2obj_{rp}: Finding objects...  '); 
if (nf > 1)
    progbar = waitbar(0, progtitle);  % will display progress
end
for j = 1:nf
    tmpobj = fo4_rp(im(:,:,j), objsize, thresh, fitstr, lsqoptions);
    tmpobj(5,:) = j;
    objs = [objs tmpobj];
    % show progress
    if mod(j,10)==0
        waitbar(j/nf, progbar, ...
            strcat(progtitle, sprintf('frame %d of %d', j, nf)));
    end
end
if nf>1
    close(progbar)
end


% -----------------------------------------------------------------------
% Check for likely multiple identifications of the same particle
%   --- find pairs within the same frame that are separated by <objsize
% If these are found, user should re-run the function with a larger
%   value of objsize.
% Won't display all close pair information -- user can re-examine
unqframes = unique(objs(5,:)); % get unique frame numbers
closepairlog = [];
for j=unqframes,
    objframe = objs(:,objs(5,:)==j);  
        % columns of the object matrix for this frame
    if size(objframe,2)>1
        % more than one particle found in this frame
        % consider all particle pairs, find any for which
        % separation < objsize
        allr = [objframe(1,:); objframe(2,:)];  % 2 x "N" matrix of x,y
        % Euclidean distance matrix.
        % uses method of "distance.m" by Roland Bunschoten (MATLAB file exchange)
        aa=sum(allr.*allr,1);
        d = sqrt(abs(aa( ones(size(aa,2),1), :)' + aa( ones(size(aa,2),1), :) - 2*allr'*allr));
        if length(objsize)<2
            isclose = (d < objsize);  % close pairs, and d for same pairs
        else
            isclose = (d < objsize(2));  % close pairs, and d for same pairs
        end
        closepair = (sum(isclose(:)) - size(objframe,2))/2.0;  
        templog = [j; closepair];
        closepairlog = [closepairlog templog];  % row 1: frame no; row 2: no. close pairs in this frame
    end
end
if ~isempty(closepairlog)
    if sum(closepairlog(2,:)>0)
        disp('Close pairs found -- recommend re-running with larger objsize.');
        clpairlogfr = closepairlog(:,closepairlog(2,:)>0);
        nclframes = size(clpairlogfr,2);
        medcl = median(clpairlogfr(2,:));
        fs = sprintf('   Number of frames with close pairs: %d out of %d total', nclframes, length(unqframes)); disp(fs);
        fs = sprintf('   In these, median number of close pairs: %d', medcl); disp(fs);
    end
end
