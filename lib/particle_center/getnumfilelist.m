% getnumfilelist.m
%
% Function to get the "base" file name and starting and ending frames for a
% series of images with filenames 'NAMExxxx', where NAME is the same for
% all images and xxxx is a sequential numbering.  
% Ideally, xxxx will be the same length for all images.  If file name
% length differs, prompt user for fbase and formatstr, 
% with defaults being shared characters and '%d' (no leading zeros)
% If the second file name is the same as the first, return empty arrays 
% for frmin and frmax (see LSseries.m for an example of dealing with this). 
% Can also be used for multi-page TIFFs, in which case the second file name
% is output as empty.  Function detects that the "first" filname chosen is
% a multipage TIFF
%
% originally extracted from TIFFseries.m
%
% Inputs:
%    FileName1, FileName2 : (Optional) FileName1 and FileName2 (empty if #1
%    is a multi-page TIFF), to bypass dialog box
% Outputs
%    fbase : base file name string (NAME) -- file name for a multipage TIFF
%    frmin : min frame number (i.e. min xxxx) -- "1" for a multipage TIFF
%    frmax : max frame number (i.e. min xxxx) -- "#frames" for a multipage TIFF
%    formatstr : formatting string to create file names
%    FileName1 : string of first file name -- not really necessary
%    FileName2 : string of last file name -- not really necessary
%    PathName1 : path of the (first) file
%    ext : file extension (e.g. 'tif')
%    ismultipage : boolean true if considering a multipage TIFF (redundant
%            with FileName2 being empty)
%
% EXAMPLE
% To use the output variables to create file names , e.g.:
%        k = 34;
%        framestr = sprintf(formatstr, k);
%        FileName = strcat(fbase, framestr, '.tif');
%        A  = imread(FileName, 'tif');  % image
% (This example doesn't use the ext variable, and assumes a TIFF format)
%
% EXAMPLE
% To out the output variables to load an image from a mulitpage TIFF:
%       k = 34;
%       A  = imread(strcat(fbase, sprintf(formatstr,k), ext));
% Raghuveer Parthasarathy
% Sept. 15, 2010
% July 10, 2011: Get rid of initial dialog, and go immediately to file 
%    selection box.Allow multi-page TIFF input
% last modified: May 17, 2012

function [fbase, frmin, frmax, formatstr, FileName1, FileName2, PathName1 ext ismultipage] = ...
    getnumfilelist(FileName1, FileName2)

firstdir = pwd;
if nargin > 0
    % filname(s) are specified, avoid dialog box.
    [path, FileName1base, ext] = fileparts(FileName1);  % to get the extension (includes ".")
    PathName1 = firstdir;
else
    % Dialog boxes
    % dialog box for filenames
    [FileName1,PathName1] = uigetfile('*.*', 'First .tif image, or multi-page TIFF...');
%    cd(PathName1)
    Nf1 = numel(imfinfo(strcat(PathName1,FileName1)));
%    cd(firstdir)
    if  Nf1 > 1
        % Multipage TIFF
        FileName2 = [];
    else
        cd(PathName1)
        [FileName2] = uigetfile('*.*', 'Last .tif image of sequence...');
        cd(firstdir)
    end
end

% The following block allows the user-input filenames to be used either for
% a series of TIFF files or a single multi-page TIFF; for the latter, leave
% Filename2 empty
if isempty(FileName2)
    ismultipage = true;
else
    ismultipage = false;
end

if ismultipage
    frmin = 1;
    frmax = Nf1;
    [pathstr, fbase, ext] = fileparts(FileName1);  % to get the extension (includes ".")
    formatstr = [];
else
    % series of TIFF files
    
    if strcmp(FileName1, FileName2)
        % Both file names are the same -- user has selected only one file
        disp('Only one (single-page) file selected!  Leaving most outputs blank.');
        pause(1);
        fbase = FileName1;
        frmin = [];
        frmax = [];
        formatstr = [];
        [pathstr, name, ext] = fileparts(FileName1);  % to get the extension (includes ".")
    else
        % the usual case
        [pathstr, name, ext] = fileparts(FileName1);  % to get the extension (includes ".")
        Nch1 = length(FileName1);  % number of characters in file name
        
        samelength = (length(FileName1)== length(FileName2));
        % determine the first character at which the two file names differ
        startk=1;
        while(strcmp(FileName1(startk),FileName2(startk)))
            startk = startk+1;
        end
        fbase = FileName1(1:(startk-1));
        if samelength
            % starting frame
            frmin = str2num(FileName1(startk:Nch1-4));  %#ok<ST2NM> % (last 4 characters are .tif, or whatever extension)
            % ending frame
            frmax = str2num(FileName2(startk:Nch1-4));  %#ok<ST2NM> % (last 4 characters are .tif, or whatever extension)
            % Number of frames (total; can skip later)
            % Nframes = frmax - frmin + 1;
            digits = ceil(log10(frmax+1));
            formatstr = strcat('%0', num2str(digits), 'd');
        else
            % File names are not the same length
            disp('Warning -- file names are not the same length!  Use dialog box:"');
            % A ridiculously convoluted way to extract the numbers at the
            % end of the first and last file names (since I can't get
            % something like a = sscanf('sdf72', '%*s%d') to work:
            a = sscanf(fliplr(FileName1(1:Nch1-4)), '%d%*s');
            frmin = str2num(fliplr(num2str(a)));
            a = sscanf(fliplr(FileName2(1:length(FileName2)-4)), '%d%*s');
            frmax = str2num(fliplr(num2str(a)));
            prompt = {'base file name', 'Min frame no.', 'Max. frame no.', 'Format string'};
            dlg_title = 'File name analysis'; num_lines= 1;
            def     = {fbase, num2str(frmin), num2str(frmax), '%d'};  % default values
            answer  = inputdlg(prompt,dlg_title,num_lines,def);
            fbase = char(answer(1));
            frmin = str2double(answer(2));
            frmax = str2double(answer(3));
            formatstr = char(answer(4));
        end
    end
    
    
end
