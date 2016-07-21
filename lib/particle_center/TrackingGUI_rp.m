% TrackingGUI_rp.m
% 
% Copyright 2012, Raghuveer Parthasarathy, The University of Oregon
%
%%
% Disclaimer / License  
%   This program is free software: you can redistribute it and/or 
%     modify it under the terms of the GNU General Public License as 
%     published by the Free Software Foundation, either version 3 of the 
%     License, or (at your option) any later version.
%   This set of programs is distributed in the hope that it will be useful, 
%   but WITHOUT ANY WARRANTY; without even the implied warranty of 
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
%   General Public License for more details.
%   You should have received a copy of the GNU General Public License 
%   (gpl.txt) along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% TrackingGUI_rp.m: 
% A graphical interface for particle tracking, allowing 
%   -- choice of tracking algorithm
%   -- linkage of objects -> tracks
%   -- display of tracking output overlayed on image frames
%
% Inputs: 
%    im : (optional) 3D array of 2D images, e.g. from TIFFseries.m
%
% Outputs:
%    *None* (!).  It's very difficult to get output from a GUI.  (See various web
%    pages.  Will just save arrays and parameters to a MAT file:
%       objsize : object size parameter used for particle localization
%       thresh : threshold parameter used for particle localization
%       objs : object matrix, output by im2obj_rp.m (see file)
%       objs_link : linked object matrix, output by nnlink_rp.m (see file)
%    Or: save 'simple' output as a text file -- just x and y positions of
%       each linked object, without additional information.
%       Row 1 = x positions (px) of object #1 in each frame in which it
%       exists; Row 2 = y positions.
%       Row 3 = x positions of object #2; Row 4 = y positions. 
%       etc.
%       Note that the columns correspond to frames, but the set of
%       frames in which each object exists need not be the same. 
%       The tab-delimited text file can be opened e.g. by WordPad or Excel
%
% TrackingGUI_rp.m  begun April 3, 2012.  
%
% The display and transparency functions as well as the object information
% selection / data table are taken from WaterGUI[1-6].m by Matt Jemielita, 2010
%
% Raghuveer Parthasarathy
% Modifications since June 2012:
% June 10, 2012 (added text  output option)
% June 27, 2012 (fixed bug for frame slider range if single image is input;
%                consistent definition of threshold default values)
% Last modified June 28, 2012


function TrackingGUI_rp(im)

% A nested function

%% Get Filenames and info for loading images
% -----------------------------------------------------------------------
programdir = cd; % present directory.

if ~exist('im', 'var') || isempty(im)
    % images were not input; get file name info
    imagesloaded = false;
    %    Can be multipage TIFF
    [fbase, frmin, frmax, formatstr, FileName1, FileName2, PathName1 ext ismultipage] = ...
        getnumfilelist;
    % If there's only one image, getnumfilelist returns empty arrays for
    % frmin, frmax, formatstr.  Replace frmin and frmax with 1s.
    if isempty(frmin)
        frmin = 1;
        frmax = 1;
    end
else
    % images were input
    imagesloaded = true;
    frmin = 1;
    frmax = size(im,3);
    FileName1 = '[input 3D array]';
end
Nframes = frmax - frmin + 1;

% -----------------------------------------------------------------------
%% GUI components
% -----------------------------------------------------------------------

%  Initialize and hide the GUI as it is being constructed.
scrsz = get(0,'ScreenSize');  % screen size

fGUI = figure('Name','TrackingGUI_rp', 'Menubar','none', ...
    'Visible','off','Position',[1 1 scrsz(3) scrsz(4)], ...
    'Color', [0 40 90]/255);  
            
% Construct the components.  Create axes, for the images
% Initialization
% Initialize file uploading controls. All these controls will be contained
% in the hupdate panel.
% hfileupdate controls what files will be analyzed by the program

% Signature
uicontrol('Style','text','BackgroundColor', [0.5 0.7 0.7], ...
    'String','TrackingGUI_rp.m.  Raghuveer Parthasarathy, 2012', ...
    'FontWeight', 'bold', 'Units', 'normalized',...
    'HorizontalAlignment', 'left', 'Position',...
    [0.67 0.96 0.25 0.03]);
%Create an exit button
hexit = uicontrol('Style','pushbutton',...
    'String','Exit','Units', 'normalized','FontWeight', 'bold',...
    'Position',[0.94 0.85 0.05 0.1], 'Callback',{@exit_Callback});

%% Panel for loading the files that will be analyzed by the GUI 
himageselect = uipanel('Title', 'Display Frame', 'FontSize', 11, ...
    'Units', 'normalized', 'Position', [0.67 0.85 0.25 0.1]);
   % Primary frame no. to load, analyze
% hFileNameText = 
uicontrol('Parent', himageselect, ...
    'Style', 'text', 'String', FileName1,...
    'FontAngle', 'italic', 'Units', 'normalized',...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.55 0.9 0.4]);
% hframenotextinit = 
uicontrol('Parent', himageselect,'Style','text','Units',...
    'normalized','String','Frame no.', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.1 0.2 0.4]);
hframenotext  = uicontrol('Parent', himageselect,'Style','edit','Units', ...
    'normalized', 'Position',[0.35, 0.1, 0.15, 0.4], ...
    'Callback',{@framenotext_Callback});
hframeno = uicontrol('Parent', himageselect,'Style','slider', ...
    'Max', frmax + 0.1*(Nframes==1), 'Min', frmin, 'Value', frmin, 'Units', 'normalized',...
    'SliderStep', [1/Nframes min([0.1, 10/Nframes])], 'Position', ...
    [0.55 0.1 0.35 0.4], 'Callback',{@frameno_Callback});

%% Tracking parameters
htrackparampanel = uipanel('Title','Parameters','FontSize',11,...
'Units',  'normalized', 'Position',[.67 .49 .32 .34]);
defaultobjsize = 7;
uicontrol('Parent', htrackparampanel,'Style','text','Units',...
    'normalized','String','bpfiltsize', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.85 0.15 0.1]);
hbpfiltsizetext  = uicontrol('Parent', htrackparampanel,'Style','edit','Units', ...
    'normalized', 'String', sprintf('%d', defaultobjsize), 'Position',[0.20, 0.85, 0.1, 0.1], ...
    'Callback',{@bpfiltsizetext_Callback});
hbpfiltsize = uicontrol('Parent', htrackparampanel,'Style','slider', ...
    'Max', 100, 'Min', 0, 'Value', defaultobjsize, 'Units', 'normalized',...
    'SliderStep', [0.01 0.1], 'Position', [0.35 0.85 0.5 0.1], ...
    'Callback',{@bpfiltsize_Callback});
uicontrol('Parent', htrackparampanel,'Style','text','Units',...
    'normalized','String','nsize', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.70 0.15 0.1]);
hnsizetext  = uicontrol('Parent', htrackparampanel,'Style','edit','Units', ...
    'normalized', 'String', sprintf('%d', defaultobjsize), 'Position',[0.20, 0.70, 0.1, 0.1], ...
    'Callback',{@nsizetext_Callback});
hnsize = uicontrol('Parent', htrackparampanel,'Style','slider', ...
    'Max', 100, 'Min', 1, 'Value', defaultobjsize, 'Units', 'normalized',...
    'SliderStep', [0.01 0.1], 'Position', [0.35 0.70 0.5 0.1], ...
    'Callback',{@nsize_Callback});
% Button to lock filtering object size and neighborood object size together
hlockobjsize  = uicontrol('Parent', htrackparampanel,'Style','toggle',...
    'String','Lock','Units', 'normalized', 'BackgroundColor', [0.8 1.0 0.6], ...
    'Position',[0.88, 0.7,0.1,0.2],'Value', 1, ...
    'Callback',{@lockobjsize_Callback});

% hthreshoptinit = 
uicontrol('Parent', htrackparampanel,'Style','text','Units',...
    'normalized','String','Threshold Option', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.5 0.3 0.12]);
hthreshopt = uicontrol('Parent', htrackparampanel, 'Style','popupmenu', ...
    'String', '1. Intens. Thresh.|2. Std. Thresh.|3. Bright N', ...
    'Units', 'normalized','Position', [0.05,0.37,0.3,0.12], ...
    'Callback',{@threshopt_Callback});

%hthreshtextinit1 = 
uicontrol('Parent', htrackparampanel, 'Style','text',...
    'String','thr (0-1)', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.38,0.52,0.13,0.12]);
defaultthresh1 = 0.99;
hthresh1 = uicontrol('Parent', htrackparampanel, 'Style','slider', ...
    'Max', 0.99999, 'Min', 0.00, ...
    'SliderStep', [0.02 0.1], 'Units', 'normalized',...
    'Value', defaultthresh1, 'Position', [0.7,0.52,0.27,0.12], ...
    'Callback',{@thresh1_Callback});
hthreshtext1  = uicontrol('Parent', htrackparampanel, 'Style','edit',...
    'Units', 'normalized', 'String', sprintf('%.4f', get(hthresh1, 'Value')), ...
    'Position',[0.53,0.52,0.15,0.12], ...
    'Callback',{@threshtext1_Callback});
%hthreshtextinit2 = 
uicontrol('Parent', htrackparampanel, 'Style','text',...
    'String','thr (>1)', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.38,0.37,0.13,0.12]);
defaultthresh2 = 3.0;
hthreshtext2  = uicontrol('Parent', htrackparampanel, 'Style','edit',...
    'String', sprintf('%.1f', defaultthresh2), 'Units', 'normalized', ...
    'Position',[0.53,0.37,0.15,0.12], ...
    'Callback',{@threshtext2_Callback});
hthresh2 = uicontrol('Parent', htrackparampanel, 'Style','slider', ...
    'Max', 100, 'Min', 1, 'Value', defaultthresh2, ...
    'SliderStep', [0.01 0.10], 'Units', 'normalized',...
    'Position', [0.7,0.37,0.27,0.12], ...
    'Callback',{@thresh2_Callback});
defaultthresh3 = 3;
%hthreshtextinit3 = 
uicontrol('Parent', htrackparampanel, 'Style','text',...
    'String','N (>=1)', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.38,0.22,0.13,0.12]);
hthreshtext3  = uicontrol('Parent', htrackparampanel, 'Style','edit',...
    'String', sprintf('%d', defaultthresh3), 'Units', 'normalized', ...
    'Position',[0.53,0.22,0.15,0.12], ...
    'Callback',{@threshtext3_Callback});
hthresh3 = uicontrol('Parent', htrackparampanel, 'Style','slider', ...
    'Max', 100, 'Min', 1, 'Value', defaultthresh3, ...
    'SliderStep', [0.01 0.10], 'Units', 'normalized',...
    'Position', [0.7,0.22,0.27,0.12], ...
    'Callback',{@thresh3_Callback});
%hfitstrtextinit = 
uicontrol('Parent', htrackparampanel, 'Style','text',...
    'String','Fit method', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.05,0.02,0.2,0.12]);
fitstrarray = {'radial', 'gaussmle', 'nonlineargauss', 'lineargauss','centroid'};
hfitstr = uicontrol('Parent', htrackparampanel, 'Style','popupmenu', ...
    'String', strcat(char(fitstrarray(1)), ' | ', ...
                     char(fitstrarray(2)), ' | ', ...
                     char(fitstrarray(3)), ' | ', ...
                     char(fitstrarray(4)), ' | ', ...
                     char(fitstrarray(5))), ...
    'Units', 'normalized','Position', [0.35,0.02,0.25,0.12],...
    'Value', 1, 'Callback',{@fitstr_Callback});
%     'String', 'radial | Gauss.MLE | nonlineargauss | lineargauss | centroid', ...

%% Displaying filtered and thresholded images
hdispprocesspanel = uipanel('Title','Display Processed Imgs','FontSize',11,...
'Units',  'normalized', 'Position',[.67 .37 .32 .1]);
% hdispimgText = 
uicontrol('Parent', hdispprocesspanel, ...
    'Style', 'text', 'String', 'Filtered for finding maxima, not for fitting',...
    'FontAngle', 'italic', 'Units', 'normalized',...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.7 0.9 0.3]);
hdispfilt  = uicontrol('Parent', hdispprocesspanel,'Style','checkbox',...
    'String','DispFilt','Units', 'normalized', ...
    'Position',[0.05, 0.2,0.3,0.4],'Value', 0,...
    'Callback',{@dispfilt_Callback});
hdispthreshfilt    = uicontrol('Parent', hdispprocesspanel,'Style','checkbox',...
    'String','Dispthreshfilt','Units', 'normalized', ...
    'Position',[0.4, 0.2,0.3,0.4],'Value', 0,...
    'Callback',{@dispthreshfilt_Callback});

%% Tracking
htrackpanel = uipanel('Title','Tracking','FontSize',11,...
'Units',  'normalized', 'Position',[.67 .14 .32 .20]);
htrackthisfr    = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Track this frame','Units', 'normalized',...
    'Position',[0.05,0.75,0.3,0.1],'value', 0, ...
    'Callback',{@trackthisfr_Callback});
htrackthisdone    = uicontrol('Parent', htrackpanel,'Style','checkbox',...
    'String','Done','Units', 'normalized',...
    'Position',[0.4,0.75,0.3,0.1],'value', 0, ...
    'Callback',{@trackthisdone_Callback});
htrackallfr    = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Track all frames','Units', 'normalized',...
    'Position',[0.05,0.55,0.3,0.1],'value', 0, ...
    'Callback',{@trackallfr_Callback});
htrackalldone    = uicontrol('Parent', htrackpanel,'Style','checkbox',...
    'String','Done','Units', 'normalized',...
    'Position',[0.4,0.55,0.3,0.1],'value', 0, ...
    'Callback',{@trackalldone_Callback});
hcleartracking    = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Clear','Units', 'normalized',...
    'Position',[0.75,0.55,0.25,0.1],'value', 0, ...
    'Callback',{@cleartracking_Callback});
hlink    = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Link Tracks','Units', 'normalized',...
    'Position',[0.05,0.05,0.3,0.1],'value', 0, ...
    'Callback',{@link_Callback});
%hlinksteptextinit = 
uicontrol('Parent', htrackpanel,'Style','text','Units',...
    'normalized','String','Link MaxStep^2', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.3 0.25 0.12]);
hlinksteptext  = uicontrol('Parent', htrackpanel,'Style','edit','Units', ...
    'normalized', 'Position',[0.35, 0.3, 0.1, 0.12], 'String', '100', ...
    'Callback',{@linksteptext_Callback});
%hlinkmemorytextinit = 
uicontrol('Parent', htrackpanel,'Style','text','Units',...
    'normalized','String','Link Memory', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.5 0.3 0.2 0.12]);
hlinkmemorytext  = uicontrol('Parent', htrackpanel,'Style','edit','Units', ...
    'normalized', 'Position',[0.75, 0.3, 0.1, 0.12], 'String', '1', ...
    'Callback',{@linkmemorytext_Callback});
hlinkdone    = uicontrol('Parent', htrackpanel,'Style','checkbox',...
    'String','Done','Units', 'normalized',...
    'Position',[0.4,0.05,0.3,0.1],'value', 0, ...
    'Callback',{@linkdone_Callback});
hclearlink = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Clear','Units', 'normalized',...
    'Position',[0.75,0.05,0.25,0.1],'value', 0, ...
    'Callback',{@clearlink_Callback});

%% Save output, in a MAT file, or load previously calculated values
hsavepanel = uipanel('Title','Save / Load','FontSize',11,...
'Units',  'normalized', 'Position',[.67 .01 .32 .10]);
uicontrol('Parent', hsavepanel,'Style','pushbutton',...
    'String','Save output','Units', 'normalized',...
    'Position',[0.05,0.1,0.25,0.8],'value', 0, ...
    'Callback',{@saveoutput_Callback});
uicontrol('Parent', hsavepanel,'Style','pushbutton',...
    'String','Save TXT','Units', 'normalized',...
    'Position',[0.35,0.1,0.25,0.8],'value', 0, ...
    'Callback',{@savetxt_Callback});
uicontrol('Parent', hsavepanel,'Style','pushbutton',...
    'String','Load obj data','Units', 'normalized',...
    'Position',[0.65,0.1,0.25,0.8],'value', 0, ...
    'Callback',{@loadoutput_Callback});

%% Messages
hmsg    = uicontrol('Style','text','BackgroundColor', [1 1 0.7],...
    'String','Message: ', 'FontWeight', 'bold', 'Units', 'normalized',...
    'HorizontalAlignment', 'left', 'Position', [0.5 0.01 0.15 0.12]);

%% Display tracks, etc.

hdisptrackspanel = uipanel('Title','Display Tracks','FontSize',11,...
'Units',  'normalized', 'Position',[.35 .01 .09 .2]);
hdispcircles    = uicontrol('Parent', hdisptrackspanel,'Style','checkbox',...
    'String','Show objs','Units', 'normalized', ...
    'Position',[0.05, 0.6,0.8,0.25],'Value', 0,...
    'Callback',{@displayimage});
hdispIDs    = uicontrol('Parent', hdisptrackspanel,'Style','checkbox',...
    'String','Show IDs','Units', 'normalized', ...
    'Position',[0.05, 0.3,0.8,0.25],'Value', 0,...
    'Callback',{@displayimage});
hdisptracks    = uicontrol('Parent', hdisptrackspanel,'Style','checkbox',...
    'String','Show Tracks','Units', 'normalized', ...
    'Position',[0.05, 0.01,0.8,0.25],'Value', 0,...
    'Callback',{@displayimage});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the button group which allows the user to either collect
% information about particular cells in the image or to have a zoom feature
% for the image shown.
hbgroup = uibuttongroup('visible','off','Position',[0.23 0.01 .1 .1]);
% Create three radio buttons in the button group.  (Could name the handles
% to the controls, but not necessary)
uicontrol('Style','Radio','String','Zoom','Units', 'normalized',...
    'pos',[.10 .05 .90 .25],'Parent',hbgroup, ...
    'HandleVisibility','off');
uicontrol('Style','Radio','String','Pan','Units',...
    'normalized', 'pos',[.10 .35 .900 .25],'Parent',hbgroup,...
    'HandleVisibility','off');
uicontrol('Style','Radio','String','Object Information','Units',...
    'normalized', 'pos',[.10 .65 .900 .25],'Parent',hbgroup,...
    'HandleVisibility','off');
% Initialize some button group properties. 
set(hbgroup,'SelectionChangeFcn',@selcbk);
set(hbgroup,'SelectedObject',[]);  % No selection
set(hbgroup,'Visible','on');
%Initializing the datacursormode
dcm_obj = datacursormode(fGUI);
%Initially the datacursormode() will be set to off. By clicking on the
%object information button it will be turned on.
set(dcm_obj,'Enable','off');
set(dcm_obj, 'Updatefcn', @ObjectInfo_Callback);

% ------------------------------------------------------------------
%% Initialize the GUI.

% ------------------------------------------------------------------
% Create all variables here

% Outputs
% NOT output in the function (very hard to do this in a GUI).  After 
% calculation, click 'save' to write to MAT file
bpfiltsize = get(hbpfiltsize, 'Value');
nsize = get(hnsize, 'Value');
lockobjsize = get(hlockobjsize, 'Value');
thresh = get(hthresh1, 'Value');
threshopt = [];
objs = [];
objs_link = [];

% Other
fitstr = char(fitstrarray(get(hfitstr,'Value')));
lsqoptions = optimset('lsqnonlin');
linkstep = round(str2double(get(hlinksteptext,'string')));
linkmem =  round(str2double(get(hlinkmemorytext,'string')));

% Load one image, or use the first frame of the 'im' array
% In general, will use 'A' for the displayed image
if imagesloaded
    A = im(:,:,1);
    iscolor = false;
else
    cd(PathName1) %Go to the directory of the images
    if ~ismultipage
        A=imread(FileName1);    % load the first file image
    else
        % multipage TIFF
        A = imread(FileName1, 1);
    end
    iscolor = (ndims(A)==3);  % TRUE if the image is color (3 layers)
    if iscolor
        prompt = {'If color: which channel to use?  (1, 2, 3): '};
        dlg_title = 'Color option'; num_lines= 1;
        % guess at which channel to use as default, based on total brightness
        b = zeros(1,3);
        for ch = 1:3
            b(ch) = sum(sum(A(:,:,ch)));
        end
        [mb, ic] = max(b);
        def     = {num2str(ic)};  % default values
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        colchannel = round(str2double(answer(1)));
        A = A(:,:,colchannel);  % select the color channel to examine
        pause(0.5);  % seems to be necessary to keep MATLAB from crashing!
    end
end

currframe = frmin;  % # of the current (primary) frame

% Show filtered, thresholded image
showfilt = false;
showthreshfilt = false;
filtA = [];  % filtered image
threshfiltA = [];  % thresholded, filtered image

cd(programdir); % Go back to the directory from which the GUI was called

% Handles to display images
% Am I going to use these?
im1 = [];  % handle to primary image display

% Status variables
trackthisdone = false(Nframes,1);  % true for frames that have been segmented
islinkdone = false;  % is linkage of objects into tracks done?

% necessary?
colormap('gray');
        
outfile = [];  % file name for saving results (MAT file)

%Creates axes on the figure. At a later time in the program these axes will
%be filled with the images we are analyzing.
axeshandle = axes('Position', [0.01 0.24 0.65 0.75]);

%Creating a table where the information about the cell we click on is
%displayed.
selectedobjectdata = cell(6,1);
t = uitable('Parent', fGUI, 'Units', 'normalized','Position',...
   [0.01 0.01 0.2 0.2]);
set(t, 'ColumnName', []);
set(t, 'RowName',{ 'Centroid x',  'Centroid y', 'Brightness',...
    'Particle ID', 'Frame no.', 'Track ID'});
set(t, 'Data', selectedobjectdata);
align([t axeshandle], 'Fixed',100,  'Top');

% Move the GUI to the center of the screen.
movegui(fGUI,'center')
% Make the GUI visible.
set(fGUI,'Visible','on');
set([htrackthisfr, htrackallfr,hcleartracking, hlink, hclearlink, ...
    hexit, hmsg], 'BackgroundColor',[0.85 1.0 0.6]);
set([hexit, hmsg], 'BackgroundColor',[1.0 0.5 0.2]);

figure(fGUI);

% Set sliders, etc.
set(hframenotext, 'String', sprintf('%d', currframe));
set(hframeno, 'Value', currframe);
% Call the threshold option function, just to highlight the default row
threshopt_Callback;

% Display the first image
displayimage

%% ----------------------------------------------------------------------
% **********************************************************************

% Callbacks and other functions

    function exit_Callback(source,eventdata)
        % Exit
        close all
    end

% Callback functions for loading and displaying frames

    function framenotext_Callback(source,eventdata)
        % set frame to view (text entry)
        currframe = round(str2double(get(source,'string')));
        % Also update slider:
        set(hframeno, 'Value', currframe);
        % Load the image
        loadimages;
        % Display
        displayimage;
        % Update the checkbox indicating whether segmentation has been done
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
    end

    function frameno_Callback(hObject, source,eventdata)
        % set primary frame to view (slider)
        currframe = round(get(hObject,'Value'));
        % Also update text entry box:
        set(hframenotext, 'String', sprintf('%d', currframe));
        % Load the image
        loadimages;
        % Display
        displayimage;
        % Update the checkbox indicating whether segmentation has been done
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
    end

    function loadimages(source,eventdata)
        set(hmsg, 'String', 'Loading images...');
        if imagesloaded
            % all images are in the 'im' array
            A = im(:,:,currframe);
        else
            % load from file
            cd(PathName1)     %Go to the directory of the images
            % load the image into "A"
            % Don't do any calculations or updates
            % Note that all file names were determined previously
            if ~ismultipage
                framestr = sprintf(formatstr, currframe);
                A  = imread(strcat(fbase, framestr, ext));
            else
                % multipage TIFF
                A  = imread(strcat(fbase, ext), currframe);
            end
            cd(programdir) %go back to the original directory
            if iscolor
                % select the color channel to examine
                A = A(:,:,colchannel);
            end
        end
        set(hmsg, 'String', 'Image loaded.');
        % if the boxes are checked, calculate the filtered and thresholded
        % images
        if showfilt || showthreshfilt
            filtA = calcfilt;
        end
        if showthreshfilt
            threshfiltA = calcthreshfilt;
        end
    end

% Object size and threshold parameters

    function bpfiltsizetext_Callback(source,eventdata)
        % set bpfiltsize, for spatial filtering (text entry)
        bpfiltsize = round(str2double(get(source,'string')));
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
        set(hbpfiltsize, 'Value', bpfiltsize);
        % if filter and neighborhood values are locked together, update
        if lockobjsize
            nsize = bpfiltsize;
            set(hnsize, 'Value', nsize);
            set(hnsizetext, 'String', sprintf('%d', nsize));
        end
    end

    function bpfiltsize_Callback(hObject, source,eventdata)
        bpfiltsize = round(get(hObject,'Value'));
        % set bpfiltsize, for spatial filtering (text entry)
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
        % Also update text entry box:
        set(hbpfiltsizetext, 'String', sprintf('%d', bpfiltsize));
        % if filter and neighborhood values are locked together, update
        if lockobjsize
            nsize = bpfiltsize;
            set(hnsize, 'Value', nsize);
            set(hnsizetext, 'String', sprintf('%d', nsize));
        end
    end

    function nsizetext_Callback(source,eventdata)
        % set neighborhood objsize (text entry)
        nsize = round(str2double(get(source,'string')));
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
        % update the slider
        set(hnsize, 'Value', nsize);
        % if filter and neighborhood values are locked together, update
        if lockobjsize
            bpfiltsize = nsize;
            set(hbpfiltsize, 'Value', bpfiltsize);
            set(hbpfiltsizetext, 'String', sprintf('%d', nsize));
        end
    end

    function nsize_Callback(hObject, source,eventdata)
        nsize = round(get(hObject,'Value'));
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
        % Also update text entry box:
        set(hnsizetext, 'String', sprintf('%d', nsize));
        % if filter and neighborhood values are locked together, update
        if lockobjsize
            bpfiltsize = nsize;
            set(hbpfiltsize, 'Value', bpfiltsize);
            set(hbpfiltsizetext, 'String', sprintf('%d', nsize));
        end
    end

    function lockobjsize_Callback(hObject, source,eventdata)
        % make the filtering and neighborhood object sizes the same
        lockobjsize = get(hlockobjsize, 'Value');
        if lockobjsize
            set(hObject,'BackgroundColor', [0.8 1.0 0.6]);
            nsize = bpfiltsize;
            set(hnsize, 'Value', nsize);
            set(hnsizetext, 'String', sprintf('%d', nsize));
            set(hObject,'BackgroundColor', [0.8 1.0 0.6]);
            % update the filtered images, if these are being displayed
            updatefiltandthresh;
        else
            set(hObject,'BackgroundColor', [1.0 0.5 0.3]);
        end
    end

    function whichthresh
        % determine which threshold value to use, based on thresholding
        % option
        threshopt = get(hthreshopt,'Value');
        switch threshopt
            case 1
                thresh = get(hthresh1,'Value');
            case 2
                thresh = -1.0*get(hthresh2,'Value');
                % fo4_rp.m interprets negative thresholds as "option 2"
                % inputs.
            case 3
                thresh = round(get(hthresh3,'Value'));
        end
    end
        
    function threshopt_Callback(source,eventdata)
        % set thresholding option
        threshopt = get(hthreshopt,'Value');
        ltgray = 0.95*[1 1 1];
        set(hthreshtext1,'BackgroundColor',ltgray);
        set(hthresh1,'BackgroundColor',ltgray);
        set(hthreshtext2,'BackgroundColor',ltgray);
        set(hthresh2,'BackgroundColor',ltgray);
        set(hthreshtext3,'BackgroundColor',ltgray);
        set(hthresh3,'BackgroundColor',ltgray);
        switch threshopt
            case 1
                set(hthreshtext1,'BackgroundColor',[1 1 0.7]);
                set(hthresh1,'BackgroundColor',[1 1 0.7]);
            case 2
                set(hthreshtext2,'BackgroundColor',[1 1 0.7]);
                set(hthresh2,'BackgroundColor',[1 1 0.7]);
            case 3
                set(hthreshtext3,'BackgroundColor',[1 1 0.7]);
                set(hthresh3,'BackgroundColor',[1 1 0.7]);
        end
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
    end

    function threshtext1_Callback(source,eventdata)
        % set thresh (text entry)
        thresh1 = str2double(get(source,'string'));
        set(hthresh1, 'Value', thresh1);
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
    end
    function threshtext2_Callback(source,eventdata)
        % set thresh (text entry)
        thresh2 = str2double(get(source,'string'));
        set(hthresh2, 'Value', thresh2);
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
    end
    function threshtext3_Callback(source,eventdata)
        % set thresh (text entry)
        thresh3 = round(str2double(get(source,'string')));
        set(hthresh3, 'Value', thresh3);
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
    end

    function thresh1_Callback(hObject, source,eventdata)
        thresh1 = get(hObject,'Value');
        % Also update text entry box:
        set(hthreshtext1, 'String', sprintf('%.4f', thresh1));
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
    end
    function thresh2_Callback(hObject, source,eventdata)
        thresh2 = get(hObject,'Value');
        % Also update text entry box:
        set(hthreshtext2, 'String', sprintf('%.4f', thresh2));
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
    end
    function thresh3_Callback(hObject, source,eventdata)
        thresh3 = round(get(hObject,'Value'));
        % Also update text entry box:
        set(hthreshtext3, 'String', sprintf('%.d', thresh3));
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updatefiltandthresh;
    end

    function fitstr_Callback(hObject, source,eventdata)
        fitstr = char(fitstrarray(get(hObject,'Value')));
    end

    % Callbacks for Displaying filtered and post-threshold images
    
    function filtA = calcfilt
        % calculate the bandpass filtered image, if the filtersize > 0
        if bpfiltsize>0
            set(hmsg, 'String', 'Calling bpass.m to filter');
            filtA = bpass(double(A),1,bpfiltsize);  % use present image
        else
            set(hmsg, 'String', 'filter size = 0; not filtering');
            filtA = double(A);
        end
    end

    function threshfiltA = calcthreshfilt
        % calculate the dilated image of local maxima that pass the
        % threshold.  
        % Calculate filtered image even if previously calculated, in case
        % objsize (bpfiltsize) has changed.
        filtA = calcfilt;  % uses bpfiltsize for filtering size
        [y, x] = calcthreshpts(filtA, get(hthreshopt,'Value'), thresh, nsize);
        threshfiltAmask = false(size(filtA));
        threshfiltAmask(sub2ind(size(filtA), round(y), round(x)))=true;
        threshfiltAmask = imdilate(threshfiltAmask, strel('disk', floor(nsize/2)));
        threshfiltA = filtA.*threshfiltAmask;
    end

    function updatefiltandthresh
        % calls functions to recalculate filtered and post-threshold
        % images, and display
        if showfilt || showthreshfilt
            filtA = calcfilt;
        end
        if showthreshfilt
            threshfiltA = calcthreshfilt;
        end
        displayimage
    end

    function dispfilt_Callback(hObject, source,eventdata)
        showfilt = get(hObject, 'Value');
        if showfilt
            filtA = calcfilt;  % calculate, even if filtA isn't empty, in case filter size has changed
            set(hdispthreshfilt, 'Value', false);  % turn off the other checkbox
        end
        displayimage;
    end

    function dispthreshfilt_Callback(hObject, source,eventdata)
        showthreshfilt = get(hObject, 'Value');
        % Could check to see if threshfiltA has already been calculated,
        % but should be careful that threshold parameter values haven't
        % changed.
        if showthreshfilt
            filtA = calcfilt;
            threshfiltA = calcthreshfilt;
            set(hdispfilt, 'Value', false);  % turn off the other checkbox
        end
        displayimage;
    end

    % Callback functions for Tracking (particle localization)
             
    function [] = trackthisfr_Callback(hObject, eventdata, handles)
        % Track a single frame (find objects)
        set(hmsg, 'String', 'Tracking started...'); 
        pause(0.1);  % does this help the display issue?
        tmpobj = fo4_rp(A, [bpfiltsize nsize], thresh, fitstr, lsqoptions);
        tmpobj(5,:) = currframe-frmin+1;
        objs = [objs tmpobj];
        trackthisdone(currframe-frmin+1) = true;  % note that tracking has been done
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
        set(hmsg, 'String', 'Tracking completed!');
    end

    function [] = trackthisdone_Callback(hObject, eventdata, handles)
        % Don't do anything if the user clicks the "tracking done" 
        % checkbox, reset its value to whatever the true array value is.
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
    end

    function [] = trackallfr_Callback(hObject, eventdata, handles)
        % Track all frames
        set(hmsg, 'String', 'Tracking of all images started...'); 
        pause(0.1);  % does this help the display issue?
        disp('Tracking of all images started...')

        progtitle = sprintf('TrackingGUI: Tracking using fo4_rp.m ...  '); 
        if (Nframes > 1)
            progbar = waitbar(0, progtitle);  % will display progress
        end
        objs = [];
        oldA = A;
        oldcurrframe = currframe;
        for j = frmin:frmax
            % Loop through all frames
            % similar code as in im2obj_rp.m
            if imagesloaded
                % all images already loaded
                tmpobj = fo4_rp(im(:,:,j), [bpfiltsize nsize], thresh, fitstr, lsqoptions);
            else
                % read from file; % use the variablea "A" and "currframe"
                currframe = j;
                loadimages;
                tmpobj = fo4_rp(A, [bpfiltsize nsize], thresh, fitstr, lsqoptions);
            end
            tmpobj(5,:) = j-frmin+1;
            objs = [objs tmpobj];
            % show progress -- not called if just one frame
            if mod(j-frmin+1,10)==0
                waitbar((j-frmin+1)/Nframes, progbar, ...
                    strcat(progtitle, sprintf('frame %d of %d', (j-frmin+1), Nframes)));
            end
            trackthisdone(j-frmin+1) = true;  % note that tracking has been done
        end
        A = oldA;
        currframe = oldcurrframe;
        if Nframes>1
            close(progbar)
        end
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
        set(htrackalldone, 'Value', sum(trackthisdone)==length(trackthisdone));
        set(hmsg, 'String', 'Tracking of all frames completed!');
        disp('Tracking of all frames completed!')
    end

    function [] = trackalldone_Callback(hObject, eventdata, handles)
        % Don't do anything if the user clicks the  
        % checkbox, reset its value to what it should be
        set(htrackalldone, 'Value', sum(trackthisdone)==length(trackthisdone));
    end

    function [] = cleartracking_Callback(hObject, eventdata, handles)
        % Clear all the tracking output, and the linking output
        objs = [];
        objs_link = [];
        trackthisdone = false(Nframes,1);
        islinkdone = false;
        set(htrackalldone, 'Value', sum(trackthisdone)==length(trackthisdone));
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
        set(hlinkdone, 'Value', false);
        set(hmsg, 'String', 'Cleared tracking output for all frames');
    end

% Parameters for linking

    function  [] = link_Callback(hObject, eventdata, handles)
        
        % Linking objects into trajectories, using nnlink_rp
        set(hmsg, 'String', 'Linking tracks ...');
        disp('Linking tracks ...')
        
        if sum(trackthisdone)==length(trackthisdone)
            objs_link = nnlink_rp(objs, linkstep, linkmem, true);
            set(hmsg, 'String', 'Linking done ...');
            islinkdone = true;
        else
            % All frames need to be tracked before linking!
            set(hmsg, 'String', 'All frames need to be tracked before linking!');
        end
        set(hlinkdone, 'Value', islinkdone);
    end

    function [] = linkdone_Callback(hObject, eventdata, handles)
        % Don't do anything if the user clicks the "link done" 
        % checkbox, reset its value to whatever the true array value is.
        set(hlinkdone, 'Value', islinkdone);
    end

    function linksteptext_Callback(source,eventdata)
        % set max step size for linkage (text entry)
        linkstep = round(str2double(get(source,'string')));
    end

    function linkmemorytext_Callback(source,eventdata)
        % set memory size for linkage (text entry)
        linkmem = round(str2double(get(source,'string')));
    end

    function [] = clearlink_Callback(hObject, eventdata, handles)
        % Clear the linkages; keep objs matrix
        objs_link = [];
        islinkdone = false;
        set(hlinkdone, 'Value', islinkdone);
        set(hmsg, 'String', 'Cleared linkage of tracks');
    end

%% Save output, or load previous calc. values
    function [] = saveoutput_Callback(hObject, eventdata, handles)
        % Save results in a MAT file
        % Dialog box for name
        prompt = {'Output File Name (*include* ".mat")'};
        dlg_title = 'Save output'; num_lines= 1;
        if isempty(outfile)
            def     = {'TrackGUIoutput.mat'};  % default value
        else
            def     = {outfile};  % default value
        end
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        outfile = char(answer(1));
        if ~isempty(outfile)
            set(hmsg, 'String', 'Starting to save variables...');
            save(outfile, 'Nframes', 'objs', 'objs_link', 'threshopt', ...
                'thresh', 'bpfiltsize', 'nsize', 'trackthisdone', 'islinkdone');  % Save these variables
            set(hmsg, 'String', 'Done saving variables.');
        end
    end

    function [] = savetxt_Callback(hObject, eventdata, handles)
        % Save results in a simple text file; positions only
        % Dialog box for name
        if islinkdone
            prompt = {'Output File Name for simple Text Output (*include* extension)'};
            dlg_title = 'Save output'; num_lines= 1;
            def     = {'TrackGUIoutputTXT.txt'};  % default value
            answer  = inputdlg(prompt,dlg_title,num_lines,def);
            outfile = char(answer(1));
            if ~isempty(outfile)
                set(hmsg, 'String', 'Starting to save variables...');
                fo = fopen(outfile, 'w');
                utrk = unique(objs_link(6,:));
                for j=1:length(utrk)
                    xj = objs_link(1,objs_link(6,:)==utrk(j));
                    yj = objs_link(2,objs_link(6,:)==utrk(j));
                    for k=1: length(xj)
                        fprintf(fo, '%.3f\t', xj(k));
                    end
                    fprintf(fo, '\n');
                    for k=1: length(xj)
                        fprintf(fo, '%.3f\t', yj(k));
                    end
                    fprintf(fo, '\n');
                end
                fclose(fo);
                set(hmsg, 'String', 'Done saving variables.');
            end
        else
            beep
            set(hmsg, 'String', 'Linking must precede simple text output.');
        end
    end

    function [] = loadoutput_Callback(hObject, eventdata, handles)
        % Load obj data from previously saved MAT file
        % Dialog box for name
        prompt = {'Input File Name (*include* ".mat")'};
        dlg_title = 'Load output'; num_lines= 1;
        if isempty(outfile)
            def     = {'TrackGUIoutput.mat'};  % default value
        else
            def     = {outfile};  % default value
        end
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        infile = char(answer(1));
        if ~isempty(infile)
            S = load(infile);  % Load variables into a structure
            % -- necessary since this is a nested function
            if S.Nframes ~= Nframes
                set(hmsg, 'String', 'ERROR!  Number of frames is different');
                warndlg('Warning: Number of frames is inconsistent.');
                pause(2)
            else
                % reassign the variables!
                objs = S.objs;
                objs_link = S.objs_link;
                threshopt = S.threshopt;
                thresh = S.thresh;
                if isfield(S, 'objsize')
                    % older version: one objsize variable instead of two
                    nsize = objsize;
                    bpfiltsize = objsize;
                else
                    bpfiltsize = S.bpfiltsize;
                    nsize = S.nsize;
                end
                islinkdone = S.islinkdone;
                trackthisdone = S.trackthisdone;
                % update the filtered images, if these are being displayed
                updatefiltandthresh;
                set(hmsg, 'String', 'Done loading variables.');
                % set sliders and checkboxes
                set(hbpfiltsize, 'Value', bpfiltsize);
                set(hbpfiltsizetext, 'String', sprintf('%d', bpfiltsize));
                set(hnsize, 'Value', nsize);
                set(hnsizetext, 'String', sprintf('%d', nsize));
                set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
                switch threshopt
                    case 1
                        set(hthresh1,'Value', thresh);
                        set(hthreshtext1, 'String', sprintf('%.4f', thresh));
                    case 2
                        set(hthresh2,'Value',-thresh);
                        set(hthreshtext2, 'String', sprintf('%.4f', -thresh));
                    case 3
                        set(hthresh3,'Value',thresh);
                        set(hthreshtext3, 'String', sprintf('%d', thresh));
                end
                set(hthreshopt, 'Value', threshopt);
                threshopt_Callback; % set the shading, etc.; a bit redundant
                set(htrackalldone, 'Value', sum(trackthisdone)==length(trackthisdone));
                set(hlinkdone, 'Value', islinkdone);  % update linkage variable
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Call back functions for image plotting

    function displayimage(hObject, eventdata, handles)
        % Display the present frame
        if exist('im1', 'var')
            delete(im1); clear im1
        end
        if showthreshfilt
            im1 = imshow(threshfiltA, [], 'Parent', axeshandle);
        elseif showfilt
            im1 = imshow(filtA, [], 'Parent', axeshandle);
        else
            im1 = imshow(A, [], 'Parent', axeshandle);
        end
        set(axeshandle, 'Visible', 'off');
        colormap('gray');
        hold on
        
        % display track information on the image, if desired
        if trackthisdone(currframe-frmin+1)
            % tracking has been done
            if islinkdone
                objsthisframe = objs_link(:,objs_link(5,:)==currframe-frmin+1);
            else
                objsthisframe = objs(:,objs(5,:)==currframe-frmin+1);
            end
            if get(hdispcircles, 'Value')
                % plot circles on top, if tracking is done
                plot(objsthisframe(1,:), objsthisframe(2,:), 'o', 'color', [0.3 0.7 0.5])
            end
            if get(hdispIDs, 'Value')
                % show IDs
                if islinkdone
                    for j=1:size(objsthisframe,2)
                        text(objsthisframe(1,j), objsthisframe(2,j),num2str(objsthisframe(6,j)), 'color', [0.7 0.4 0.1])
                    end
                else
                    for j=1:size(objsthisframe,2)
                        text(objsthisframe(1,j), objsthisframe(2,j),num2str(objsthisframe(4,j)), 'color', [0.7 0.4 0.1])
                    end
                end
            end
            if get(hdisptracks, 'Value')
                % plot lines corresponding to the tracks of each object
                if islinkdone
                    % necessary
                    for j=1:size(objsthisframe,2)
                        allx = objs_link(1,objs_link(6,:)==objsthisframe(6,j));
                        ally = objs_link(2,objs_link(6,:)==objsthisframe(6,j));
                        plot(allx, ally, '-', 'color', [0.7 0.3 0.5])
                    end
                end
            end
        end
    end

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function output_txt = ObjectInfo_Callback(obj,eventdata)
        % Fill in the Object Properties table.
        % obj          Currently not used (empty)
        % event_obj    Object containing event data
        % output_txt   Data cursor text (string or cell array of strings)
        pos = get(eventdata,'Position');
        
        % Get the tracking information of the selected object
        if islinkdone
            objsthisframe = objs_link(:,objs_link(5,:)==currframe-frmin+1);
        else
            objsthisframe = objs(:,objs(5,:)==currframe-frmin+1);
        end
        % Euclidean distance matrix -- calculation method from Roland
        % Bunschoten's distance.m on the File Exchange
        aa=sum(pos'.*pos',1); 
        bb=sum([objsthisframe(1,:); objsthisframe(2,:)].*[objsthisframe(1,:); objsthisframe(2,:)],1);  
        d = sqrt(abs(aa( ones(size(bb,2),1), :)' + bb( ones(size(aa,2),1), :) - 2*pos*[objsthisframe(1,:); objsthisframe(2,:)]));
        % [avoid an extra function call] d = distance(pos', [objsthisframe(1,:); objsthisframe(2,:)]);
        [mind, imind] = min(d);
        thisobj = objsthisframe(:,imind);
        if size(thisobj,1)>1
            % not likely, but two closest objects; take the first one (arbitrary)
            thisobj = thisobj(:,1);
        end
        
        % Get the track or particle ID of the selected object
        
        % Properties of this object in this frame
        % When clicking on the object the only information in output_txt is 
        % shown on the figure
        % The rest of the information is sent to the data table.
        if islinkdone
            % Linkage across frames is done
            output_txt = thisobj(6);
        else
            % Not done, so use particle ID
            output_txt = thisobj(4);
        end
        set(t, 'Data', thisobj);  % for data table
    end



%Callback for the button group

    function selcbk(hObject,eventdata)
        switch get(eventdata.NewValue,'String') % Get Tag of selected object.
            case 'Zoom'
                % enable zoom feature.
                zoom(fGUI, 'on');
            case 'Object Information'
                %Enabling data cursor mode-This will allow us to collect data about
                %individual cells in the image.
                zoom(fGUI, 'off');
                set(dcm_obj,'Enable','on');
            case 'Pan'
                % enable pan feature
                pan(fGUI, 'on');
                % Continue with more cases as necessary. Might want to add a "grab'
                % feature to the image.
            otherwise
                % Code for when there is no match.
        end
        
    end

 end
