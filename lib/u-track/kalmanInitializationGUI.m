function varargout = kalmanInitializationGUI(varargin)
% KALMANINITIALIZATIONGUI M-file for kalmanInitializationGUI.fig
%      KALMANINITIALIZATIONGUI, by itself, creates a new KALMANINITIALIZATIONGUI or raises the existing
%      singleton*.
%
%      H = KALMANINITIALIZATIONGUI returns the handle to a new KALMANINITIALIZATIONGUI or the handle to
%      the existing singleton*.
%
%      KALMANINITIALIZATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KALMANINITIALIZATIONGUI.M with the given input arguments.
%
%      KALMANINITIALIZATIONGUI('Property','Value',...) creates a new KALMANINITIALIZATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kalmanInitializationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kalmanInitializationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kalmanInitializationGUI

% Last Modified by GUIDE v2.5 10-Dec-2010 14:22:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kalmanInitializationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @kalmanInitializationGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before kalmanInitializationGUI is made visible.
function kalmanInitializationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.gapclosingFig = kalmanInitializationGUI{procID}('mainFig', handles.figure1, procID);
%
% userData.mainFig
% userData.procID
% userData.handles_main
% userData.userData_main
% userData.crtProc
% userData.parameters

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

handles.output = hObject;
set(handles.uipanel_1, 'SelectionChangeFcn', @uipanel_1_SelectionChangeFcn);
userData = get(handles.figure1, 'UserData');

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);
userData.userData_main = get(userData.handles_main.figure1, 'UserData');
userData.crtProc = userData.userData_main.crtProc;

u = get(userData.handles_main.popupmenu_kalman_initialize, 'UserData');
userData.kalmanInitParam = u{userData.procID};
kalmanInitParam = userData.kalmanInitParam;

% Parameter Setup
if ~isempty(kalmanInitParam)
   
    
    if ~isempty(kalmanInitParam.initVelocity)% Initial Valocity Estimate
       
        set(handles.edit_v_1, 'String', num2str(kalmanInitParam.initVelocity(1)))
        set(handles.edit_v_2, 'String', num2str(kalmanInitParam.initVelocity(2)))
        set(handles.edit_v_3, 'String', num2str(kalmanInitParam.initVelocity(3)))
        
    elseif ~isempty(kalmanInitParam.convergePoint) % Reference Point for Initial Estimate
        
        set(handles.edit_1, 'String', num2str(kalmanInitParam.convergePoint(1)))
        set(handles.edit_2, 'String', num2str(kalmanInitParam.convergePoint(2)))
        set(handles.edit_3, 'String', num2str(kalmanInitParam.convergePoint(3)))   
        
        set(handles.radiobutton_2, 'Value', 1)
        
        arrayfun(@(x)eval(['set(handles.text_v_',num2str(x),', ''Enable'', ''off'')']), 1:3)
        arrayfun(@(x)eval(['set(handles.edit_v_',num2str(x),', ''Enable'', ''off'')']), 1:3)
        
        arrayfun(@(x)eval(['set(handles.text_',num2str(x),', ''Enable'', ''on'')']), 1:3)
        arrayfun(@(x)eval(['set(handles.edit_',num2str(x),', ''Enable'', ''on'')']), 1:3)          
        
    else
        error('User-defined: kalman initialization parameter error.')
       
    end
    
    set(handles.edit_radius, 'String', num2str(kalmanInitParam.searchRadiusFirstIteration))
end

% Get icon infomation
userData.questIconData = userData.userData_main.questIconData;
userData.colormap = userData.userData_main.colormap;

% ----------------------Set up help icon------------------------

% Set up help icon
set(hObject,'colormap',userData.colormap);
% Set up package help. Package icon is tagged as '0'
axes(handles.axes_help);
Img = image(userData.questIconData); 
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
if openHelpFile
    set(Img, 'UserData', struct('class', 'kalmanInitializationGUI'))
else
    set(Img, 'UserData', 'Please refer to help file.')
end



set(handles.figure1, 'UserData', userData)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kalmanInitializationGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = kalmanInitializationGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~get(handles.radiobutton_1, 'Value') && ~get(handles.radiobutton_2, 'Value')
    
    delete(handles.figure1);
    return
end

userData = get(handles.figure1, 'UserData');
kalmanInitParam = userData.kalmanInitParam;

v{1} = get(handles.edit_v_1, 'String');
v{2} = get(handles.edit_v_2, 'String');
v{3} = get(handles.edit_v_3, 'String');

c{1} = get(handles.edit_1, 'String');
c{2} = get(handles.edit_2, 'String');
c{3} = get(handles.edit_3, 'String');

searchRadiusFirstIteration = get(handles.edit_radius, 'String');

e(1) = all(cellfun(@(x)isempty(x), v));
e(2) = all(cellfun(@(x)isempty(x), v));
e(3) = isempty(searchRadiusFirstIteration);

if ~all(e)
   
    if get(handles.radiobutton_1, 'Value')
        
        % Initial Velocity Estimate
        temp = cellfun(@(x)isempty(x), v);
        if any(temp)
            errordlg('Please provide all three values to "vX", "vY" and "vZ" respectively.','Error','modal')
            return

        elseif any(cellfun(@(x)isnan(str2double(x)), v)) || any(cellfun(@(x)(str2double(x) < 0), v))
            errordlg('Please provide a valid value to parameter "Initial Velocity Estimate".','Error','modal')
            return

        else
            v = cellfun(@(x)str2double(x), v, 'UniformOutput', true);
        end

    else
        
        % Reference Point for Initial Velocity Estimate
        temp = cellfun(@(x)isempty(x), c);
        if any(temp)
            errordlg('Please provide all three values to "X", "Y" and "Z" respectively.','Error','modal')
            return

        elseif any(cellfun(@(x)isnan(str2double(x)), c)) || any(cellfun(@(x)(str2double(x) < 0), c))
            errordlg('Please provide a valid value to parameter "Initial Velocity Estimate".','Error','modal')
            return

        else
            c = cellfun(@(x)str2double(x), c, 'UniformOutput', true);
        end           
    end
    
    if isempty( searchRadiusFirstIteration )
        errordlg('Parameter "Search Radius for Iteration" is requied by the algorithm.','Error','modal')
        return

    elseif isnan(str2double(searchRadiusFirstIteration)) || str2double(searchRadiusFirstIteration) < 0
        errordlg('Please provide a valid value to parameter "Search Radius for Iteration".','Error','modal')
        return

    else
        searchRadiusFirstIteration = str2double(searchRadiusFirstIteration);
    end    
    
    % Set Parameters
    
    kalmanInitParam.searchRadiusFirstIteration = searchRadiusFirstIteration;
    
    if get(handles.radiobutton_1, 'Value')
        
        kalmanInitParam.initVelocity = v;
        kalmanInitParam.convergePoint = [];
    else
        kalmanInitParam.initVelocity = [];
        kalmanInitParam.convergePoint = c;        
    end

else
    
    kalmanInitParam = [];
end

u = get(userData.handles_main.popupmenu_kalman_initialize, 'UserData');
u{userData.procID} = kalmanInitParam;

set(userData.handles_main.popupmenu_kalman_initialize, 'UserData', u)   

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);

% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

arrayfun(@(x)eval(['set(handles.edit_',num2str(x),', ''String'', [])']), 1:3)
arrayfun(@(x)eval(['set(handles.edit_v_',num2str(x),', ''String'', [])']), 1:3)
set(handles.edit_radius, 'String', [])

function edit_radius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_radius as text
%        str2double(get(hObject,'String')) returns contents of edit_radius as a double


% --- Executes during object creation, after setting all properties.
function edit_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_v_3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_v_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_v_3 as text
%        str2double(get(hObject,'String')) returns contents of edit_v_3 as a double


% --- Executes during object creation, after setting all properties.
function edit_v_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_v_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_v_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_v_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_v_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_v_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_v_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_v_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_v_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_v_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_v_1 as text
%        str2double(get(hObject,'String')) returns contents of edit_v_1 as a double


% --- Executes during object creation, after setting all properties.
function edit_v_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_v_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_3 as text
%        str2double(get(hObject,'String')) returns contents of edit_3 as a double


% --- Executes during object creation, after setting all properties.
function edit_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_1 as text
%        str2double(get(hObject,'String')) returns contents of edit_1 as a double


% --- Executes during object creation, after setting all properties.
function edit_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uipanel_1_SelectionChangeFcn(hObject, eventdata)
% Call back function of ration button group uipanel_1
handles = guidata(hObject); 

% Highlight the content under new radiobutton
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton_1'
        
        arrayfun(@(x)eval(['set(handles.text_v_',num2str(x),', ''Enable'', ''on'')']), 1:3)
        arrayfun(@(x)eval(['set(handles.edit_v_',num2str(x),', ''Enable'', ''on'')']), 1:3)
        
        arrayfun(@(x)eval(['set(handles.text_',num2str(x),', ''Enable'', ''off'')']), 1:3)
        arrayfun(@(x)eval(['set(handles.edit_',num2str(x),', ''Enable'', ''off'')']), 1:3)
        
    case 'radiobutton_2'

        arrayfun(@(x)eval(['set(handles.text_v_',num2str(x),', ''Enable'', ''off'')']), 1:3)
        arrayfun(@(x)eval(['set(handles.edit_v_',num2str(x),', ''Enable'', ''off'')']), 1:3)
        
        arrayfun(@(x)eval(['set(handles.text_',num2str(x),', ''Enable'', ''on'')']), 1:3)
        arrayfun(@(x)eval(['set(handles.edit_',num2str(x),', ''Enable'', ''on'')']), 1:3)        
        
    otherwise
       disp('User-defined Warning: No radio button tag is ',...
           'found when SelectionChangeFcn is triggered.');
       return;
end

%updates the handles structure
guidata(hObject, handles);
