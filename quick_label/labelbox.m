function varargout = labelbox(varargin)
% LABELBOX MATLAB code for labelbox.fig
%      LABELBOX, by itself, creates a new LABELBOX or raises the existing
%      singleton*.
%
%      H = LABELBOX returns the handle to a new LABELBOX or the handle to
%      the existing singleton*.
%
%      LABELBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LABELBOX.M with the given input arguments.
%
%      LABELBOX('Property','Value',...) creates a new LABELBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before labelbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to labelbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help labelbox

% Last Modified by GUIDE v2.5 27-Jun-2013 16:48:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @labelbox_OpeningFcn, ...
                   'gui_OutputFcn',  @labelbox_OutputFcn, ...
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

% --- Executes just before labelbox is made visible.
function labelbox_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to labelbox (see VARARGIN)

% change what clicking "x" does.
% it now calls the imdone_Callback function
set(handles.figure1,'CloseRequestFcn',@(hObject,eventdata)labelbox(...
    'xbutton_Callback',hObject,eventdata,guidata(hObject)))

% change window name from "labelbox" 
set(handles.figure1,'Name','Awesome box!');

% get data from command line and put in handles structure
if isempty(varargin); 
    varargin{1} = padarray(zeros([64 64]),[0 0 1],1);
    varargin{2} = [1 1 1];
end
    
handles.im = varargin{1};
handles.mid = ceil(size(handles.im,3)/2);
handles.top = size(handles.im,3);
handles.ax_pixdim = varargin{2};
handles.currsl = handles.mid;

% Choose default command line output for labelbox
handles.mask = false(size(handles.im));

% Update handles structure
guidata(hObject, handles);

% plot stuff
axes(handles.axes1); 
    imshow(handles.im(:,:,handles.currsl)); daspect(handles.ax_pixdim); axis off; 
set(handles.slicefield,'String', ...
    sprintf('%3i / %3i',handles.currsl,handles.top));

% fix slider increment
step = [1, 1] / (handles.top - 1);
set(handles.slider1,'Min',1,'Max',handles.top,'SliderStep',step,...
    'Value',handles.currsl);

% UIWAIT makes labelbox wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = labelbox_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles); varargout{1} = 'empty';
else varargout{1} = handles.mask;
end
closereq;

% --- Executes on button press in sliceup.
function sliceup_Callback(hObject, ~, handles)
% hObject    handle to sliceup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.currsl<handles.top
    % redraw & increment & update
    axes(handles.axes1); 
        imshow(handles.im(:,:,handles.currsl)); 
        daspect(handles.ax_pixdim); axis off;
    handles.currsl = handles.currsl+1; 
    set(handles.slicefield,'String', ...
        sprintf('%3i / %3i',handles.currsl,handles.top));
    set(handles.slider1,'Value',handles.currsl);
    guidata(hObject, handles);
else
    fprintf('\nTOP OF IMAGE. NO MORE UPPY FOR YOU!');
end

% --- Executes on button press in slicedown.
function slicedown_Callback(hObject, ~, handles)
% hObject    handle to slicedown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.currsl>1
    % redraw & increment & update
    axes(handles.axes1); 
        imshow(handles.im(:,:,handles.currsl)); 
        daspect(handles.ax_pixdim); axis off;
    handles.currsl = handles.currsl-1; 
    set(handles.slicefield,'String', ...
        sprintf('%3i / %3i',handles.currsl,handles.top));
    set(handles.slider1,'Value',handles.currsl);
    guidata(hObject, handles);
else
    fprintf('\nBOTTOM OF IMAGE. NO MORE DOWNY FOR YOU!');
end

% --- Executes on slider movement.
function slider1_Callback(hObject, ~, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.currsl = round(get(handles.slider1,'Value'));
axes(handles.axes1); 
    imshow(handles.im(:,:,handles.currsl)); 
    daspect(handles.ax_pixdim); axis off;
set(handles.slicefield,'String', ...
    sprintf('%3i / %3i',handles.currsl,handles.top));
guidata(hObject, handles);

% --- Executes on button press in clearslice.
function clearslice_Callback(~, ~, ~)
% hObject    handle to clearslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('pushbutton3 pressed');

% --- Executes on button press in imdone.
function imdone_Callback(~, ~, handles) %#ok<*DEFNU>
% hObject    handle to imdone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('imdone pressed')
if ~isequal(zeros(size(handles.mask)),handles.mask)
    g = questdlg('You have entered data on the mask. Ready to close?',...
        'Close','YES','NO','NO');
    if strcmp(g,'YES')
        uiresume; % activates labelbox_OutputFcn
    else
        %nothing
    end
else
    uiresume; % activates labelbox_OutputFcn
end

% --- Executes on x button (close) press.
function xbutton_Callback(~, ~, ~) 
% hObject    handle to imdone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('xbutton pressed')
g = questdlg(sprintf('Ready to close?\nEntered data will be lost'),...
    'Close','YES','NO','NO');
if strcmp(g,'YES')
    closereq;
    return; % 
else
    %nothing
end

% --- Executes on button press in logisticizeme.
function logisticizeme_Callback(~, ~, ~)
% hObject    handle to logisticizeme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('logisticizeme_Callback pressed');


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, ~)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
