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

% Last Modified by GUIDE v2.5 01-Jul-2013 17:33:10

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


%%% KEY I/O FUNCTIONS
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
% this if statement only occurs if you call labelbox w/ no args.
if isempty(varargin); 
    varargin{1} = rand(100,100,5);
    varargin{2} = [1 1 1];
    warning('labelbox called with no args. dummy.')
end
handles.im = varargin{1};
handles.mid = ceil(size(handles.im,3)/2);
handles.top = size(handles.im,3);
handles.ax_pixdim = varargin{2};
handles.currsl = handles.mid;

% Choose default command line output for labelbox
handles.mask = false(size(handles.im));
handles.rgb = zeros([size(handles.im) 3]); 
handles.rgbvalid = false(size(handles.im,3),1); 

% plot stuff
handles = replotim(handles);

% fix slider increment
step = [1, 1] / (handles.top - 1);
set(handles.slider1,'Min',1,'Max',handles.top,'SliderStep',step,...
    'Value',handles.currsl);

% make done flag
handles.done = 1; 

% Update handles structure
guidata(hObject, handles);

% keeps figure open until imdone clicked
waitfor(handles.figure1,'Visible');
% UIWAIT makes labelbox wait for user response (see UIRESUME)
% uiwait(handles.figure1)

% --- Outputs from this function are returned to the command line.
function varargout = labelbox_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles); varargout{1} = 'empty'; disp('empty passed');
else varargout{1} = handles.mask; disp('mask passed');
end
closereq;

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
        set(handles.figure1,'Visible','off') % activates labelbox_OutputFcn
    else
        %nothing
    end
else
    set(handles.figure1,'Visible','off'); % activates labelbox_OutputFcn
end

% --- Executes on x button (close) press.
function xbutton_Callback(~, ~, ~) 
% hObject    handle to imdone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('xbutton pressed')
% g = questdlg(sprintf('Ready to close?\nEntered data will be lost'),...
%     'Close','YES','NO','YES');
% if strcmp(g,'YES')
%     closereq;
%     return; % 
% else
%     %nothing
% end
closereq;


%%% MOVEMENT FUNCTIONS %%%
% --- Executes on button press in sliceup.
function sliceup_Callback(hObject, ~, handles)
% keyboard;
% hObject    handle to sliceup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.currsl<handles.top
    % increment, redraw & update
    handles.currsl = handles.currsl+1;
    handles = replotim(handles); 
    changetxt(handles.mssgbox,'');
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
    handles.currsl = handles.currsl-1;
    handles = replotim(handles);
    changetxt(handles.mssgbox,'');
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
handles = replotim(handles);
changetxt(handles.mssgbox,'');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, ~)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%%% LABELLING FUNCTIONS %%%
% --- Executes on button press in labelme.
function labelme_Callback(hObject, ~, handles)
% hObject    handle to labelme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

changetxt(handles.mssgbox,'IN LABEL MODE');
    
% create impoly object
h = impoly(handles.axes1);
wait(h);
handles.mask(:,:, handles.currsl) = createMask(h); 
% on completion
changetxt(handles.mssgbox,'IN LABEL MODE\nMask Accepted!');
delete(h)

% draw ROI
handles.rgbvalid(handles.currsl) = false;
handles = replotim(handles);

% update handles (mask updated)
guidata(hObject,handles)

% --- Executes on button press in clearslice.
function clearslice_Callback(hObject, ~, handles)
% hObject    handle to clearslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask(:,:,handles.currsl) = zeros; 
changetxt(handles.mssgbox,'CLEARED MASK!');
handles = replotim(handles);
guidata(hObject, handles)

%%% MATH HEAVY FUNCTIONS %%%
% --- Executes on button press in logisticizeme.
function logisticizeme_Callback(~, ~, handles)
% hObject    handle to logisticizeme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('logisticizeme_Callback pressed');
keyboard;


%%% UTILITY/ FUNCTIONS (TOOLS) %%%
function changetxt(handle, str)
% CHANGETXT(handle,str)
set(handle,'String',sprintf(str));

function handles = replotim(handles)
% changeable constants 
alpha = .25; % change if you want to

% plot
axes(handles.axes1); 
tmp = handles.mask(:,:,handles.currsl);
% keyboard;
if ~any(tmp(:)); % no mask
    imshow(handles.im(:,:,handles.currsl));
elseif handles.rgbvalid(handles.currsl) % already done this mask
    imagesc(squeeze(handles.rgb(:,:,handles.currsl,:))); 
else % new mask 
    tmp = plot_segmentation_overlay(handles.im,handles.mask, alpha,cool,[0 1], ...
        handles.currsl);
    handles.rgb(:,:,handles.currsl,:) = permute(tmp,[1 2 4 3]); 
    handles.rgbvalid(handles.currsl)=true;
end
daspect(handles.ax_pixdim); axis off;

% change slice no string.
set(handles.slicefield,'String', ...
    sprintf('%3i / %3i',handles.currsl,handles.top));


% --- Executes on button press in render3d.
function render3d_Callback(~, ~, handles)
% hObject    handle to render3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('Generating 3D rendering. Please wait.\n');
g = figure; 

opts.cr_buffer = -1; 
opts.resdims = handles.ax_pixdim.^-1;
opts.xslices = round(size(handles.im,1)/2); 
opts.yslices = round(size(handles.im,2)/2); 
opts.zslices = round(size(handles.im,3)/2); 
opts.fignum = g; 
opts.slicealpha = .3;

render_3D_labels(handles.im, handles.mask,opts); 

