function varargout = vuWatershedSegmentationGUI(varargin)
% VUWATERSHEDSEGMENTATIONGUI M-file for vuWatershedSegmentationGUI.fig
%      VUWATERSHEDSEGMENTATIONGUI, by itself, creates a new VUWATERSHEDSEGMENTATIONGUI or raises the existing
%      singleton*.
%
%      H = VUWATERSHEDSEGMENTATIONGUI returns the handle to a new VUWATERSHEDSEGMENTATIONGUI or the handle to
%      the existing singleton*.
%
%      VUWATERSHEDSEGMENTATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUWATERSHEDSEGMENTATIONGUI.M with the given input arguments.
%
%      VUWATERSHEDSEGMENTATIONGUI('Property','Value',...) creates a new VUWATERSHEDSEGMENTATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuWatershedSegmentationGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuWatershedSegmentationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuWatershedSegmentationGUI

% Last Modified by GUIDE v2.5 11-Jul-2007 09:06:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuWatershedSegmentationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @vuWatershedSegmentationGUI_OutputFcn, ...
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


% --- Executes just before vuWatershedSegmentationGUI is made visible.
function vuWatershedSegmentationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuWatershedSegmentationGUI (see VARARGIN)

% Choose default command line output for vuWatershedSegmentationGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if nargin < 1
  error('MATLAB:vuWatershedSegmentationGUI:NotEnoughInputs', 'Not enough input arguments.');
end

% Get Figure initialized
X = squeeze(varargin{1});

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuWatershedSegmentationGUI:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 2)
  error('MATLAB:vuWatershedSegmentationGUI:UnknownDims', 'vuWatershedSegmentationGUI can only handle images of 2 dimensions.');
end

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

setappdata(gcf,'image',X);
setappdata(gcf,'seg_image',X);

% Starting position
set(handles.level_slide,'Value',0.001);
set(handles.threshold_slide,'Value',0.15);
set(handles.rowSlide,'Min',1);
set(handles.rowSlide,'Max',X.Dims(1));
set(handles.rowSlide,'SliderStep',[1./X.Dims(1) 0.1]);
set(handles.rowSlide,'Value',round((X.Dims(1))/2));
rowNumber = get(handles.rowSlide,'Value');
axes(handles.mainFig);
imagesc(X.Data)
h = line([1 X.Dims(2)],[X.Dims(1)-rowNumber X.Dims(1)-rowNumber],'Color','r');
axis image
setappdata(gcf,'line_h',h);

rowData = X.Data(X.Dims(1)-rowNumber,:);
maxX = max(X.Data(:));
minX = min(X.Data(:));
setappdata(gcf,'maxX',maxX);
setappdata(gcf,'minX',minX);
rangeX = maxX - minX;
rowData(rowData<(0.15*(rangeX) + minX)) = 0.15*(rangeX) + minX;
newMin = 0.15*(rangeX) + minX;
level = 0.001*(maxX - newMin) + newMin;
axes(handles.watershedFig)
ylim([minX maxX]);
plot(rowData,'Color','r');
line([1 X.Dims(2)],[level level],'Color','b');
% UIWAIT makes vuWatershedSegmentationGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vuWatershedSegmentationGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
seg_im = vuWatershedSegmentation(getappdata(gcf,'image'),'level',get(handles.level_slide,'Value'),'threshold',get(handles.threshold_slide,'Value'));
h = getappdata(gcf,'line_h');
delete(h);
setappdata(gcf,'seg_image',seg_im);
rowNumber = get(handles.rowSlide,'Value');
axes(handles.mainFig);
imagesc(uint8(seg_im.Data))
h = line([1 seg_im.Dims(2)],[seg_im.Dims(1)-rowNumber seg_im.Dims(1)-rowNumber],'Color','r');
axis image
setappdata(gcf,'line_h',h);
assignin('base','WatershedSementationGUI_Output',seg_im);

% --- Executes on slider movement.
function level_slide_Callback(hObject, eventdata, handles)
% hObject    handle to level_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(hObject,'Value',round(get(hObject,'Value')*1000)/1000)
set(handles.LevelText,'String',['Watershed Level = ' num2str(get(hObject,'Value')*100) '%']);

% Watershed Plot
X = getappdata(gcf,'image');
rowNumber = round(get(handles.rowSlide,'Value'));
maxX = getappdata(gcf,'maxX');
minX = getappdata(gcf,'minX');
level = get(handles.level_slide,'Value');
threshold = get(handles.threshold_slide,'Value');
rowData = X.Data(X.Dims(1)-rowNumber,:);
rangeX = maxX - minX;
rowData(rowData<(threshold*(rangeX) + minX)) = threshold*(rangeX) + minX;
newMin = threshold*(rangeX) + minX;
level = level*(maxX - newMin) + newMin;
axes(handles.watershedFig)
plot(rowData,'Color','r');
line([1 X.Dims(2)],[level level],'Color','b');
ylim([minX maxX]);

% --- Executes during object creation, after setting all properties.
function level_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to level_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function threshold_slide_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(hObject,'Value',round(get(hObject,'Value')*1000)/1000)
set(handles.ThresholdText,'String',['Watershed Threshold = ' num2str(get(hObject,'Value')*100) '%']);

% Watershed Plot
X = getappdata(gcf,'image');
rowNumber = round(get(handles.rowSlide,'Value'));
maxX = getappdata(gcf,'maxX');
minX = getappdata(gcf,'minX');
rangeX = maxX - minX;
level = get(handles.level_slide,'Value');
threshold = get(handles.threshold_slide,'Value');
rowData = X.Data(X.Dims(1)-rowNumber,:);
rowData(rowData<(threshold*(rangeX) + minX)) = threshold*(rangeX) + minX;
newMin = threshold*(rangeX) + minX;
level = level*(maxX - newMin) + newMin;
axes(handles.watershedFig)
plot(rowData,'Color','r');
line([1 X.Dims(2)],[level level],'Color','b');
ylim([minX maxX]);

% --- Executes during object creation, after setting all properties.
function threshold_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function rowSlide_Callback(hObject, eventdata, handles)
% hObject    handle to rowSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
X = getappdata(gcf,'image');
h = getappdata(gcf,'line_h');
rowNumber = round(get(handles.rowSlide,'Value'));
axes(handles.mainFig);
delete(h);
h = line([1 X.Dims(2)],[X.Dims(1)-rowNumber X.Dims(1)-rowNumber],'Color','r');
axis image
setappdata(gcf,'line_h',h);

% Watershed Plot
maxX = getappdata(gcf,'maxX');
minX = getappdata(gcf,'minX');
level = get(handles.level_slide,'Value');
threshold = get(handles.threshold_slide,'Value');
rowData = X.Data(X.Dims(1)-rowNumber,:);
rangeX = maxX - minX;
rowData(rowData<(threshold*(rangeX) + minX)) = threshold*(rangeX) + minX;
newMin = threshold*(rangeX) + minX;
level = level*(maxX - newMin) + newMin;
axes(handles.watershedFig)
plot(rowData,'Color','r');
line([1 X.Dims(2)],[level level],'Color','b');
ylim([minX maxX]);
% --- Executes during object creation, after setting all properties.
function rowSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rowSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in resetImage.
function resetImage_Callback(hObject, eventdata, handles)
% hObject    handle to resetImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = getappdata(gcf,'image');
h = getappdata(gcf,'line_h');
rowNumber = get(handles.rowSlide,'Value');
delete(h);
axes(handles.mainFig);
imagesc(X.Data)
h = line([1 X.Dims(2)],[X.Dims(1)-rowNumber X.Dims(1)-rowNumber],'Color','r');
axis image
setappdata(gcf,'line_h',h);

