function varargout = vuSigmoidFunctionGUI(varargin)
% VUSIGMOIDFUNCTIONGUI M-file for vuSigmoidFunctionGUI.fig
%      VUSIGMOIDFUNCTIONGUI, by itself, creates a new VUSIGMOIDFUNCTIONGUI or raises the existing
%      singleton*.
%
%      H = VUSIGMOIDFUNCTIONGUI returns the handle to a new VUSIGMOIDFUNCTIONGUI or the handle to
%      the existing singleton*.
%
%      VUSIGMOIDFUNCTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUSIGMOIDFUNCTIONGUI.M with the given input arguments.
%
%      VUSIGMOIDFUNCTIONGUI('Property','Value',...) creates a new VUSIGMOIDFUNCTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuSigmoidFunctionGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuSigmoidFunctionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuSigmoidFunctionGUI

% Last Modified by GUIDE v2.5 10-Jul-2007 16:48:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuSigmoidFunctionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @vuSigmoidFunctionGUI_OutputFcn, ...
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


% --- Executes just before vuSigmoidFunctionGUI is made visible.
function vuSigmoidFunctionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuSigmoidFunctionGUI (see VARARGIN)

% Choose default command line output for vuSigmoidFunctionGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if nargin < 1
  error('MATLAB:vuSigmoidFunctionGUI:NotEnoughInputs', 'Not enough input arguments.');
end

% Get Figure initialized
X = squeeze(varargin{1});

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuSigmoidFunctionGUI:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 2)
  error('MATLAB:vuSigmoidFunctionGUI:UnknownDims', 'vuSigmoidFunctionGUI can only handle images of 2 dimensions.');
end

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

setappdata(gcf,'image',X);
setappdata(gcf,'sig_image',X);

minX = min(X.Data(:));
maxX = max(X.Data(:));
if (minX == maxX)
    error('MATLAB:vuSigmoidFunctionGUI:NoData','The image is blank!');
end

% Starting position
axes(handles.mainFig);
imagesc(X.Data)
axis image
colorbar
set(handles.alphaSlide,'Min',-1*maxX/2);
set(handles.alphaSlide,'Max',maxX/2);
set(handles.alphaSlide,'SliderStep',[1./500 0.1]);
set(handles.alphaSlide,'Value',1);
set(handles.alphaText,'String','Alpha Value = 1');
set(handles.betaSlide,'Min',minX);
set(handles.betaSlide,'Max',maxX);
set(handles.betaSlide,'SliderStep',[1./(maxX-minX+1) 0.1]);
set(handles.betaSlide,'Value',round((maxX-minX)/2+minX));
set(handles.betaText,'String',['Beta Value = ' num2str(round((maxX-minX)/2+minX))]);

sig_fx = (1).*(1./(1+exp(-1.*((minX:(maxX-minX)/1000:maxX)-get(handles.betaSlide,'Value'))./get(handles.alphaSlide,'Value')))) + 0;
axes(handles.sigFig)
plot(minX:(maxX-minX)/1000:maxX,sig_fx),axis tight

% UIWAIT makes vuSigmoidFunctionGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vuSigmoidFunctionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function betaSlide_Callback(hObject, eventdata, handles)
% hObject    handle to betaSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minX = get(handles.betaSlide,'Min');
maxX = get(handles.betaSlide,'Max');

set(handles.betaText,'String',['Beta Value = ' num2str(round(get(hObject,'Value')*10)/10)]);

sig_fx = (1).*(1./(1+exp(-1.*((minX:(maxX-minX)/1000:maxX)-get(handles.betaSlide,'Value'))./get(handles.alphaSlide,'Value')))) + 0;
axes(handles.sigFig)
plot(minX:(maxX-minX)/1000:maxX,sig_fx),axis tight

% --- Executes during object creation, after setting all properties.
function betaSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to betaSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function alphaSlide_Callback(hObject, eventdata, handles)
% hObject    handle to alphaSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minX = get(handles.betaSlide,'Min');
maxX = get(handles.betaSlide,'Max');

set(handles.alphaText,'String',['Alpha Value = ' num2str(round(get(hObject,'Value')*10)/10)]);

sig_fx = (1).*(1./(1+exp(-1.*((minX:(maxX-minX)/1000:maxX)-get(handles.betaSlide,'Value'))./get(handles.alphaSlide,'Value')))) + 0;
axes(handles.sigFig)
plot(minX:(maxX-minX)/1000:maxX,sig_fx),axis tight

% --- Executes during object creation, after setting all properties.
function alphaSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in update.
function update_Callback(hObject, eventdata, handles)
% hObject    handle to update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sig_im = vuSigmoidFunction(getappdata(gcf,'image'),'alpha',get(handles.alphaSlide,'Value'),'beta',get(handles.betaSlide,'Value'));
setappdata(gcf,'sig_image',sig_im);
axes(handles.mainFig);
imagesc(sig_im.Data)
axis image
colorbar

% --- Executes on button press in resetImage.
function resetImage_Callback(hObject, eventdata, handles)
% hObject    handle to resetImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sig_im = getappdata(gcf,'image');
setappdata(gcf,'sig_image',sig_im);
axes(handles.mainFig);
imagesc(sig_im.Data)
axis image
colorbar