function h1 = labelbox_export()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
% 
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.


load labelbox_export.mat


appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', [], ...
    'taginfo', struct(...
    'figure', 2, ...
    'pushbutton', 19, ...
    'text', 17, ...
    'edit', 6, ...
    'frame', 4, ...
    'radiobutton', 19, ...
    'uipanel', 11, ...
    'axes', 2, ...
    'uitoolbar', 2, ...
    'uitoggletool', 7, ...
    'uipushtool', 2), ...
    'override', 1, ...
    'release', 13, ...
    'resize', 'simple', ...
    'accessibility', 'callback', ...
    'mfile', 0, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastFilename', '/home-local/bryanfw/fwbsrc/quick_label/labelbox.fig');
appdata.lastValidTag = 'figure1';
appdata.GUIDELayoutEditor = [];
appdata.initTags = struct(...
    'handle', [], ...
    'tag', 'figure1');

h1 = figure(...
'Color',[0.701960784313725 0.701960784313725 0.701960784313725],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'DockControls','off',...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','labelbox',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[544 55 742 570],...
'HandleVisibility','callback',...
'UserData',[],...
'Tag','figure1',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'calculate';

h2 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Callback',@(hObject,eventdata)labelbox('calculate_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'ListboxTop',0,...
'Position',[0.788409703504043 0.1 0.172506738544474 0.0526315789473684],...
'String','Logisticize Me!',...
'UserData',[],...
'Tag','calculate',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'axes1';

h3 = axes(...
'Parent',h1,...
'Position',[0.0202156334231806 0.0333333333333333 0.749326145552561 0.93859649122807],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'LooseInset',[0.136488201867296 0.121219953001468 0.0997413782876395 0.0826499679555462],...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','axes1',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h4 = get(h3,'title');

set(h4,...
'Parent',h3,...
'Units','data',...
'FontUnits','points',...
'BackgroundColor','none',...
'Color',[0 0 0],...
'DisplayName',blanks(0),...
'EdgeColor','none',...
'EraseMode','normal',...
'DVIMode','auto',...
'FontAngle','normal',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[0.49910071942446 1.01401869158878 1.00005459937205],...
'Rotation',0,...
'String',blanks(0),...
'Interpreter','tex',...
'VerticalAlignment','bottom',...
'ButtonDownFcn',[],...
'CreateFcn', {@local_CreateFcn, [], ''} ,...
'DeleteFcn',[],...
'BusyAction','queue',...
'HandleVisibility','off',...
'HelpTopicKey',blanks(0),...
'HitTest','on',...
'Interruptible','on',...
'SelectionHighlight','on',...
'Serializable','on',...
'Tag',blanks(0),...
'UserData',[],...
'Visible','on',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'IncludeRenderer','on',...
'Clipping','off');

h5 = get(h3,'xlabel');

set(h5,...
'Parent',h3,...
'Units','data',...
'FontUnits','points',...
'BackgroundColor','none',...
'Color',[0 0 0],...
'DisplayName',blanks(0),...
'EdgeColor','none',...
'EraseMode','normal',...
'DVIMode','auto',...
'FontAngle','normal',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[0.49910071942446 -0.0476635514018693 1.00005459937205],...
'Rotation',0,...
'String',blanks(0),...
'Interpreter','tex',...
'VerticalAlignment','cap',...
'ButtonDownFcn',[],...
'CreateFcn', {@local_CreateFcn, [], ''} ,...
'DeleteFcn',[],...
'BusyAction','queue',...
'HandleVisibility','off',...
'HelpTopicKey',blanks(0),...
'HitTest','on',...
'Interruptible','on',...
'SelectionHighlight','on',...
'Serializable','on',...
'Tag',blanks(0),...
'UserData',[],...
'Visible','on',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'IncludeRenderer','on',...
'Clipping','off');

h6 = get(h3,'ylabel');

set(h6,...
'Parent',h3,...
'Units','data',...
'FontUnits','points',...
'BackgroundColor','none',...
'Color',[0 0 0],...
'DisplayName',blanks(0),...
'EdgeColor','none',...
'EraseMode','normal',...
'DVIMode','auto',...
'FontAngle','normal',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[-0.0494604316546763 0.498130841121495 1.00005459937205],...
'Rotation',90,...
'String',blanks(0),...
'Interpreter','tex',...
'VerticalAlignment','bottom',...
'ButtonDownFcn',[],...
'CreateFcn', {@local_CreateFcn, [], ''} ,...
'DeleteFcn',[],...
'BusyAction','queue',...
'HandleVisibility','off',...
'HelpTopicKey',blanks(0),...
'HitTest','on',...
'Interruptible','on',...
'SelectionHighlight','on',...
'Serializable','on',...
'Tag',blanks(0),...
'UserData',[],...
'Visible','on',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'IncludeRenderer','on',...
'Clipping','off');

h7 = get(h3,'zlabel');

set(h7,...
'Parent',h3,...
'Units','data',...
'FontUnits','points',...
'BackgroundColor','none',...
'Color',[0 0 0],...
'DisplayName',blanks(0),...
'EdgeColor','none',...
'EraseMode','normal',...
'DVIMode','auto',...
'FontAngle','normal',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','right',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[-0.0278776978417266 1.02710280373832 1.00005459937205],...
'Rotation',0,...
'String',blanks(0),...
'Interpreter','tex',...
'VerticalAlignment','middle',...
'ButtonDownFcn',[],...
'CreateFcn', {@local_CreateFcn, [], ''} ,...
'DeleteFcn',[],...
'BusyAction','queue',...
'HandleVisibility','off',...
'HelpTopicKey',blanks(0),...
'HitTest','on',...
'Interruptible','on',...
'SelectionHighlight','on',...
'Serializable','on',...
'Tag',blanks(0),...
'UserData',[],...
'Visible','off',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'IncludeRenderer','on',...
'Clipping','off');

appdata = [];
appdata.lastValidTag = 'pushbutton9';

h8 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Callback',@(hObject,eventdata)labelbox('pushbutton9_Callback',hObject,eventdata,guidata(hObject)),...
'Position',[0.784366576819407 0.733333333333333 0.0444743935309973 0.236842105263158],...
'String','+',...
'Tag','pushbutton9',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'Untitled_1';

h9 = uimenu(...
'Parent',h1,...
'Callback',@(hObject,eventdata)labelbox('Untitled_1_Callback',hObject,eventdata,guidata(hObject)),...
'Label','Untitled 1',...
'Tag','Untitled_1',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uitoolbar1';

h10 = uitoolbar(...
'Parent',h1,...
'Tag','uitoolbar1',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Exploration.Pan';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'uitoggletool5';

h11 = uitoggletool(...
'Parent',h10,...
'ClickedCallback','%default',...
'CData',mat{1},...
'TooltipString','Pan',...
'Tag','uitoggletool5',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Exploration.ZoomOut';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'uitoggletool2';

h12 = uitoggletool(...
'Parent',h10,...
'ClickedCallback','%default',...
'CData',mat{2},...
'TooltipString','Zoom Out',...
'Tag','uitoggletool2',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Exploration.ZoomIn';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'uitoggletool3';

h13 = uitoggletool(...
'Parent',h10,...
'ClickedCallback','%default',...
'CData',mat{3},...
'TooltipString','Zoom In',...
'Tag','uitoggletool3',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Exploration.DataCursor';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'uitoggletool1';

h14 = uitoggletool(...
'Parent',h10,...
'ClickedCallback','%default',...
'CData',mat{4},...
'TooltipString','Data Cursor',...
'Tag','uitoggletool1',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Annotation.InsertColorbar';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'uitoggletool4';

h15 = uitoggletool(...
'Parent',h10,...
'ClickedCallback','%default',...
'CData',mat{5},...
'TooltipString','Insert Colorbar',...
'Tag','uitoggletool4',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton13';

h16 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Callback',@(hObject,eventdata)labelbox('pushbutton13_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'ListboxTop',0,...
'Position',[0.788409703504043 0.166666666666667 0.172506738544474 0.0526315789473684],...
'String','Clear Slice',...
'UserData',[],...
'Tag','pushbutton13',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton14';

h17 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Callback',@(hObject,eventdata)labelbox('pushbutton14_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'ListboxTop',0,...
'Position',[0.788409703504043 0.0333333333333333 0.172506738544474 0.0526315789473684],...
'String','I''m done.',...
'UserData',[],...
'Tag','pushbutton14',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton16';

h18 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Callback',@(hObject,eventdata)labelbox('pushbutton16_Callback',hObject,eventdata,guidata(hObject)),...
'Position',[0.78167115902965 0.471929824561404 0.0444743935309973 0.236842105263158],...
'String','-',...
'Tag','pushbutton16',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text15';

h19 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Position',[0.862533692722372 0.908771929824561 0.101078167115903 0.0368421052631579],...
'String','Slice:',...
'Style','text',...
'Tag','text15',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text16';

h20 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'FontWeight','bold',...
'Position',[0.86522911051213 0.863157894736842 0.0956873315363881 0.0473684210526316],...
'String','NOTESASD',...
'Style','text',...
'Tag','text16',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );



% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   if isa(createfcn,'function_handle')
       createfcn(hObject, eventdata);
   else
       eval(createfcn);
   end
end