function [dicom_filename,dicomdir_filepath] = getfileDICOMDIR(dicomdir_filename)

%% Initialize dicom_filename to empty
dicom_filename='';

%% Allow user to select a file if input FILENAME is not provided or is empty
if nargin<1,
    dicomdir_filename=pwd;
end
if length(dicomdir_filename)==0 | isdir(dicomdir_filename),
    [fn, pn] = uigetfile({'*.*','DICOMDIR files'},'Select a DICOMDIR file',dicomdir_filename);
    if fn~=0,
        dicomdir_filename = sprintf('%s%s',pn,fn);
    else
        %disp('getfileDICOMDIR cancelled');
        dicomdir_filepath = '';
        return;
    end
end

%% Assign dicomdir_filepath
[pathstr,name,ext,versn] = fileparts(dicomdir_filename);
dicomdir_filepath = pathstr;

%% Load dicominfo and create dicomdir structure
dicomdir = dicomdir_struct(dicomdir_filename);

%% iconpath parent to custom icons
custom_icons_subpath = 'svn\My Documents\matlab';
custom_icons_parentpath = sprintf('C:\\%s',custom_icons_subpath);
if ~isdir(custom_icons_parentpath),
    % this is within an MCR distributable
    mcr_path = pwd;
    [pathstr,mcr_name,ext,versn] = fileparts(mcr_path);
    custom_icons_parentpath = sprintf('%s\\%s_mcr\\%s', mcr_path, mcr_name, custom_icons_subpath);
end
%disp( sprintf('custom_icons_parentpath = %s', custom_icons_parentpath) );

%% Initialize figure
h = figure('Units', 'Pixels', 'Color', [1 1 1], 'Position', [50 50 1080 800]);
set(h,'NumberTitle','off','Name',sprintf('getfileDICOMDIR - %s ::: DOUBLE CLICK SERIES TO SELECT', dicomdir_filename) );
set(h,'MenuBar','none','ToolBar','none');
set(h,'Color', [1 1 1]);
root = uitreenode('root', dicomdir.patient_header_str.str, [custom_icons_parentpath, '\folderclosed.gif'], false);
mtree = uitree(h,'Root',root,'ExpandFcn', @myExpfcn);
set(mtree, 'NodeExpandedCallback', {@myNodeExpanded, mtree, @myExpfcn});
set(mtree, 'NodeWillExpandCallback', @myNodeWillExpand);
set(mtree, 'NodeSelectedCallback', @myNodeSelected);
set(mtree, 'Units', 'normalized', 'position', [0 0 1 1]);

jtree = mtree.Tree;
set(jtree,'ShowsRootHandles','on');
set(jtree,'RightSelectionEnabled','on');

%% Store figure handle with dicomdir struct
dicomdir.h = h;

%% Initialize dicomdir.dicom_filename to empty
dicomdir.dicom_filename = '';

%% Set dicominfo structure to UserData of tree.FigureComponent
set(mtree.FigureComponent,'UserData',dicomdir);

%% Set UITREE font to courier for fixed spacing
font_courier = javaObject('javax.swing.plaf.FontUIResource','Courier', 0, 11);
jtree.setFont(font_courier);

%%
set(jtree,'MouseClickedCallback',{@myMouseClicked, mtree});

%% Expand patient nodes by default
patient_nodes = autoexpand(mtree,root);

%% Expand first patient node by default
if length(patient_nodes)>=1,
    study_nodes = autoexpand(mtree,patient_nodes(1));   
end

%% Wait for user completion
% Make the GUI blocking.
uiwait(h);

%% Resume and return dicom_filename
try
    dicomdir = get(mtree.FigureComponent,'UserData');
    dicom_filename = dicomdir.dicom_filename;
    delete(h);
catch
    try
        delete(h);
    catch
        % nothing to do
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cNode = myNodeSelected(tree,evd)

tmp = tree.FigureComponent;
dicomdir = get(tmp,'UserData');
evdnode  = evd.getCurrentNode;
dicomdir.selected_node_value = evdnode.getValue;
set(tmp,'UserData',dicomdir);
%disp( sprintf('nodeSelected : %s', dicomdir.selected_node_value) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cNode = myNodeWillExpand(tree,ev)
cNode = ev.getCurrentNode;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodes = myExpfcn(tree, value)

tmp = tree.FigureComponent;
dicomdir = get(tmp,'UserData');

%% iconpath parent to custom icons
custom_icons_subpath = 'svn\My Documents\matlab';
custom_icons_parentpath = sprintf('C:\\%s',custom_icons_subpath);
if ~isdir(custom_icons_parentpath),
    % this is within an MCR distributable
    mcr_path = pwd;
    [pathstr,mcr_name,ext,versn] = fileparts(mcr_path);
    custom_icons_parentpath = sprintf('%s\\%s_mcr\\%s', mcr_path, mcr_name, custom_icons_subpath);
end
%disp( sprintf('custom_icons_parentpath = %s', custom_icons_parentpath) );

% try
        
      if strcmp(value,'root'), % root level expands to patients
                        
            childnode_iconpath = [custom_icons_parentpath, '/folderclosed.gif'];
            childnode_isleaf = false;
            
            patients = fieldnames(dicomdir.tree);
            
            count = 0;
            for k=1:length(patients),
                childnode_value = patients{k};
                childnode_description = dicomdir.tree.(patients{k}).patient_description_str;
                count=count+1;
                nodes(count) = uitreenode(childnode_value, childnode_description, childnode_iconpath, childnode_isleaf);
            end

      % levels are: PATIENT -> STUDY -> SERIES -> IMAGE      
            
      elseif findstr(value,'SERIES'), % series expands to image
              
            childnode_iconpath = [matlabroot, '/toolbox/matlab/icons/greenarrowicon.gif'];
            childnode_isleaf = true; % images are not expandable
                        
            matches = regexp(value,'\w*','match');
            patient_name = matches{1};
            study_name = matches{2};
            series_name = matches{3};
                   
            image_count=0;
            images={};
            
            matches = regexp(fieldnames( dicomdir.tree.(patient_name).(study_name).(series_name) ),'^IMAGE_.+$','match');
            for k=1:length(matches),
                if ~isempty(matches{k}),
                    image_count=image_count+1;
                    images{image_count} = char(matches{k});
                end
            end

            matches = regexp(fieldnames( dicomdir.tree.(patient_name).(study_name).(series_name) ),'^PRIVATE_.+$','match');
            for k=1:length(matches),
                if ~isempty(matches{k}),
                    image_count=image_count+1;
                    images{image_count} = char(matches{k});
                end
            end

            matches = regexp(fieldnames( dicomdir.tree.(patient_name).(study_name).(series_name) ),'^PRESENTATION_.+$','match');
            for k=1:length(matches),
                if ~isempty(matches{k}),
                    image_count=image_count+1;
                    images{image_count} = char(matches{k});
                end
            end             
            
            count = 0;
            for k=1:length(images),
                childnode_value = sprintf('%s %s %s %s',patient_name, study_name, series_name, images{k});
                childnode_description = dicomdir.tree.(patient_name).(study_name).(series_name).(images{k}).image_description_str;
                count=count+1;
                nodes(count) = uitreenode(childnode_value, childnode_description, childnode_iconpath, childnode_isleaf);
            end
            
      elseif findstr(value,'STUDY'), % study expands to series
              
            childnode_iconpath = [custom_icons_parentpath, '/mrseries.gif'];
            childnode_isleaf = true; % do not allow series to expand for getfileDICOMDIR
                        
            matches = regexp(value,'\w*','match');
            patient = matches{1};
            study = matches{2};
            
            matches = regexp(fieldnames( dicomdir.tree.(patient).(study) ),'^SERIES_.+$','match');
            series_count=0;
            series={};
            for k=1:length(matches),
                if ~isempty(matches{k}),
                    series_count=series_count+1;
                    series{series_count} = char(matches{k});
                end
            end
            
            count = 0;
            for k=1:length(series),
                if isfield(dicomdir.tree.(patient).(study).(series{k}),'SeriesDescription'),
                    series_number = dicomdir.tree.(patient).(study).(series{k}).SeriesNumber;
                    if series_number>0,
                        childnode_value = sprintf('%s %s %s',patient,study, series{k});
                        childnode_description = dicomdir.tree.(patient).(study).(series{k}).series_description_str;
                        count=count+1;
                        nodes(count) = uitreenode(childnode_value, childnode_description, childnode_iconpath, childnode_isleaf);
                    end
                elseif isfield(dicomdir.tree.(patient).(study).(series{k}),'ProtocolName'),
                    series_number = dicomdir.tree.(patient).(study).(series{k}).SeriesNumber;
                    if series_number>0,
                        childnode_value = sprintf('%s %s %s',patient, study, series{k});
                        childnode_description = dicomdir.tree.(patient).(study).(series{k}).series_description_str;
                        count=count+1;
                        nodes(count) = uitreenode(childnode_value, childnode_description, childnode_iconpath, childnode_isleaf);
                    end
                end
            end
            
      elseif findstr(value,'PATIENT'), % patient expands to study
              
            childnode_iconpath = [matlabroot, '/toolbox/matlab/icons/bookicon.gif'];
            childnode_isleaf = false;
            
            patient = value;
            matches = regexp(fieldnames( dicomdir.tree.(patient) ),'^STUDY_.+$','match');
            study_count=0;
            studies={};
            for k=1:length(matches),
                if ~isempty(matches{k}),
                    study_count=study_count+1;
                    studies{study_count} = char(matches{k});
                end
            end
            
            count = 0;
            for k=1:length(studies),
                childnode_value = sprintf('%s %s',patient, studies{k});
                childnode_description = dicomdir.tree.(patient).(studies{k}).study_description_str;
                count=count+1;
                nodes(count) = uitreenode(childnode_value, childnode_description, childnode_iconpath, childnode_isleaf);
            end
                        
      else
            disp(sprintf('Unknown node type : %s', value));
            nodes = [];
            return;
      end
           
% catch
%   error('uitree node type not recognized, may need to define an ExpandFcn for the nodes');
% end

if (count == 0)
  nodes = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dicomdir = dicomdir_struct(filename)

%% Load dicom information
try
    info = dicominfo(filename);
catch
    errordlg( sprintf('Error : %s',lasterr) );
end

%% Remove 'DirectoryRecordSequence' field, it will be rebuilt
dicomdir = rmfield(info,'DirectoryRecordSequence');

[pathstr,name,ext,version] = fileparts(dicomdir.Filename);
dicomdir.dicompath = pathstr;

%% Store patient header string parameters in dicominfo structure
dicomdir.patient_header_str.columns.name.max = 32;
dicomdir.patient_header_str.columns.birth_date.max = 10;
dicomdir.patient_header_str.columns.registration_id.max = 28;
dicomdir.patient_header_str.columns.sex.max = 3;
dicomdir.patient_header_str.columns.exam_name.max = 28;
dicomdir.patient_header_str.columns.exam_date.max = 10;
dicomdir.patient_header_str.columns.exam_time.max = 10;

dicomdir.patient_header_str.columns.name.str = 'PATIENT NAME';
dicomdir.patient_header_str.columns.birth_date.str = 'BIRTH DATE';
dicomdir.patient_header_str.columns.registration_id.str = 'REGISTRATION ID';
dicomdir.patient_header_str.columns.sex.str = 'SEX';
dicomdir.patient_header_str.columns.exam_name.str = 'EXAM NAME';
dicomdir.patient_header_str.columns.exam_date.str = 'EXAM DATE';
dicomdir.patient_header_str.columns.exam_time.str = 'EXAM TIME';

dicomdir.patient_header_str.columns.name.fieldlist = {'PatientName.FamilyName'};
dicomdir.patient_header_str.columns.birth_date.fieldlist = {'PatientBirthDate'};
dicomdir.patient_header_str.columns.registration_id.fieldlist = {'PatientID'};
dicomdir.patient_header_str.columns.sex.fieldlist = {'PatientSex'};
dicomdir.patient_header_str.columns.exam_name.fieldlist = {'STUDY_1.StudyDescription'};
dicomdir.patient_header_str.columns.exam_date.fieldlist = {'STUDY_1.StudyDate'};
dicomdir.patient_header_str.columns.exam_time.fieldlist = {'STUDY_1.StudyTime'};

dicomdir.patient_header_str.columns.name.format_func = 'format_str';
dicomdir.patient_header_str.columns.birth_date.format_func = 'format_date';
dicomdir.patient_header_str.columns.registration_id.format_func = 'format_str';
dicomdir.patient_header_str.columns.sex.format_func = 'format_str';
dicomdir.patient_header_str.columns.exam_name.format_func = 'format_str';
dicomdir.patient_header_str.columns.exam_date.format_func = 'format_date';
dicomdir.patient_header_str.columns.exam_time.format_func = 'format_time';

%% Create Series Column Header Strings
dicomdir.patient_header_str.str = '';
f = fieldnames(dicomdir.patient_header_str.columns);
for k=1:length(f),
    control_str = sprintf('%%s | %%-%ds', ...
        dicomdir.patient_header_str.columns.(f{k}).max );
    dicomdir.patient_header_str.str = sprintf(control_str, ...
        dicomdir.patient_header_str.str, ...
        dicomdir.patient_header_str.columns.(f{k}).str );
end
dicomdir.patient_header_str.str(1:2) = [];
% pad to align column headers
dicomdir.patient_header_str.str = sprintf('   %s',dicomdir.patient_header_str.str); 

%% Store series header string parameters in dicominfo structure
dicomdir.series_header_str.columns.sc.max = 4;
dicomdir.series_header_str.columns.rec.max = 3;
dicomdir.series_header_str.columns.scan_name.max = 20;
dicomdir.series_header_str.columns.orient.max = 3;
dicomdir.series_header_str.columns.stacks.max = 3;
dicomdir.series_header_str.columns.technique.max = 7;
dicomdir.series_header_str.columns.volsel.max = 0;
dicomdir.series_header_str.columns.tr.max = 8;
dicomdir.series_header_str.columns.te.max = 8;
dicomdir.series_header_str.columns.total.max = 0;
dicomdir.series_header_str.columns.allpars.max = 0;
dicomdir.series_header_str.columns.sl.max = 4;
dicomdir.series_header_str.columns.ec.max = 2;
dicomdir.series_header_str.columns.dyn.max = 4;
dicomdir.series_header_str.columns.ph.max = 2;
dicomdir.series_header_str.columns.cs.max = 0;
dicomdir.series_header_str.columns.scan_date.max = 10;
dicomdir.series_header_str.columns.scan_time.max = 10;

dicomdir.series_header_str.columns.sc.str = 'SCAN';
dicomdir.series_header_str.columns.rec.str = 'REC';
dicomdir.series_header_str.columns.scan_name.str = 'SCAN NAME';
dicomdir.series_header_str.columns.orient.str = 'ORI';
dicomdir.series_header_str.columns.stacks.str = 'STK';
dicomdir.series_header_str.columns.technique.str = 'TECHNIQ';
dicomdir.series_header_str.columns.volsel.str = 'VOL SEL';
dicomdir.series_header_str.columns.tr.str = 'TR';
dicomdir.series_header_str.columns.te.str = 'TE';
dicomdir.series_header_str.columns.total.str = 'TOTAL';
dicomdir.series_header_str.columns.allpars.str = 'ALLPARS';
dicomdir.series_header_str.columns.sl.str = 'SL';
dicomdir.series_header_str.columns.ec.str = 'EC';
dicomdir.series_header_str.columns.dyn.str = 'DYN';
dicomdir.series_header_str.columns.ph.str = 'PH';
dicomdir.series_header_str.columns.cs.str = 'CS';
dicomdir.series_header_str.columns.scan_date.str = 'SCAN DATE';
dicomdir.series_header_str.columns.scan_time.str = 'SCAN TIME';

dicomdir.series_header_str.columns.sc.fieldlist = {'SeriesNumber'};
dicomdir.series_header_str.columns.rec.fieldlist = {'SeriesNumber'};
dicomdir.series_header_str.columns.scan_name.fieldlist = {'ProtocolName','SeriesDescription'};
dicomdir.series_header_str.columns.orient.fieldlist = {'Private_2001_105f.Item_1.Private_2005_1081'};
dicomdir.series_header_str.columns.stacks.fieldlist = {'Private_2001_1060'};
dicomdir.series_header_str.columns.technique.fieldlist = {'Private_2001_1020'};
dicomdir.series_header_str.columns.volsel.fieldlist = {'VolSel'};
dicomdir.series_header_str.columns.tr.fieldlist = {'Private_2005_1030'};
dicomdir.series_header_str.columns.te.fieldlist = {'Private_2001_1025'};
dicomdir.series_header_str.columns.total.fieldlist = {'TOTAL'};
dicomdir.series_header_str.columns.allpars.fieldlist = {'AllPars'};
dicomdir.series_header_str.columns.sl.fieldlist = {'Private_2001_1018'};
dicomdir.series_header_str.columns.ec.fieldlist = {'Private_2001_1014'};
dicomdir.series_header_str.columns.dyn.fieldlist = {'Private_2001_1081'};
dicomdir.series_header_str.columns.ph.fieldlist = {'Private_2001_1017'};
dicomdir.series_header_str.columns.cs.fieldlist = {'CS'};
dicomdir.series_header_str.columns.scan_date.fieldlist = {'SeriesDate'};
dicomdir.series_header_str.columns.scan_time.fieldlist = {'SeriesTime'};

%Private_2001_107B, MRSeries Acquisition Number
%Private_2001_101D, Reconstruction Number
%Private_2001_1071, Private_2001_1072, Private_2001_1073 - Stack Angulation [AP,FH,RL]
%Private_2001_1074, Private_2001_1075, Private_2001_1076 - Stack FOV [AP,FH,RL]
%Private_2001_1078, Private_2001_1079, Private_2001_107a - Stack Offcenter [AP,FH,RL]
%Private_2005_10a3, Stack Coil ID
%Private_2005_10a4, Stack CBB Coil 1
%Private_2005_10a5, Stack CBB Coil 2
%Private_2005_10a6, Stack Channel Combi
%Private_2005_10a7, Stack Coil Conn
%Private_2005_1390, Stack Coil Function

dicomdir.series_header_str.columns.sc.format_func = 'format_sc';
dicomdir.series_header_str.columns.rec.format_func = 'format_rec';
dicomdir.series_header_str.columns.scan_name.format_func = 'format_str';
dicomdir.series_header_str.columns.orient.format_func = 'format_orient';
dicomdir.series_header_str.columns.stacks.format_func = 'num2str';
dicomdir.series_header_str.columns.technique.format_func = 'format_str';
dicomdir.series_header_str.columns.volsel.format_func = 'format_str';
dicomdir.series_header_str.columns.tr.format_func = 'format_tr';
dicomdir.series_header_str.columns.te.format_func = 'format_str';
dicomdir.series_header_str.columns.total.format_func = 'num2str';
dicomdir.series_header_str.columns.allpars.format_func = 'format_str';
dicomdir.series_header_str.columns.sl.format_func = 'num2str';
dicomdir.series_header_str.columns.ec.format_func = 'num2str';
dicomdir.series_header_str.columns.dyn.format_func = 'num2str';
dicomdir.series_header_str.columns.ph.format_func = 'num2str';
dicomdir.series_header_str.columns.cs.format_func = 'format_str';
dicomdir.series_header_str.columns.scan_date.format_func = 'format_date';
dicomdir.series_header_str.columns.scan_time.format_func = 'format_time';

%% Create Series Column Header Strings
dicomdir.series_header_str.str = '';
f = fieldnames(dicomdir.series_header_str.columns);
for k=1:length(f),
    if (dicomdir.series_header_str.columns.(f{k}).max>0 ),
        control_str = sprintf('%%s | %%-%ds', ...
            dicomdir.series_header_str.columns.(f{k}).max );
        dicomdir.series_header_str.str = sprintf(control_str, ...
            dicomdir.series_header_str.str, ...
            dicomdir.series_header_str.columns.(f{k}).str );
    end
end
dicomdir.series_header_str.str(1:2) = [];
% pad to align column headers
dicomdir.series_header_str.str = sprintf('  %s',dicomdir.series_header_str.str); 

%% Create IMAGE description strings
dicomdir.image_header_str.columns.type.max = 16;
dicomdir.image_header_str.columns.path.max = 80;
dicomdir.image_header_str.columns.size.max = 20;

dicomdir.image_header_str.columns.type.fieldlist = {'ReferencedSOPClassUIDInFile'};
dicomdir.image_header_str.columns.path.fieldlist = {'fullpath'};
dicomdir.image_header_str.columns.size.fieldlist = {'fullpath'};

dicomdir.image_header_str.columns.type.format_func = 'format_sop';
dicomdir.image_header_str.columns.path.format_func = 'format_str';
dicomdir.image_header_str.columns.size.format_func = 'format_size';

items = fieldnames(info.DirectoryRecordSequence);

patient_count=0;
study_count=0;
series_count=0;
presentation_count=0;
image_count=0;
private_count=0;
item=1;
n_items = length(items);
seeking = 'PATIENT';

% levels are: PATIENT -> STUDY -> SERIES -> (PRESENTATION | IMAGE | PRIVATE)
while item<=n_items,
    
    switch seeking,
        
        case 'PATIENT',
              
            if strcmp('PATIENT',info.DirectoryRecordSequence.(items{item}).DirectoryRecordType),
                patient_count = patient_count + 1;
                patient_field_name = sprintf('PATIENT_%d',patient_count);
                dicomdir.tree.(patient_field_name) = info.DirectoryRecordSequence.(items{item});                                         
                seeking='STUDY';
                study_count=0;
                item=item+1;
            else
                warning( sprintf('Unexpected DirectoryRecordType : %s (item=%d)', info.DirectoryRecordSequence.(items{item}).DirectoryRecordType, item) );
                item=item+1;
            end

        case 'STUDY',
            
            if strcmp('STUDY',info.DirectoryRecordSequence.(items{item}).DirectoryRecordType),
                study_count = study_count + 1;
                study_field_name = sprintf('STUDY_%d',study_count);
                dicomdir.tree.(patient_field_name).(study_field_name) = info.DirectoryRecordSequence.(items{item});               
                seeking='SERIES';
                series_count=0;
                item=item+1;
            else
                % seeking study, but found something else
                seeking='PATIENT';
            end
 
        case 'SERIES',
            
            if strcmp('SERIES',info.DirectoryRecordSequence.(items{item}).DirectoryRecordType),
                series_count = series_count + 1;
                series_number = info.DirectoryRecordSequence.(items{item}).SeriesNumber;
                series_field_name = sprintf('SERIES_x%05d_%03d', series_number, series_count);
                dicomdir.tree.(patient_field_name).(study_field_name).(series_field_name) = info.DirectoryRecordSequence.(items{item});                             
                seeking='PRESENTATION_or_IMAGE_or_PRIVATE';
                presentation_count=0;
                image_count=0;
                private_count=0;
                item=item+1;
            else
                % seeking SERIES, but found something else
                seeking='STUDY';
            end            

        case 'PRESENTATION_or_IMAGE_or_PRIVATE',
            
            if strcmp('PRESENTATION',info.DirectoryRecordSequence.(items{item}).DirectoryRecordType),
                presentation_count = presentation_count + 1;
                presentation_field_name = sprintf('PRESENTATION_%d',presentation_count);
                dicomdir.tree.(patient_field_name).(study_field_name).(series_field_name).(presentation_field_name) = info.DirectoryRecordSequence.(items{item});
                dicomdir.tree.(patient_field_name).(study_field_name).(series_field_name).(presentation_field_name).fullpath = strrep( sprintf('%s\\%s', ...
                    dicomdir.dicompath, ...
                    info.DirectoryRecordSequence.(items{item}).ReferencedFileID), '\\', '\');
                item=item+1; 
            elseif strcmp('IMAGE',info.DirectoryRecordSequence.(items{item}).DirectoryRecordType),
                image_count = image_count + 1;
                image_field_name = sprintf('IMAGE_%d',image_count);
                dicomdir.tree.(patient_field_name).(study_field_name).(series_field_name).(image_field_name) = info.DirectoryRecordSequence.(items{item});
                dicomdir.tree.(patient_field_name).(study_field_name).(series_field_name).(image_field_name).fullpath = strrep( sprintf('%s\\%s', ...
                    dicomdir.dicompath, ...
                    info.DirectoryRecordSequence.(items{item}).ReferencedFileID), '\\', '\');
                item=item+1;
            elseif strcmp('PRIVATE',info.DirectoryRecordSequence.(items{item}).DirectoryRecordType),
                private_count = private_count + 1;
                private_field_name = sprintf('PRIVATE_%d',private_count);
                dicomdir.tree.(patient_field_name).(study_field_name).(series_field_name).(private_field_name) = info.DirectoryRecordSequence.(items{item});
                dicomdir.tree.(patient_field_name).(study_field_name).(series_field_name).(private_field_name).fullpath = strrep( sprintf('%s\\%s', ...
                    dicomdir.dicompath, ...
                    info.DirectoryRecordSequence.(items{item}).ReferencedFileID), '\\', '\');
                item=item+1;  
            else
                % seeking IMAGE, but found something else
                seeking='SERIES';
            end              

        otherwise,
           % found something else
           warning( sprintf('Unexpected DirectoryRecordType : %s (item=%d)', info.DirectoryRecordSequence.(items{item}).DirectoryRecordType, item) );
           item=item+1;
           
    end
    
end

%% Order patients alphabetically
patients = fieldnames(dicomdir.tree);
for k=1:length(patients),
    patient_names{k} = upper(dicomdir.tree.(patients{k}).PatientName.FamilyName);
end

[dummy patient_sort_idx] = sort(patient_names);

for k=1:length(patients),
    sorted_patient = patients{patient_sort_idx(k)};
    sortable_patient_name = sprintf('PATIENT_x%03d',k);
    dicomdir.tree.(sortable_patient_name) = dicomdir.tree.(sorted_patient);
end

dicomdir.tree = rmfield(dicomdir.tree,patients);
dicomdir.tree = orderfields(dicomdir.tree);

%% Order series by SeriesNumber
patients = fieldnames(dicomdir.tree);
for p=1:length(patients),
    
    patient = patients{p};
    matches = regexp(fieldnames( dicomdir.tree.(patient) ),'^STUDY_.+$','match');
    study_count=0;
    studies={};
    for k=1:length(matches),
        if ~isempty(matches{k}),
            study_count=study_count+1;
            studies{study_count} = char(matches{k});
        end
    end
    
    for s=1:length(studies),
        study = studies{s};
        % order fields alphabetically
        dicomdir.tree.(patient).(study) = orderfields(dicomdir.tree.(patient).(study));
    end
    
end

%% Create description strings
patients = fieldnames(dicomdir.tree);
for p=1:length(patients),
    patient_name = patients{p};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Create patient description string
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = fieldnames(dicomdir.patient_header_str.columns);

    for k=1:length(f),
        patient_column_entries{k} = '*';
        fieldlist = dicomdir.patient_header_str.columns.(f{k}).fieldlist;
        for n=1:length(fieldlist),
            if( my_isfield(dicomdir.tree.(patient_name),fieldlist{n}) ),
                patient_column_entries{k} = my_getfield(dicomdir.tree.(patient_name),fieldlist{n});
                if isempty(dicomdir.patient_header_str.columns.(f{k}).format_func),
                    patient_column_entries{k} = my_getfield(dicomdir.tree.(patient_name),fieldlist{n});
                else
                    tmp = my_getfield(dicomdir.tree.(patient_name),fieldlist{n});
                    eval( sprintf('patient_column_entries{k} = %s(tmp);', dicomdir.patient_header_str.columns.(f{k}).format_func) );
                end
                break;
            end
        end
        tmpstr = patient_column_entries{k};
        patient_column_entries{k} = tmpstr(1:min(end,dicomdir.patient_header_str.columns.(f{k}).max));
    end

    dicomdir.tree.(patient_name).patient_description_str = '';
    for k=1:length(f),
        if ~isempty(patient_column_entries{k}),
            control_str = sprintf('%%s | %%-%ds', ...
                dicomdir.patient_header_str.columns.(f{k}).max );
            dicomdir.tree.(patient_name).patient_description_str = sprintf(control_str, ...
                dicomdir.tree.(patient_name).patient_description_str, ...
                patient_column_entries{k} );
        end
    end
    dicomdir.tree.(patient_name).patient_description_str(1:2) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    matches = regexp(fieldnames( dicomdir.tree.(patient_name) ),'^STUDY_.+$','match');
    study_count=0;
    studies={};
    for k=1:length(matches),
        if ~isempty(matches{k}),
            study_count=study_count+1;
            studies{study_count} = char(matches{k});
        end
    end

    for s=1:length(studies),
        study_name = studies{s};
        
        dicomdir.tree.(patient_name).(study_name).study_description_str = dicomdir.series_header_str.str;
        
        matches = regexp(fieldnames( dicomdir.tree.(patient_name).(study_name) ),'^SERIES_.+$','match');
        series_count=0;
        series={};
        for k=1:length(matches),
            if ~isempty(matches{k}),
                series_count=series_count+1;
                series{series_count} = char(matches{k});
            end
        end

        for ss=1:length(series),
            
            series_name = series{ss};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Create series description string
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f = fieldnames(dicomdir.series_header_str.columns);

            for k=1:length(f),
                series_column_entries{k} = '*';
                fieldlist = dicomdir.series_header_str.columns.(f{k}).fieldlist;
                for n=1:length(fieldlist),
                    if( my_isfield(dicomdir.tree.(patient_name).(study_name).(series_name),fieldlist{n}) ),
                        if isempty(dicomdir.series_header_str.columns.(f{k}).format_func)
                            series_column_entries{k} = my_getfield(dicomdir.tree.(patient_name).(study_name).(series_name),fieldlist{n});
                        else
                            tmp = my_getfield(dicomdir.tree.(patient_name).(study_name).(series_name),fieldlist{n});
                            eval( sprintf('series_column_entries{k} = %s(tmp);', dicomdir.series_header_str.columns.(f{k}).format_func) );
                        end
                        break;
                    end
                end
                tmpstr = series_column_entries{k};
                series_column_entries{k} = tmpstr(1:min(end,dicomdir.series_header_str.columns.(f{k}).max));
            end

            dicomdir.tree.(patient_name).(study_name).(series_name).series_description_str = '';
            for k=1:length(f),
                if ~isempty(series_column_entries{k}),
                    control_str = sprintf('%%s | %%-%ds', ...
                        dicomdir.series_header_str.columns.(f{k}).max );
                    dicomdir.tree.(patient_name).(study_name).(series_name).series_description_str = ...
                        sprintf(control_str, ...
                        dicomdir.tree.(patient_name).(study_name).(series_name).series_description_str, ...
                        series_column_entries{k} );
                end
            end
            dicomdir.tree.(patient_name).(study_name).(series_name).series_description_str(1:2) = [];               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            image_count=0;
            images={};
            
            matches = regexp(fieldnames( dicomdir.tree.(patient_name).(study_name).(series_name) ),'^IMAGE_.+$','match');
            for k=1:length(matches),
                if ~isempty(matches{k}),
                    image_count=image_count+1;
                    images{image_count} = char(matches{k});
                end
            end
 
            matches = regexp(fieldnames( dicomdir.tree.(patient_name).(study_name).(series_name) ),'^PRIVATE_.+$','match');
            for k=1:length(matches),
                if ~isempty(matches{k}),
                    image_count=image_count+1;
                    images{image_count} = char(matches{k});
                end
            end

            matches = regexp(fieldnames( dicomdir.tree.(patient_name).(study_name).(series_name) ),'^PRESENTATION_.+$','match');
            for k=1:length(matches),
                if ~isempty(matches{k}),
                    image_count=image_count+1;
                    images{image_count} = char(matches{k});
                end
            end            
            
            for ii=1:length(images),
                
                image_name = images{ii};
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Create IMAGE description string
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                f = fieldnames(dicomdir.image_header_str.columns);

                for k=1:length(f),
                    image_column_entries{k} = '*';
                    fieldlist = dicomdir.image_header_str.columns.(f{k}).fieldlist;
                    for n=1:length(fieldlist),
                        if( my_isfield(dicomdir.tree.(patient_name).(study_name).(series_name).(image_name),fieldlist{n}) ),
                            tmp = my_getfield(dicomdir.tree.(patient_name).(study_name).(series_name).(image_name),fieldlist{n});
                            eval( sprintf('image_column_entries{k} = %s(tmp);', dicomdir.image_header_str.columns.(f{k}).format_func) );
                            break;
                        end
                    end
                    tmpstr = image_column_entries{k};
                    image_column_entries{k} = tmpstr(1:min(end,dicomdir.image_header_str.columns.(f{k}).max));
                end

                dicomdir.tree.(patient_name).(study_name).(series_name).(image_name).image_description_str = '';
                for k=1:length(f),
                    if ~isempty(image_column_entries{k}),
                        control_str = sprintf('%%s | %%-%ds', ...
                            dicomdir.image_header_str.columns.(f{k}).max );
                        dicomdir.tree.(patient_name).(study_name).(series_name).(image_name).image_description_str = ...
                            sprintf(control_str, ...
                            dicomdir.tree.(patient_name).(study_name).(series_name).(image_name).image_description_str, ...
                            image_column_entries{k} );
                    end
                end
                dicomdir.tree.(patient_name).(study_name).(series_name).(image_name).image_description_str(1:2) = [];
                
            end % end image
        end % end series
    end % end study
end % end patient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function childnodes = autoexpand(tree,parentnode)

childnodes = myExpfcn(tree,parentnode.getValue);
if (length(childnodes)==1),
    % Then we dont have an array of nodes. Create an array.
    tmpnodes = childnodes;
    childnodes = javaArray('com.mathworks.hg.peer.UITreeNode', 1);
    childnodes(1) = java(tmpnodes);
end
tree.add(parentnode, childnodes);
tree.setLoaded(parentnode, true);
tree.expand(parentnode);

% autoexpand study
if (length(childnodes)>=1),
    value = childnodes(1).getValue;
    if ( findstr(value,'STUDY') & isempty(findstr(value,'SERIES')) ),
        series_nodes = autoexpand(tree,childnodes(1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = my_isfield(s,fieldname)
try
    [x,y]=strtok(fieldname(end:-1:1),'.');
    if ~isempty(y),
        lower_fieldname = x(end:-1:1);
        upper_fieldname = y(end:-1:1); 
        eval( sprintf('result = isfield(s.%s,''%s'');',upper_fieldname(1:end-1),lower_fieldname) );
    else
        result = isfield(s,fieldname);
    end
catch
    result=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = my_getfield(s,fieldname)
[x,y]=strtok(fieldname(end:-1:1),'.');
if ~isempty(y),
    lower_fieldname = x(end:-1:1);
    upper_fieldname = y(end:-1:1); 
    eval( sprintf('result = getfield(s.%s,''%s'');',upper_fieldname(1:end-1),lower_fieldname) );
else
    result = getfield(s,fieldname);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_str(value)
str = value;
if isempty(str),
    str = '*';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_sc(value)
str = num2str(floor(value/100));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_rec(value)
str = num2str(mod(value,100));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_orient(value)
switch value,
    case 'FH',
        str = 'TRA';
    case 'RL',
        str = 'SAG';
    case 'AP',
        str = 'COR';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_date(value)
str = sprintf('%s/%s/%s', value(5:6), value(7:8), value(1:4) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_tr(value)
str='';
for k=1:length(value),
    str = sprintf('%s/%s',str,num2str(value(k)));
end
str(1)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_sop(value)

switch value,
    case '1.2.840.10008.5.1.4.1.1.4',
        str = 'MR CLASSIC';
    case '1.2.840.10008.5.1.4.1.1.4.1',
        str = 'MR ENHANCED';
    case '1.2.840.10008.5.1.4.1.1.4.2',
        str = 'MR SPECTRO';
    case '1.2.840.10008.5.1.4.1.1.7',
        str = 'SECON CAPT';
    case '1.2.840.10008.5.1.4.1.1.11.1',
        str = 'GRAYSC PS';
    case '1.2.840.10008.5.1.4.1.1.66',
        str = 'SCAN PARS';
    case '1.3.46.670589.11.0.0.12.1',
        str = 'PRIV SPECT';
    case '1.3.46.670589.11.0.0.12.2',
        str = 'PRIV DATA';
    otherwise,
        str = 'UNKNOWN';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_size(value)
d = dir(value);
if ~isempty(d),
    byte_str = num2str(d.bytes);
    count=0;
    byte_str = byte_str(end:-1:1);
    for k=1:length(byte_str),
        count=count+1;
        comma_byte_str(count) = byte_str(k);
        if mod(k,3)==0,
            count=count+1;
            comma_byte_str(count) = ',';
        end
    end
    str = comma_byte_str(end:-1:1);
    if str(1)==',',
        str(1) = [];
    end
    str = sprintf('%s bytes',str);
    str = sprintf('%+20s', str);
else
    str = 'NOT FOUND';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = format_time(value)
if length(value)>6,
    str = sprintf('%s:%s:%s.%s', value(1:2), value(3:4), value(5:6), value(8) );
else
    str = sprintf('%s:%s:%s', value(1:2), value(3:4), value(5:6) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function childnodes = myNodeExpanded(src, evd, tree, expfcn)                           

evdnode  = evd.getCurrentNode;

if ~tree.isLoaded(evdnode)
    value = evdnode.getValue;

    % <call a user function(value) which returns uitreenodes>;
    cbk = expfcn;
    if iscell(cbk)
        childnodes = feval(cbk{1}, tree, value, cbk{2:end});
    else
        childnodes = feval(cbk, tree, value);
    end

    if (length(childnodes) == 1)
        % Then we dont have an array of nodes. Create an array.
        chnodes = childnodes;
        childnodes = javaArray('com.mathworks.hg.peer.UITreeNode', 1);
        childnodes(1) = java(chnodes);
    end

    tree.add(evdnode, childnodes);
    tree.setLoaded(evdnode, true);
    
    % autoexpand study
    if (length(childnodes)>=1),
        value = childnodes(1).getValue;
        if ( findstr(value,'STUDY') & isempty(findstr(value,'SERIES')) ),
            series_nodes = autoexpand(tree,childnodes(1));
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function type = nodetype(value)   

type = '';
if strcmp(value,'root'),
    type = 'root';
elseif findstr(value,'PRIVATE'),
    type = 'PRIVATE';     
elseif findstr(value,'PRESENTATION'),
    type = 'PRESENTATION';    
elseif findstr(value,'IMAGE'),
    type = 'IMAGE';    
elseif findstr(value,'SERIES'),
    type = 'SERIES';
elseif findstr(value,'STUDY'),
    type = 'STUDY';
elseif findstr(value,'PATIENT'),
    type = 'PATIENT';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = myMouseClicked(src, evd, tree)

mtree = tree;
jtree = mtree.Tree;
dicomdir = get(mtree.FigureComponent,'UserData');

try
    selected_node = jtree.getLastSelectedPathComponent();
    selected_node_value = selected_node.getValue();
catch
    return;
end

matches = regexp(selected_node_value,'\w*','match');
if length(matches)==3,
    patient_name = matches{1};
    study_name = matches{2};
    series_name = matches{3};
else % not clicking on a series
    return;
end

matches = regexp(fieldnames( dicomdir.tree.(patient_name).(study_name).(series_name) ),'^IMAGE_.+$','match');
image_count=0;
for k=1:length(matches),
    if ~isempty(matches{k}),
        image_count=image_count+1;
        images{image_count} = char(matches{k});
    end
end

dicom_filename = '';
for k=1:image_count,
    sop = dicomdir.tree.(patient_name).(study_name).(series_name).(images{k}).ReferencedSOPClassUIDInFile;
    if strcmp(sop,'1.2.840.10008.5.1.4.1.1.4.1'),
        dicom_filename = dicomdir.tree.(patient_name).(study_name).(series_name).(images{k}).fullpath;
    end
end
dicomdir.dicom_filename = dicom_filename;
set(mtree.FigureComponent,'UserData',dicomdir);

click_count = get(evd,'ClickCount');
if click_count>1,
    uiresume(dicomdir.h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
