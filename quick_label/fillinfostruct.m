% for quick_label abdomen project

datapath = '~/mounts/fs2/masi/bryanfw/abdominal_data/';
origdatapath = [datapath 'ImgRaw/'];
wadepath = [datapath 'wade/'];
rebepath = [datapath 'rebe/'];

origdatadir = dir([origdatapath '*.nii.gz']);
wadedir = dir([wadepath '*.nii.gz']);
rebedir = dir([rebepath '*.nii.gz']);

rebedir = struct2cell(rebedir); rebedir = rebedir(1,:);
wadedir = struct2cell(wadedir); wadedir = wadedir(1,:);

for ii = 1:length(orig_data);
    obs.intensity{ii} = [origdatapath origdatadir(ii).name];
    id = regexp(origdatadir(ii).name,'_');
    id = origdatadir(ii).name(1:(id(1)-1));
    
    % find corresponding rebecca-rated labels 
    ind = find(not(cellfun('isempty',strfind(rebedir,id))));
    if ~isempty(ind); 
        obs.rebetruth{ii} = [rebepath rebedir{ind}];
    end
    
    % find corresponding wade-rated labels
    % wade rated some twice, indicated by wade1 or wade2 suffix
    ind = not(cellfun('isempty',strfind(wadedir,id))); 
    ind1 = find(ind & not(cellfun('isempty',strfind(wadedir,'wade1'))));
    ind2 = find(ind & not(cellfun('isempty',strfind(wadedir,'wade2'))));
    
    if ~isempty(ind1)
        obs.wade1truth{ii} = [wadepath wadedir{ind}];
    end
    if ~isempty(ind2)
        obs.wade2truth{ii} = [wadepath wadedir{ind2}];
    end   
end

if length(obs.wade1truth)<length(obs.intensity)
    obs.wade2truth{length(obs.intensity)}=[];
end
if length(obs.wade2truth)<length(obs.intensity)
    obs.wade1truth{length(obs.intensity)}=[];
end
if length(obs.rebetruth)<length(obs.intensity)
    obs.rebetruth{length(obs.intensity)}=[];
end

havewade1 = not(cellfun('isempty',obs.wade1truth));
havewade2 = not(cellfun('isempty',obs.wade2truth)); 
haverebe  = not(cellfun('isempty',obs.rebetruth));  

obs.havetruth = sum([havewade1;havewade2;haverebe]);

