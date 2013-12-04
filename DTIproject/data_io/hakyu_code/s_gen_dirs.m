
%[s_gen_dirs] generate directories for saving recon related files
%
% USAGE:
%   s_gen_dirs
%
% OUTPUT:
%   saveDataDir_s:  Post processed files are saved intermediate to recon
%   saveReconDir_s: Final recon files are saved 
%
%
% Last modified
%  2010.08.09.
%  2012.05.17.
%   Add another term at the end of dataDir_s name usign
%   DATAFLAGparams.addTermAtDataDir.
%
% Ha-Kyu



%% Generate saveDir_s and saveReconDir_s and sharedDir_s

% Get saveDir_s in which common data is saved for each DWI scan.
saveDirName_s = sprintf('%s_R%dE%dS%d', ...
    DWIparams.filename,DWIparams.SENSE_FACTOR,DWIparams.EPI_FACTOR(1), ...
    DWIparams.nSHOT);
savePathName_s = [loadDir_s];

if exist([savePathName_s,saveDirName_s],'dir')==0
    saveDataDir_s = [savePathName_s saveDirName_s];
    success_id = mkdir(saveDataDir_s);
    if success_id==1
        fprintf('   [%s]\n',saveDataDir_s)
        fprintf('   is generated successfully.\n')
    else
        error('s_gen_dirs:main','Check pathname.')
    end
else
    saveDataDir_s = [savePathName_s saveDirName_s];
end

if ~isempty(DATAFLAGparams.addTermAtDataDir)
    saveDataDir_s = [saveDataDir_s '_' DATAFLAGparams.addTermAtDataDir];
end


% Get saveReconDir_s in which all recon related data are saved.
saveReconDir_s = sprintf('%s%srecon',saveDataDir_s,filesep);
if exist(saveReconDir_s,'dir')==0
    success_id = mkdir(saveReconDir_s);
    if success_id==1
        fprintf('   [%s]\n',saveReconDir_s)
        fprintf('   is generated successfully for recon data.\n')
    else
        error('s_gen_dirs:main','Check pathname.')
    end    
end

% Generate directory for sharedDir_s, if it doesn't exist.
if exist(sharedDir_s,'dir')==0
    success_id = mkdir(sharedDir_s);
    if success_id==1
        fprintf('   [%s]\n',sharedDir_s)
        fprintf('   is generated successfully for shared data.\n')
    else
        error('s_gen_dirs:main','Check pathname.')
    end    
end

% Report.
fprintf('[saveDataDir_s] is,\n')
fprintf('    [%s]\n',saveDataDir_s)
fprintf('[saveReconDir_s] is,\n')
fprintf('    [%s]\n',saveReconDir_s)
fprintf('\n')













