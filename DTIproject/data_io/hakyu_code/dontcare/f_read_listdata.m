
function [info,varargout] = f_read_listdata(name,info,varargin)
%[f_read_listdata] reads .LIST/.DATA data.
%
% USAGE:
%   [info,STD,PHX,FRX,NOI,PHC] =
%   f_read_listdata('C:\data\20050308\PROPELLER DATA\raw_777.list'); 
%
% INPUT:
%   Full raw file or raw file name to read.
%   Loading options can be used such as,
%       ky      : ky position
%       kz      : kz position
%       loca    : slice or stack position
%       card    : cardiac phase
%       echo    : echoes
%       mix     : mix
%       dyn     : dynamics
%       aver    : averages or NSA
%       chan    : coils or channels
%       extr1   : dw or not (0,1)
%       extr2   : dw orientation
%
% OUTPUT:
%   info:   Tag for all data
%   STD:    Standard type data
%   PHX:    PHX type data
%   FRX:    FRX type data
%   NOI:    NOI type data
%   PHC:    PHC type data
%
% NOTE:
%   Bodycoil is 27th (3T, SENSE-Head-8) and 26th (7T, RX-Intf-1).
%
%
% Last modified
% 2010.08.06.
%   Add start_sign.PHC_start_sign and start_sign.STD_start_sign as output.
%   But it is commented.
% 2010.08.09.
%   This function is generated from [readListData_HKJ4.m] which is
%   originally developed by EB Welch.
% 2010.09.07.
%   Selectively load LIST/DATA data based on loadopts.
% 2010.09.08.
%   Remove unused idx***.
% 2010.09.09.
%   Remove unused comments and scripts.
% 2010.09.20.
%   Final version of [test_f_read_listdata.m]. It is saved as
%   [f_read_listdata.m] here. This can use loadopts for loading options.
% 2010.12.15.
%   Modify 'Check point #5' for faster routine.
% 2010.12.16.
%   Get listtext as input.



%% Read .LIST/.DATA files
disp('  ________________________')

% Get list/data filename.
toks = regexp(name,'^(.*?)(\.list|\.data)?$','tokens');
prefix = toks{1}{1};
% listname = sprintf('%s.list',prefix);
dataname = sprintf('%s.data',prefix);

% Read list file.
% fid = fopen(listname,'r');
% if fid~=-1,
%     listtext = fread(fid,inf,'uint8=>char')';
%     fclose(fid);
% else
%     error( sprintf('cannot open %s for reading', listname) );
% end
% fprintf('      LIST file is read\n')

% Set attributes.
attr = {'mix','dyn','card','echo','loca','chan','extr1','extr2','ky', ...
    'kz','n.a.','aver','sign','rf','grad','enc','rtop','rr','size','offset'};

% Clean attr names (will be used as variablenames and fieldnames)
for k=1:length(attr),
    attr{k} = cleanFieldname( attr{k} );
end

% Set pattern for tags.
% pattern = '(?<typ>\w+)';
% for k=1:length(attr),
%     pattern = sprintf('%s\\s+(?<%s>-?\\d+)',pattern,attr{k});
% end


%----- Check point #2 ----- start
% tic
% info = regexp(listtext, pattern, 'names');
% t = toc;
% fprintf('      CP#2: %.3f sec\n',t)
%----- Check point #2 ----- end



%% Read STD, PHX, FRX, NOI, PHC data
nout = nargout - 1;
for ind_typ = 1:nout   % 1:STD, 2:PHX, 3:FRX, 4:NOI, 5:PHC
    
    if ind_typ==1        
        idxTYP = find( strcmp({info(:).typ},'STD') );        
        typname_s = 'STD';
    elseif ind_typ==2
        idxTYP = find( strcmp({info(:).typ},'PHX') );
        typname_s = 'PHX';
    elseif ind_typ==3
        idxTYP = find( strcmp({info(:).typ},'FRX') );
        typname_s = 'FRX';
    elseif ind_typ==4
        idxTYP = find( strcmp({info(:).typ},'NOI') );
        typname_s = 'NOI';
    elseif ind_typ==5
        idxTYP = find( strcmp({info(:).typ},'PHC') );
        typname_s = 'PHC';
    end
    fprintf('      READING [%s] data\n',typname_s)
    fprintf('        length(idxTYP) = %d\n',length(idxTYP))
    
    
    % Pass ind_typ when there is no corresponding data.
    if isempty(idxTYP)
        continue;
    else
        
        %----- Check point #4 ----- start
        tic
        for k=1:length(attr)
            eval( sprintf('%s = sort( str2num( char( unique( {info(idxTYP).%s} ) ) ) );', ...
                attr{k},attr{k}) );
        end
        t = toc;
        fprintf('      CP#4: %.3f sec\n',t)
        %----- Check point #4 ----- end
        
        
        % DATA is a multi-dimensional array organized as
        if ind_typ==1
            order = {'kx','ky','kz','loca','dyn','card','echo','mix', ...
                'aver','chan','extr1','extr2'};
        elseif ind_typ==2
            order = {'kx','loca','echo','mix','chan','sign'};
        elseif ind_typ==3
            order = {'kx','loca','echo','mix','chan','sign'};
        elseif ind_typ==4
            order = {'kx','loca','chan'};
        elseif ind_typ==5
            order = {'kx','ky','kz','loca','dyn','card','echo','mix', ...
                'chan','aver','sign','grad'};
        end
        
        
        %------------------------------------------------- start
        % 2010.09.07.
        % Use only optional load. This must be done only for 'STD'
        if nargin > 2
            if strcmpi(typname_s,'STD')
                % Parse loadopts.
                p = inputParser;
                p.StructExpand = true;
                p.CaseSensitive = true;
                p.KeepUnmatched = false; % throw an error for unmatched inputs
                p.addRequired('name', @ischar);
                for k=1:length(order)
                    p.addParamValue(order{k}, [], @isnumeric);
                end
                p.parse(name, varargin{:});
                loadopts = rmfield(p.Results,'name');
                
                fprintf('      loadopts parsed\n')
                
                
                %----- Check point #5 ----- start
                tic
                % Update idxTYP corresponding to loadopts.                
                for k = 1:length(order)
                    if ~isempty(loadopts.(order{k}))
                        
                        % Slower routine.
%                         idxtyp = [];                        
%                         eval(sprintf('%s_hold = %s;',order{k},order{k}))
%                         eval(sprintf('val_first = %s_hold(1);',order{k}))
%                         tmp = val_first;
%                         for ind=1:length(idxTYP)
%                             val = str2double(info(idxTYP(ind)).(order{k}));
%                             if any(val==(loadopts.(order{k})+tmp))
%                                 idxtyp = [idxtyp,idxTYP(ind)];
%                             end
%                         end
%                         eval(sprintf('idxtyp_%s = idxtyp;',order{k}))
%                         clear  idxtyp
                        
                        % Faster routine: 2010.12.15
                        tt = clock;
                        eval(sprintf('%s_hold = %s;',order{k},order{k}))
                        eval(sprintf('val_first = %s_hold(1);',order{k}))
                        tmp = val_first;
                        val_v = str2double({info(idxTYP).(order{k})});
                        ind_v = (val_v==(loadopts.(order{k})+tmp)); % logical indexing
                        idxtyp = idxTYP(ind_v);
                        eval(sprintf('idxtyp_%s = idxtyp;',order{k}))
                        fprintf('        calculated idxtyp for %s in %f sec\n', ...
                            order{k},etime(clock,tt))
                        clear  idxtyp  tt  val_v  ind_v     
                        
                    end                    
                end
                idxtyp_c = {};
                ind = 0;
                for k = 1:length(order)
                    if ~isempty(loadopts.(order{k}))
                        ind = ind+1;
                        idxtyp_c{ind} = eval(sprintf('idxtyp_%s',order{k}));
                    end
                end
                idxtyp = idxtyp_c{1};
                if length(idxtyp_c) > 1
                    for ind = 2:length(idxtyp_c)
                        idxtyp = intersect(idxtyp,idxtyp_c{ind});
                    end
                end
                idxTYP = idxtyp;                
                clear idxtyp;
                t = toc;
                fprintf('      CP#5: %.3f sec\n',t)
                
                % Check if updated idxTYP is ok.
                %disp({info(idxTYP).(order{k})) must show
                %(loadopts.(order{k})+1)th value in order{k} of loadopts.
                
                %----- Check point #5 ----- end
                
                
                % Regenerate attr values based on loadopts.
                for k=1:length(attr)
                    ind_cell = strfind(order,attr{k});
                    for n = 1:length(ind_cell)
                        if ind_cell{n}==1
                            indfound = n;
                            break
                        else
                            indfound = [];
                        end
                    end
                    if ~isempty(indfound)
                        if ~isempty(loadopts.(order{indfound}))
                            % Take attr values only intersecting with loadopts.
                            eval(sprintf('val = %s;',order{indfound}));
                            if val(1) < 0
                                % Index starts from -ve ~ +ve such as ky,
                                % kz etc.
                                eval(sprintf('%s = val(loadopts.(order{indfound})+1);', ...
                                    order{indfound}))
                            else
                                % Index starts from 0 ~ +ve such as chan.
                                eval(sprintf('%s = loadopts.(order{indfound});', ...
                                    order{indfound}))
                            end
                        end
                    end
                end % for k=1
                
            end
        end
        %------------------------------------------------- end
        
        
        % Prepare the size of raw data to read.
        tmpstr = '';
        for k=1:length(order),
            if strcmp(order{k},'kx')==1,
                tmpstr = sprintf('%s max(size)/4/2',tmpstr);
            else
                %tmpstr = sprintf('%s length(%s)',tmpstr,order{k});
                if strcmp(order{k},'sign')~=1
                    
                    %--------------------------------------- start
                    if (exist('loadopts','var')==1) && ...
                            (strcmpi(typname_s,'STD')) && ...
                            (~isempty(loadopts.(order{k})))
                        
                        if strcmpi(order{k},'chan')==1
                            % bodycoil is 27th(3T) or 26th(7T), then take
                            % care of this by using length().
                            tmpstr = sprintf('%s length(chan)',tmpstr);
                        else
                            %tmpstr = sprintf('%s max(%s)-min(%s)+1',tmpstr,order{k},order{k});
                            tmpstr = sprintf('%s length(%s)',tmpstr,order{k});
                        end
                    else
                        if strcmpi(order{k},'chan')==1
                            % bodycoil is 27th(3T) or 26th(7T), then take
                            % care of this by using length().
                            tmpstr = sprintf('%s length(chan)',tmpstr);
                        else
                            tmpstr = sprintf('%s max(%s)-min(%s)+1',tmpstr,order{k},order{k});
                        end
                    end
                    %--------------------------------------- end
                    
                else
                    tmpstr = sprintf('%s max(%s)-min(%s)',tmpstr,order{k},order{k});
                end
            end
        end
        
        
        % Preserve raw data space.
        % tmpdata is complex then data_%s must be complex.
        %eval( sprintf('data_%s = zeros([%s],''single'');', typname_s,tmpstr) );
        eval( sprintf('data_%s = complex(zeros([%s],''single''),zeros([%s],''single''));', ...
            typname_s,tmpstr,tmpstr) ); 
        clear  size
        fprintf('      data_%s is generated in size of,\n',typname_s)
        fprintf('        ')
        eval(sprintf('disp(size(data_%s))',typname_s))
        pack
        
        
        %----- Check point #7 ----- start
        % Change the index to start from 1.
        tic        
        for k=1:length(order),
            if strcmp(order{k},'kx')==0,
                
                %--------------------------------------- start
                if (exist('loadopts','var')==1) && ...
                        (strcmpi(typname_s,'STD')) && ...
                        (~isempty(loadopts.(order{k})))
                    
%                     if strcmpi(order{k},'chan')
%                         % bodycoil is 27th(3T) or 26th(7T), then take care
%                         % of this.
%                         tmp = 1;
%                         eval( sprintf('%s(:) = %s(:) + tmp;',order{k},order{k}) );
%                     else
%                         eval( sprintf('tmp = 1 - min(%s(:));',order{k}) );
%                         eval( sprintf('%s(:) = %s(:) + tmp;',order{k},order{k}) );
%                     end
                    
                    eval( sprintf('tmp = 1 - min(%s);',sprintf('%s_hold',order{k})) );
                else
                    eval( sprintf('tmp = 1 - min(%s(:));',order{k}) );
                    eval( sprintf('%s(:) = %s(:) + tmp;',order{k},order{k}) );
                end
                %--------------------------------------- end
                
                for j=1:length(idxTYP)
                    if strcmp(order{k},'sign')==1
                        signval = str2double(info(idxTYP(j)).(order{k}));
                        if signval == -1
                            tmp = 3;
                        elseif signval == 1
                            tmp = 0;
                        else
                            error('Unknown sign value.')
                        end
                    end
                    
                    % From this on, every order{k} index starts from 1.
                    info(idxTYP(j)).(order{k}) = ...
                        str2double( info(idxTYP(j)).(order{k}) ) + tmp;
                end
            end
        end
        t = toc;
        fprintf('      CP#7: %.3f sec\n',t)
        
        % Check if info(idxTYP).(order{k}) is generated correctly.
        %info(idxTYP).(order{k}) must show (loadopts.(order{k})+1)th
        %element in kz+tmp array. For example, in kz case with k = 3,
        %{info(idxTYP).(order{k})} must show all 3 or 5 in case
        %kz=[-3,-2,-1,0,1,2], tmp=1-(-3)=4 and loadopts.kz = [2,4].
                
        %----- Check point #7 ----- end
        
        
        
        %----- Check point #8 ----- start        
        fid = fopen(dataname,'r','ieee-le');
        if fid==-1,
            error( sprintf('cannot open %s for reading', listname) );
        end
                
        if strcmpi(computer,'PCWIN') || strcmpi(computer,'PCWIN64')
            flag_pcwin=1;
            hwait = waitbar(0,'=============================================');
            set( get( findobj(hwait,'type','axes'),'Title') ,'Interpreter','none');
            set( get( findobj(hwait,'type','axes'),'Title') ,'String', ...
                sprintf('Reading raw data file, typ %s ...', typname_s) );
        else
            flag_pcwin=0;
        end
        
        tic
        N = length(idxTYP);
        fprintf('        length(idxTYP) = %d at CP#8\n',N)
        for n=1:N,
            if ( fseek(fid, str2num( info(idxTYP(n)).offset ) ,'bof') == 0)
                
                % Read raw data.
                tmpdata = fread(fid, str2num( info(idxTYP(n)).size )/4 ,'float32');
                if n==1
                    t0=clock;
                elseif n==4000
                    t1=etime(clock,t0);
                    fprintf('          4000 tmpdata read in %f sec\n',t1)
                end                
                
                tmpdata = tmpdata(1:2:end) + 1i*tmpdata(2:2:end);
                if ind_typ==2 || ind_typ==3 || ind_typ==5
                    if info(idxTYP(n)).sign==1
                        signval = 1;
                    elseif info(idxTYP(n)).sign==2
                        signval = -1;
                    else
                        error('Unknown signval.')
                    end
                    tmpdata = tmpdata * signval;
                else
                    tmpdata = tmpdata * str2num( info(idxTYP(n)).sign );
                end
                
                % Make single complex data: size of double.
                tmpdata = single(tmpdata);
                
                
                % Generate index to read into the reserved raw data space
                % using updated idxTYP.
                tmpstr='';
                for k=1:length(order),
                    if strcmp(order{k},'kx')==1,
                        tmpstr = sprintf('%s,1:%d', tmpstr, length(tmpdata));
                    else
                        idx = info(idxTYP(n)).(order{k});
                        
                        %---------------------------------------- start
                        if (exist('loadopts','var')==1) && ...
                                (strcmpi(typname_s,'STD')) && ...
                                (~isempty(loadopts.(order{k})))
                            
                            % idx must be one of val.
                            val = loadopts.(order{k})+1; % index from 1
                            ind = find(idx==val);
                            idx = ind;
                            tmpstr = sprintf('%s,%d', tmpstr, idx);
                        else
                            % 'chan' case.
                            % This is to take care of bodycoil index, when
                            % the index is larger than actual number of
                            % coils. For example, 3T 'SENSE-Head-8' coil
                            % has 8 coils (0~7) and the bodycoil index is
                            % 27 (0~). Then the 'idx' is 28 and the
                            % reserved datasize is 28 along the coil
                            % dimension, though there are 9 coils actually
                            % including surface and bodycoil.
                            if strcmpi(order{k},'chan')
                                nchan = length(chan);
                                if idx==chan(nchan)                                    
                                    idx = chan(nchan-1)+1;
                                end
                            end
                            tmpstr = sprintf('%s,%d', tmpstr, idx);
                        end
                        %---------------------------------------- end
                        
                    end
                end                
                tmpstr(1)=[]; % Delete initial comma
                                
                % Read raw data: data_%s is double and tmpdata is single(complex), hence double.
                eval( sprintf('data_%s(%s) = tmpdata;', typname_s,tmpstr) );
                
            else
                error('Cannot FSEEK to offset=[%f] in data file', ...
                    info(idxTYP(k)).offset);
            end
            if n==4000, fprintf('        '), end
            if mod(n,4000)==0
                fprintf('.')
                if flag_pcwin==1, waitbar(n/N,hwait); end
            end            
        end % for n=1:N        
        fclose(fid);
        if flag_pcwin==1, close(hwait); end
        t = toc;
        fprintf('\n')
        fprintf('      CP#8: %.3f sec\n',t)
        %----- Check point #8 -----
        
        % Output.        
        varargout{ind_typ} = eval(sprintf('data_%s',typname_s));
        eval(sprintf('clear  data_%s',typname_s))
        
    end % if isempty(idxTYP)
    
    % Clear.
    clear  idxTYP
        
end % for ind_typ
disp('  ________________________')


%% NOTE
% # Complex data vector types:
% # --------------------------
% # STD = Standard data vector (image data or spectroscopy data)
% # REJ = Rejected standard data vector
% #       (only for scans with arrhythmia rejection)
% # PHX = Correction data vector for EPI/GraSE phase correction
% # FRX = Correction data vector for frequency spectrum correction
% # NOI = Preparation phase data vector for noise determination
% # NAV = Phase navigator data vector
% #
% # Other attributes of complex data vectors:
% # -----------------------------------------
% # mix    = mixed sequence number
% # dyn    = dynamic scan number
% # card   = cardiac phase number
% # echo   = echo number
% # loca   = location number
% # chan   = synco channel number
% # extr1  = extra attribute 1 (semantics depend on type of scan)
% # extr2  = extra attribute 2 (   ''       ''   ''  ''  ''  '' )
% # kx,ky  = k-space coordinates in 1st and 2nd preparation direction (spectroscopy data)
% # ky,kz  = k-space coordinates in 1st and 2nd preparation direction (image data)
% # aver   = sequence number of this signal average
% # sign   = sign of measurement gradient used for this data vector (1 = positive, -1 = negative)
% # rf     = sequence number of this rf echo (only for TSE, TFE, GraSE)
% # grad   = sequence number of this gradient echo (only for EPI/GraSE)
% # enc    = encoding time (only for EPI/GraSE)
% # rtop   = R-top offset in ms
% # rr     = RR-interval length in ms
% # size   = data vector size   in number of bytes (1 complex element = 2 floats = 8 bytes)
% # offset = data vector offset in number of bytes (first data vector starts at offset 0)
% #
% # The complex data vectors are represented as binary data in little endian single precision IEEE float format.
% #
% # Please note that complex data vector attributes which are not relevant for a certain type of vector
% # may have arbitrary values!

% # Identifying attributes of complex data vectors:
% # -----------------------------------------------
% # The next table specifies the identifying attributes for each type of complex data vector:
% #
% # typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    aver  sign  rf    grad  enc   rtop  rr    size   offset
% # --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ------ ------
% #
% # STD   *     *     *     *     *     *     *     *     *     *     *     *     *     *
% # REJ   *     *     *     *     *     *     *     *     *     *     *     *     *     *
% # PHX   *                 *     *     *                                   *     *     *
% # FRX   *                 *     *     *                                   *            
% # NOI                           *     *                                                
% # NAV   *     *     *     *     *     *     *     *     *     *     *     *     *     *



%% Subfunctions
function [s] = cleanFieldname(s)

illegal_chars = {...
    '+','-','*','.',...
    '^','\','/','.',...
    '=','~','<','>',...
    '&','|',':',';',...
    '(',')','{','}',...
    '[',']','{','}',...
    '''','%',' ','!', ...
    '@','#','$','`',...
    '?',',','"',...
    };

general_replacement_char = '_';
firstchar_replacement_char = 'x'; % cannot be an underscore

for k=1:length(illegal_chars),
    s = strrep(s,illegal_chars{k},general_replacement_char);
end

% first character cannot be a number
firstchar_code = double(s(1));
if ( (firstchar_code>=double('0')) && (firstchar_code<=double('9')) )
    s(1) = firstchar_replacement_char;
end

% first character cannot be an underscore
if(s(1)=='_'),
    s(1) = firstchar_replacement_char;
end









