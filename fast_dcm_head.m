function sparsehead = fast_dcm_head(dcmdir,tags,N)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% fact_dcm_head, version 1.1
%
% Software provided with no warranty or guarantee
%
% Created by: Jon Dudley (jonathan.dudley@cchmc.org)
% Last modified: 2021-07-06
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%
% Function to scan a directory of DICOM files and return a structure
% "sparsehead" containing information from a user defined set of DICOM
% tags. Basically, it's MATLAB's built-in "dicominfo.m" but ~10-100x faster
% with a less verbose output (and far, far less development time behind it,
% so it might not, you know, work all the time).
%
% Example: you just need the ImagePositionPatient information from each of
% 672 DICOM files from a given scan in directory FOLDER
%
% >> tic,
% >> DCMS = dir(fullfile(FOLDER,'*.dcm'));
% >> info = cellfun(@dicominfo,fullfile({DCMS.folder},{DCMS.name}));
% >> toc
% Elapsed time is 59.867770 seconds.
%
% vs. 
%
% >> tic
% >> info = fast_dcm_head(FOLDER,'ImagePositionPatient');
% >> toc
% Elapsed time is 0.435852 seconds.
%
% Most of this elapsed time is overhead associated with opening individual
% DICOM files for reading. In other words, specifying additional tags
% only adds a nominal amount of time. It also means that the performance gains
% are greater when parsing files stored on a device with low latency (e.g. your
% local hard drive) compared to devices with high latency (e.g. network storage 
% locations)
%
% This code was inspired by Brunno Machado de Campos's "alt_dicominfo"
% which you should check out as it might suit your needs better:
% *********************************************************************************** %
% https://www.mathworks.com/matlabcentral/fileexchange/65606-alt_dicominfo-img-bitrange
% *********************************************************************************** %
% fast_dcm_head uses the same premise to quickly scan DICOM headers for
% tags of interest, but leverages MATLAB's default DICOM dictionary to
% allow user to specify/select (almost) any tag, thus it offers more
% flexibility than alt_dicominfo at the cost of requiring more explicit
% input. fast_dcm_head also supports both implicit and explicit value
% representation in DICOM headers.
%
% fast_dcm_head does not (yet?) fully support "enhanced-mode" DICOM files
% (i.e. 1 file per N-Dimensional acquisition vs. the 1 file per 2D slice of
% the "classic" DICOM format). In this case it should work for tags that
% are common across all slices/timepoints/dynamics - like "PatientName" or
% "Manufacturer" - but not for tags with variable information like
% "DiffusionBValue" or "ImagePositionPatient"
%
% Usage:
%   >> head = fast_dcm_head(folder,tags,N);
%       Scans the first "N" DICOM files contained in directory "folder"
%       for DICOM header information specified in "tags". Returns
%       structure "head" with fieldnames matching those in "tags"
%
%       - folder must be a string or character vector specifying a
%       valid filepath containing DICOM files
%       - tags can be a string (e.g. 'DiffusionBValue') or a
%       one-dimensional cell array of strings (e.g.
%       {'PatientBirthName';'PatientAge'}). Entries that do not match
%       tag names listed in MATLAB's default DICOM dictionary will be
%       skipped (and the user will be alerted to this fact)
%       - N is a positive integer. If "N" is greater than the number of
%       DICOM files found in "folder" then fast_dcm_head will scan all
%       DICOM files in "folder"
%
%   >> head = fast_dcm_head(folder,tags);
%       Same as above, but necessarily scans all DICOM files in
%       "folder"
%
%   >> head = fast_dcm_head(folder);
%       Same as above, but will prompt user to select tags from a(n
%       alphabetized) list generated from MATLAB's default DICOM
%       dictionary. NOTE: the more tags that are selected, the longer
%       fast_dcm_head will take to run. At some point, MATLAB's
%       dicominfo will perform faster.
%
%   >> head = fast_dcm_head;
%       Same as above, but prompts user to select a directory
%       containing DICOM files via MATLAB's uigetdir function
%
% Changes in version 1.1
%   > changed way MATLAB's DICOM dictionary was read to support Windows machines

    %------------------------------------------------%
    % Read MATLAB's DICOM dictionary into cell array %
    %------------------------------------------------%
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 4);
    % Specify range and delimiter
    opts.DataLines = [15, Inf];
    opts.Delimiter = "\t";
    % Specify column names and types
    opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4"];
    opts.VariableTypes = ["char", "char", "char", "char"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    % Specify variable properties
    opts = setvaropts(opts, ["VarName1", "VarName2", "VarName3", "VarName4"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["VarName1", "VarName2", "VarName3", "VarName4"], "EmptyFieldRule", "auto");
    % Import the data
    try
        D = readtable(fullfile(matlabroot,'toolbox','images','iptformats',...
            'dicom-dict.txt'), opts);
    catch
        fprintf(2,['Fatal Error:\n'...
            '\tCould not read MATLAB''s default dicom-dict.txt file\n'...
            '\tfast_dcm_head aborted.\n']);
        return
    end
    % Convert to output type
    D = table2cell(D);
    
    % Omit meta-information from options (meta-information always stored as
    % explicit little-endian data and the way this code is currently
    % written can't handle switching between implicit/explicit or
    % big/little endian
    D(strncmp(D(:,1),'(0002',5),:) = [];
    % Omit dictionary entries with ambiguous group numbers from options
    D(contains(D(:,1),'xx'),:) = [];
    % Omit pixel data tags from options
    D(strncmp(D(:,1),'(7FE',4),:) = [];
    
    
    
    % Make properties valid structure fieldnames for output
    D(:,3) = cellfun(@(f) strrep(strrep(f,' ','_'),'-','_'),D(:,3),'Uni',0);

    %-------------------------------------------%
    % Prompt user for input variables if needed %
    %-------------------------------------------%
    if ~exist('dcmdir','var')
        dcmdir = uigetdir(pwd,'Select directory containing DICOM files');
        if ~dcmdir
            return
        end
    end
    dcms = dir(fullfile(dcmdir,'*.dcm'));
    [~,di] = sort(str2double(regexp({dcms.name},'\d+','match','once' )));
    dcms = dcms(di);
    dcms = fullfile(dcmdir,{dcms.name});
    
    if isempty(dcms)
       fprintf(2,['Fatal Error:\n'...
           '\tNo *.dcm files found in ' strrep(dcmdir,'\','\\') '\n'...
           '\tfast_dcm_head aborted.\n']);
       return
    end
    if ~exist('tags','var')
        % Prompt user to select fields to populate
        [onames,o] = sort(D(:,3));  % sort names for list dialog
        D = D(o,:);
        [sel,ok] = listdlg('ListString',onames,'ListSize',[300 500],...
            'PromptString','Select DICOM tags to populate.');
        if ~ok
            fprintf('No selections made\n')
            return
        end
        D = D(sel,:);
    else
        if ~iscell(tags)
            tags = {tags};
        end
        iname = ~ismember(D(:,3),tags);
        itag = ~ismember(D(:,1),tags);
        D(iname & itag,:) = [];
        ntag = ~ismember(tags,D(:,1)) & ~ismember(tags,D(:,3));
        tags(~ntag) = [];
        ni = size(D,1)+1:size(D,1)+sum(ntag);
        for n=1:length(ni)
            if ~contains(tags{n},',')
                continue
            else
                grpelm = strsplit(strrep(strrep(tags{n},'(',''),')',''),',');
                if isnan(str2double(grpelm))
                    continue
                end
            end
            D(ni,1) = tags(n);
            D(ni,2) = {'Unknown'};
            D(ni,3) = {['tag_' strrep(strrep(strrep(tags{n},',','_'),'(',''),')','')]};
        end
        
        if isempty(D)
            fprintf(2,'No matching entries found in dicom dictionary. fast_dcm_head aborted.\n');
            return            
        end
    end
    
    if ~exist('N','var')
        N = length(dcms);
    elseif N>length(dcms)
        N = length(dcms);
    end

    
    sparsehead = struct();
    for f=1:N
        sparsehead(f).Filename = dcms{f};
        % Read first 50k bytes as 8-bit unsigned integers, don't convert
        fa = fopen(dcms{f},'r');
        u8 = fread(fa,5000,'uint8=>uint8');
        
        if f==1
            % Check VR encoding for first file (hopefully it's a safe
            % assumption that all DICOM files of a given scan have the VR
            % encoding
            tag = [2 0 16 0 85 73];
            si = strfind(u8',tag)-1;
            fseek(fa,si+6,'bof');
            l = fread(fa,1,'uint16');
            % Find transfer syntaxUID
            VRenc = deblank(fread(fa,l,'*char')');
            switch VRenc
                case '1.2.840.10008.1.2'
                    % Implicit little-endian VR encoding
                    exp = false;
                    lp = 'uint32';
                    tagL = 4;
                case '1.2.840.10008.1.2.1'
                    exp = true;
                    lp = 'uint16';
                    tagL = 6;
                case '1.2.840.10008.1.2.2'
                    % Explicit Big-endian encoding, rewind and re-read?
                    frewind(fa);
                    u8 = fread(fa,50000,'uint8=>uint8',0,'b');
                    % I don't have any Big-endian encoded DICOM files to test
                    % this out.
                    fprintf(['Warning:\n'...
                       '\tBig-endian encoded data found\n'...
                       '\tfast_dcm_head has not been tested for this condition and may not work properly...\n']);
                otherwise
                    % Some other transfer syntax, let's assume explicit
                    exp = true;
                    lp = 'uint16';
                    tagL = 6;
                    fprintf(['Warning:\n'...
                       '\t\n'...
                       '\tfast_dcm_head has not been tested for this condition and may not work properly...\n']);
            end
        end
        
        %----------------------------------%
        % Loop through selected DICOM tags %
        %----------------------------------%
        for h=1:size(D,1)
            if exp
                tag = [hex2dec(D{h,1}(4:5)) hex2dec(D{h,1}(2:3))...
                    hex2dec(D{h,1}(9:10)) hex2dec(D{h,1}(7:8))...
                    double(D{h,2})];
            else
                tag = [hex2dec(D{h,1}(4:5)) hex2dec(D{h,1}(2:3))...
                    hex2dec(D{h,1}(9:10)) hex2dec(D{h,1}(7:8))];
            end
            si = strfind(u8',tag)-1;
            if isempty(si)  % tag not found
                headval = 'Tag not found';
            else
                si = si(1); % take only first instance
                fseek(fa,si+tagL,'bof');
                l = fread(fa,1,lp);
                if l>512
                    % if length is >512 bytes, probably is the case that
                    % multiple matches to the search tag were found and the
                    % first instance was a "false positive"
                    headval = 'Tag not found';
                else
                    switch D{h,2}
                        case 'FD'
                            headval = fread(fa,l/8,'double');
                        case 'FL'
                            headval = fread(fa,l/4,'single');
                        case 'US'
                            headval = fread(fa,l/2,'uint16');
                        case 'UL'
                            headval = fread(fa,l/4,'uint32');
                        case 'SS'
                            headval = fread(fa,l/2,'int16');
                        case 'SL'
                            headval = fread(fa,l/4,'int32');
                        case {'DS','IS'}
                            headval = str2double(strsplit(deblank(...
                                fread(fa,l,'*char')'),'\'));
                        case {'AE','AS','CS','DA','DT','LO','LT','OD',...
                                'OF','OW','PN','SH','ST','TM','UI','UT'}
                            headval = deblank(fread(fa,l,'*char')');
                        otherwise
                            headval = fread(fa,l/4,'single');
                    end
                end
            end            
            sparsehead(f).(D{h,3}) = headval;
        end
    end
    fclose('all');

end