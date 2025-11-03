function [cmdline,a,b] = dcm2niix(dcmdir,niidir,niiname)
    if strcmp(dcmdir(end),filesep)
        dcmdir = dcmdir(1:end-1);
    end
    [path,folder] = fileparts(dcmdir);
    if ~exist('niidir','var')
        niidir = path;
    end
    if ~isfolder(niidir)
        mkdir(niidir)
    end
    if ~exist('niiname','var')
        niiname = folder;
    end
    dcm2niix = 'C:\Program Files\MRIcroGL\Resources\dcm2niix.exe';
    cmdline = sprintf('"%s" -f "%s" -p y -z y -ba n -o "%s" "%s"',...
        dcm2niix,niiname,niidir,dcmdir);
    [a,b] = system(cmdline);
end