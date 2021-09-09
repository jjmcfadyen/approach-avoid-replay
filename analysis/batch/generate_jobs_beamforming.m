function generate_jobs_beamforming(filename,cluster_type)
% This script gets the template and duplicates them, according to the names
% of the 'opts' files in 'batch'

%% Settings

% Directories on cluster
switch cluster_type
    case 'myriad'
        datadir = '/home/skgtjm6/Scratch';
        scriptdir = '~/Scratch/2020_RiskyReplay/scripts';
    case 'holly'
        datadir = '/data/holly-host/jmcfadyen/2020_RiskyReplay/meg';
        scriptdir = '/data/holly-host/jmcfadyen/2020_RiskyReplay';
end

% Directories on work PC
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';

% Job settings
timechar = '2:00:00'; % max job duration (hours : minutes : seconds)
functionName = 'cluster_beamforming'; % the function that will be called on the cluster
RAM = '32G'; % RAM allocated to each job

%% Read in jobs

% Read in template file
fid = fopen(fullfile(dir_batch,['template_beamforming_' cluster_type '.sh']));
tline = fgetl(fid);
template = {};
while ischar(tline)
    template = [template; tline];
    tline = fgetl(fid);
end
fclose(fid);

% Duplicate the template
F = template;

% what to call copy of template
[FILEPATH,NAME,EXT] = fileparts(filename);
subject = strsplit(NAME,'oat_');
subject = subject{2};

fname = fullfile(dir_batch,['job_' subject '_beamforming.sh']);
jobname = ['s' subject '_beamforming'];

for i = 1:length(template) % for each line...

    L = template{i}; % get this line
    
    % see if there are any tags to be replaced with variables
    L = strrep(L,'[TIME]',timechar);
    L = strrep(L,'[RAM]',RAM);
    L = strrep(L,'[DATADIR]',datadir);
    L = strrep(L,'[SCRIPTDIR]',scriptdir);

    L = strrep(L,'[FILENAME]',jobname);
    L = strrep(L,'[FUNCTION]',functionName);
    L = strrep(L,'[ARGS]',filename);

    F{i} = L;

    % write file
    fid = fopen(fname,'w');
    fprintf(fid, '%s\n',F{:}) ;
    fclose(fid) ;
end

end
