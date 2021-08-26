function generate_jobs_timefrequency(filename)
% This script gets the template and duplicates them, according to the names
% of the 'opts' files in 'batch'

%% Settings

% Directories on cluster
datadir = '/home/skgtjm6/Scratch';
scriptdir = '~/Scratch/2020_RiskyReplay/scripts';

% Directories on work PC
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';

% Job settings
timechar = '3:30:00'; % max job duration (hours : minutes : seconds)
functionName = 'cluster_timefrequency'; % the function that will be called on the cluster
RAM = '8G'; % RAM allocated to each job

%% Read in jobs

% Read in template file
fid = fopen(fullfile(dir_batch,'template_timefrequency.sh'));
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
filesplit = strsplit(filename, '_');
subject = filesplit{2}(end-5:end);
waveletwidth = filesplit{3};

fname = fullfile(dir_batch,['job_' subject '_' waveletwidth '_timefrequency.sh']);
jobname = ['s' subject '_' waveletwidth '_timefrequency'];

for i = 1:length(template) % for each line...

    L = template{i}; % get this line
    
    % see if there are any tags to be replaced with variables
    L = strrep(L,'[TIME]',timechar);
    L = strrep(L,'[RAM]',['mem=' RAM]);
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
