function generate_jobs_replay(subject,traintime,nulldata,lambda)
% This script gets the template and duplicates them, according to the names
% of the 'opts' files in 'batch'
% Made for HOLLY (not myriad)

%% Settings

% Directories on cluster
datadir = '/data/holly-host/jmcfadyen';
scriptdir = '/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts';
batchdir = '/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/batch/';

% Directories on work PC
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';

% Job settings
functionName = 'cluster_replay'; % the function that will be called on the cluster
RAM = '16G'; % RAM allocated to each job

%% Read in jobs

% Read in template file
fid = fopen('template_replay_holly.sh');
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
fname = fullfile(dir_batch,['job_' subject '_t' num2str(traintime) '_n' num2str(nulldata) '_replay.sh']);
jobname = ['s' subject '_t' num2str(traintime) '_n' num2str(nulldata) '_replay'];

for i = 1:length(template) % for each line...

    L = template{i}; % get this line

    % see if there are any tags to be replaced with variables
    L = strrep(L,'[RAM]',[RAM]);
    L = strrep(L,'[DATADIR]',datadir);
    L = strrep(L,'[SCRIPTDIR]',scriptdir);

    L = strrep(L,'[FILENAME]',jobname);
    L = strrep(L,'[FUNCTION]',functionName);
    
    L = strrep(L,'[SUBJECT]',subject);
    L = strrep(L,'[TRAINTIME]',num2str(traintime));
    L = strrep(L,'[NULLDATA]',num2str(nulldata));
    L = strrep(L,'[LAMBDA]',num2str(lambda));

    F{i} = L;

    % write file
    fid = fopen(fname,'w');
    fprintf(fid, '%s\n',F{:}) ;
    fclose(fid) ;
end

end
