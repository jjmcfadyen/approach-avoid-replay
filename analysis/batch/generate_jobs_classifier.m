function generate_jobs_classifier(subject,traintime,nulldata,cluster_type)
% This script gets the template and duplicates them, according to the names
% of the 'opts' files in 'batch'

jobarray = false;
if ischar(traintime)
    if any(contains(traintime,'-')) && strcmp(cluster_type,'holly')
        jobarray = true;
    end
end

%% Settings

% Directories on cluster
switch cluster_type
    case 'myriad'
        datadir = '/home/skgtjm6/Scratch';
        scriptdir = '~/Scratch/2020_RiskyReplay/scripts';
        batchdir = '/lustre/scratch/scratch/skgtjm6/2020_RiskyReplay/scripts/';
    case 'holly'
        datadir = '/data/holly-host/jmcfadyen';
        scriptdir = '/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts';
        batchdir = '/data/holly-host/jmcfadyen/2020_RiskyReplay/scripts/batch/';
end

% Directories on work PC
dir_batch = 'D:\2020_RiskyReplay\approach-avoid-replay\analysis\batch';

% Job settings
timechar = '2:30:00'; % max job duration (hours : minutes : seconds)
functionName = 'crossvalidate_classifier'; % the function that will be called on the cluster
RAM = '16G'; % RAM allocated to each job

%% Read in jobs

% Read in template file
fid = fopen(['template_classifier_' cluster_type '.sh']);
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
if jobarray && strcmp(cluster_type,'holly')
    fname = fullfile(dir_batch,subject,['job_' subject '_idx' num2str(traintime) '_n' num2str(nulldata) '_classifier.sh']);
    jobname = ['s' subject '_idx' num2str(traintime) '_n' num2str(nulldata) '_classifier'];
    argname = ['[],''' batchdir ''',''' subject ''',$SGE_TASK_ID,' num2str(nulldata)];
else
    fname = fullfile(dir_batch,subject,['job_' subject '_t' num2str(traintime) '_n' num2str(nulldata) '_classifier.sh']);
    jobname = ['s' subject '_t' num2str(traintime) '_n' num2str(nulldata) '_classifier'];
    argname = ['''' batchdir,subject,'/data_',subject '_t' num2str(traintime) '_n' num2str(nulldata) '.mat'''];
end

for i = 1:length(template) % for each line...

    L = template{i}; % get this line

    % see if there are any tags to be replaced with variables
    L = strrep(L,'[TIME]',timechar);
    L = strrep(L,'[RAM]',RAM);
    L = strrep(L,'[DATADIR]',datadir);
    L = strrep(L,'[SCRIPTDIR]',scriptdir);

    L = strrep(L,'[FILENAME]',jobname);
    L = strrep(L,'[FUNCTION]',functionName);
    
    L = strrep(L,'[ARGS]',argname);

    if jobarray
        L = strrep(L,'[JOBIDX]',traintime);
    else
        L = strrep(L,'[JOBIDX]',num2str(1));
    end

    F{i} = L;

    % write file
    fid = fopen(fname,'w');
    fprintf(fid, '%s\n',F{:}) ;
    fclose(fid) ;
end

% copyfile(fullfile(dir_batch,'jobController.sh'),fullfile(dir_batch,subject,'jobController.sh'));

end
