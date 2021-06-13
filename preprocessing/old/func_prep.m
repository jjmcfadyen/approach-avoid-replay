function func_prep(datafile,params)

tic

[subject,session,event,in] = pp_osl_cfg(datafile);
disp(['----- FILENAME: ' datafile])
disp(['----- SUBJECT:  ' subject])
disp(['----- SESSION:  ' session])

% update output directory
params.dir_output = fullfile('/Volumes/SOFTWARE BA/processed',subject);
if ~exist(params.dir_output)
    mkdir(params.dir_output);
end
addpath(params.dir_output);

% update the location of the .dat file
filename = fullfile(params.dir_data,datafile);

load(filename);

if isempty(D.fiducials)
    error(['No fiducials for ' datafile])
end

D.path = params.dir_data;
D.data.fname = [filename(1:end-4) '.dat'];

% fix outdated fieldnames
for g = 1:3
    try
        D.sensors.meg.balance.(['G' num2str(g) 'BR']).labelorg = D.sensors.meg.balance.(['G' num2str(g) 'BR']).labelold;
        D.sensors.meg.balance.(['G' num2str(g) 'BR']).chantypeorg = D.sensors.meg.balance.(['G' num2str(g) 'BR']).chantypeold;
        D.sensors.meg.balance.(['G' num2str(g) 'BR']).chanunitorg = D.sensors.meg.balance.(['G' num2str(g) 'BR']).chanunitold;
        D.sensors.meg.balance.(['G' num2str(g) 'BR']) = rmfield(D.sensors.meg.balance.(['G' num2str(g) 'BR']),'labelold');
        D.sensors.meg.balance.(['G' num2str(g) 'BR']) = rmfield(D.sensors.meg.balance.(['G' num2str(g) 'BR']),'chantypeold');
        D.sensors.meg.balance.(['G' num2str(g) 'BR']) = rmfield(D.sensors.meg.balance.(['G' num2str(g) 'BR']),'chanunitold');
    catch
    end
end

% save
save(filename,'D');

% add in response events (forgot to do this in the code that does
% the first half of the preprocessing on Windows PC)
D = spm_eeg_load(filename);
E = D.events;

resp = squeeze(D(find(strcmp(D.chanlabels,params.resp_chan)),:,1));
peaks = find(diff(resp > 0.5) > 0);

cc = length(E);
for r = 1:length(peaks)
    cc = cc + 1;
    E(cc).type = 'response';
    E(cc).value = round(max(resp(peaks(r):peaks(r)+20)));
    E(cc).time = peaks(r)/D.fsample;
    E(cc).duration = 0.2; % in seconds (?)
    E(cc).offset = 0;
end

[~,sortidx] = sort([E.time]);
E = E(sortidx);

D = events(D,1,E);
save(filename);

end