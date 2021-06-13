function parameters = get_parameters(dir_raw)

%% Parameter table

% Create column headings
vN = {'subjectID','schar','rawfile','task','block'};
parameters = [];

% Get subjects from raw data directory
subjects = dir(dir_raw);
subjects = extractfield(subjects(3:end),'name'); % ignore initial '.' and '..'

N = length(subjects);
disp(['Detected ' num2str(N) ' subjects in raw data folder.'])

% Loop through subjects and add raw filenames
for s = 1:N
   
    filelist = extractfield(dir(fullfile(dir_raw,subjects{s},'*.ds')),'name');
    
    nRaw = length(filelist);
    disp(['--- ' num2str(nRaw) ' files for ' subjects{s}])
    
    thisp = array2table(cell(nRaw,length(vN)),'variablenames',vN);
    thisp.subjectID = repmat(str2double(subjects(s)),nRaw,1);
    thisp.schar = repmat(subjects(s),nRaw,1);
    thisp.rawfile = filelist';
    
    thisp.task(1:4) = {'FL'};
    thisp.block(1:4) = num2cell(1:4);
    
    thisp.task(5:end) = {'task'};
    
    if strcmp(subjects{s},'012882')
        thisp.block(5:end) = num2cell([0:7 9:(nRaw-4)]);
    else
        thisp.block(5:end) = num2cell(0:(nRaw-5));
    end
    
    thisp.block = cell2mat(thisp.block);
    
    if isempty(parameters)
        parameters = thisp;
    else
        parameters = [parameters; thisp];
    end
    
end

end