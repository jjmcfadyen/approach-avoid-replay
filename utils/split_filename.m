function [subject, task, run] = split_filename(filename)

filesplit = strsplit(filename, '_');
for i = 1:length(filesplit)
    if ~isnan(str2double(filesplit{i}))
        break
    end
end

subject = filesplit{i};
task = filesplit{i+1};
run = strsplit(filesplit{i+2}(2:end),'.mat');
run = str2double(run{1});

end