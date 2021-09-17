function SPM = update_spm_filepaths(SPM,newdir)

cd(newdir)
load('SPM.mat');

for i = 1:length(SPM.xCon)
    
   obj = SPM.xCon(i).Vcon.private.cdata;
   fname = obj.fname;
   if ~contains(fname,newdir)
       [FILEPATH,NAME,EXT] = fileparts(fname);
       fname = fullfile(newdir,[NAME EXT]);
       obj.fname = fname;
       SPM.xCon(i).Vcon.private.cdata = obj;
   end
   
   obj = SPM.xCon(i).Vspm.private.cdata;
   fname = obj.fname;
   if ~contains(fname,newdir)
       [FILEPATH,NAME,EXT] = fileparts(fname);
       fname = fullfile(newdir,[NAME EXT]);
       obj.fname = fname;
       SPM.xCon(i).Vspm.private.cdata = obj;
   end
   
end

% NOTE: you may need to copy the file from the cluster for SPM.xVol.G
fname = SPM.xVol.G;
if ~contains(fname,newdir)
    [FILEPATH,NAME,EXT] = fileparts(fname);
    SPM.xVol.G = fullfile(newdir,[NAME EXT]);
end

save('SPM.mat','SPM')

end
