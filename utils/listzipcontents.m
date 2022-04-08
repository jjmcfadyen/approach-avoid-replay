function filelist = listzipcontents(zipFilename)
% Create a Java file of the ZIP filename.
zipJavaFile  = java.io.File(zipFilename);
% Create a Java ZipFile and validate it.
zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);
% Extract the entries from the ZipFile.
entries = zipFile.getEntries;
% Initialize the file list.
filelist={};
% Loop through the entries and add to the file list.
while entries.hasMoreElements
    filelist = cat(1,filelist,char(entries.nextElement));
end
end