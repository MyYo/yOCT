function awsWriteJSON(json,fp)
%This function write a JSON file from AWS or locally
%json - configuration file
%fp - file path

if (awsIsAWSPath(fp))
    %Load Data from AWS
    isAWS = true;
    fpToSave = 'tmp.json';
else
    isAWS = false;
    fpToSave = fp;
end

%Encode and save
txt = jsonencode(json);
txt = strrep(txt,'"',[newline '"']);
txt = strrep(txt,[newline '":'],'":');
txt = strrep(txt,[':' newline '"'],':"');
txt = strrep(txt,[newline '",'],'",');
txt = strrep(txt,[newline '"}'],['"' newline '}']);
fid = fopen(fpToSave,'w');
fprintf(fid,'%s',txt);
fclose(fid);

if (isAWS)
    %Upload if required
    awsCopyFileFolder(fpToSave,fp); 
    delete(fpToSave);
end