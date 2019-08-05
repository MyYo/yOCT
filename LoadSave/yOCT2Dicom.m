function yOCT2Dicom (data, filepath, dimensions, c)
%This function saves a grayscale version of data to a Dicom stack file.
%Dimensions are (z,x) and each frame is y
%INPUTS
%   filpath - filepath of output Dicom file
%   data - scan data (dimensions are (z,x,y)
%   c - [min, max] of the grayscale

%% Input check
if ~exist('c','var') || isempty(c)
    c = [prctile(data(:),20), prctile(data(:),99.999)]; %min value is at percentile 20 because most of volume has no signals
end

switch dimensions.x.units
    case 'microns'
        SizeX=0.001*(dimensions.x.values(2)-dimensions.x.values(1));
    case 'mm'
        SizeX=dimensions.x.values(2)-dimensions.x.values(1);
    case 'NA'
        SizeX=0;
        disp("error: couldn't assign X pixel size to Dicom file")
end

switch dimensions.y.units
    case 'microns'
        SizeY=0.001*(dimensions.y.values(2)-dimensions.y.values(1));
    case 'mm'
        SizeY=dimensions.y.values(2)-dimensions.y.values(1);
    case 'NA'
        SizeY=0;
        disp("error: couldn't assign Y pixel size to Dicom file")
end

if isfield(dimensions,'z')
    switch dimensions.z.units
        case 'microns'
            SizeZ=0.001*(dimensions.z.values(2)-dimensions.z.values(1));
        case 'mm'
            SizeZ=dimensions.z.values(2)-dimensions.z.values(1);
        case 'NA'
            SizeZ=0;
            disp("error: couldn't assign Z pixel size to Dicom file")
    end
else
    SizeZ=0;
    disp("No Z pixel size was written to Dicom file")
end

if ~isfield(dimensions.aux,'OCTSystem')
    dimensions.aux.OCTSystem = 'Unknown';
end

%% Do we need AWS?
if awsIsAWSPath(filepath)
    %Load Data from AWS
    isAWS = true;
    awsSetCredentials(1); %Use the advanced version as uploading is more challenging
    awsFilePath = filepath;
    awsFilePath = awsModifyPathForCompetability(awsFilePath,true); %We will use this path for AWS CLI
    filepath = [tempname '.dcm'];
else
    isAWS = false;
end

%% Preform writing to the file
for yi=1:size(data,3)
    color(:,:,:) = uint16((squeeze(data(:,:,yi))-c(1))/(c(2)-c(1))*65535);
    color(color>65535) = 65535;
    color(color<0) = 0;
    color4D(:,:,:,yi)=color;
end

dicomwrite(color4D,filepath,'Modality',['OCT: ' dimensions.aux.OCTSystem], ...
    'ObjectType','Secondary Capture Image Storage')

% add resolution parameters
info=dicominfo(filepath);
info.PixelSpacing = [SizeX SizeZ];
info.SpacingBetweenSlices = SizeY;
info.EchoTime = [c(1) c(2)]; %saves previous pixel values
dicomwrite(color4D,filepath, info,'CreateMode','Copy')

%% Upload file to cloud if required
if (isAWS)
    [~,~] = system(['aws s3 cp "' filepath '" "' awsFilePath '"']);
    delete(filepath); %Cleanup
end


