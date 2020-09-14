function yOCTProcessTiledScan(varargin)
%This function Processes Tiled scan, my assumption is that scan size is
%very big, so the processed volume will not be returned directly to Matlab,
%but will be saved directly to disk (or cloud).
%For speed purposes, make sure that input and output folder are either both
%local or both on the cloud. In case, both are on the cloud - run using
%cluster.
%USAGE:
%   yOCTProcessTiledScan(tiledScanInputFolder,outputPath,[params])
%   yOCTProcessTiledScan({parameters})
%INPUTS:
%   - tiledScanInputFolder - where tiled scan is saved. Make sure the
%       ScanInfo.json is present in the folder
%   - outputPath - where products of the processing are saved. can be a
%       a string (path) to tif file or folder (for tif folder). If you
%       input a cell array with both file and folder, will save both
%   - params can be any processing parameters used by
%     yOCTLoadInterfFromFile or yOCTInterfToScanCpx or any of those below
%NAME VALUE INPUTS:
%   Parameter           Default Value   Notes
%Z position stitching:
%   focusSigma                  20      If stitching along Z axis (multiple focus points), what is the size of each focus in z [pixel]
%   focusPositionInImageZpix    NaN     Z position [pix] of focus in each scan (one number)
%   zSetOriginAsFocusOfZDepth0  true    When saving volume how to set z coordinate system. 
%                                       When set to true, will consider the focus position of the volume that was scanned at zDepths=0.
%                                       When set to false, will take the top of the OCT scan for zDepths=0.
%Save some Y planes in a debug folder:
%   yPlanesOutputFolder         ''      If set will save some y planes for debug purpose in that folder
%   howManyYPlanes              3       How many y planes to save (if yPlanesOutput folder is set)
%Other parameters:
%   applyPathLengthCorrection true  Apply path link correction, if probe ini has the information.
%   v                   true        verbose mode      
%
%OUTPUT:
%   No output is returned. Will save mag2db(scan Abs) to outputPath, and
%   debugFolder

%% Parameters
% In case of using focusSigma, how far from focus to go before cutting off
% the signal, avoiding very low values (on log scale)
cuttoffSigma = 3; 

%% Processing of input
p = inputParser;
addRequired(p,'tiledScanInputFolder',@isstr);

% Define the outputs
addRequired(p,'outputPath');

% Z position stitching
addParameter(p,'dispersionQuadraticTerm',79430000,@isnumeric);
addParameter(p,'focusSigma',20,@isnumeric);
addParameter(p,'focusPositionInImageZpix',NaN,@isnumeric);
addParameter(p,'zSetOriginAsFocusOfZDepth0',true);

% Save some Y planes in a debug folder
addParameter(p,'yPlanesOutputFolder','',@isstr);
addParameter(p,'howManyYPlanes',3,@isnumeric);

% Debug
addParameter(p,'v',true,@islogical);

%TODO(yonatan) shift this parameter to ProcessScanFunction
addParameter(p,'applyPathLengthCorrection',true);

p.KeepUnmatched = true;
if (~iscell(varargin{1}))
    parse(p,varargin{:});
else
    parse(p,varargin{1}{:});
end

% Gather unmatched varibles, we will use them as passing inputs
vals = struct2cell(p.Unmatched);
nams = fieldnames(p.Unmatched);
tmp = [nams(:)'; vals(:)']; 
reconstructConfig = tmp(:)';

in = p.Results;
v = in.v;

% Fix input path
tiledScanInputFolder = awsModifyPathForCompetability([fileparts(in.tiledScanInputFolder) '/']);

% Fix output path
outputPath = in.outputPath;
if ischar(outputPath)
    outputPath = {outputPath};
end

% Set credentials
if any(cellfun(@(x)(awsIsAWSPath(x)),outputPath))
    % Any of the output folders is on the cloud
    awsSetCredentials(1);
elseif awsIsAWSPath(in.tiledScanInputFolder)
    % Input folder is on the cloud
    awsSetCredentials();
end

zSetOriginAsFocusOfZDepth0 = in.zSetOriginAsFocusOfZDepth0;
if (zSetOriginAsFocusOfZDepth0 && isnan(in.focusPositionInImageZpix))
    warming('Because no focus position was set, zSetOriginAsFocusOfZDepth0 cannot be "true", changed to "false". See help of yOCTProcessTiledScan function.');
    zSetOriginAsFocusOfZDepth0 = false;
end

%% Load configuration file & set parameters
json = awsReadJSON([tiledScanInputFolder 'ScanInfo.json']);

%Figure out dispersion parameters
if ~isfield(in,'dispersionQuadraticTerm') || isempty(in.dispersionQuadraticTerm)
    if isfield(json.octProbe,'DefaultDispersionParameterA')
        warning('octProbe has dispersionParameterA which is beeing depriciated in favor of dispersionQuadraticTerm. Please adjust probe.ini');
        in.dispersionParameterA = json.octProbe.DefaultDispersionParameterA;
        reconstructConfig = [reconstructConfig {'dispersionParameterA', in.dispersionParameterA}];
    else
        in.dispersionQuadraticTerm = json.octProbe.DefaultDispersionQuadraticTerm;
        reconstructConfig = [reconstructConfig {'dispersionQuadraticTerm', in.dispersionQuadraticTerm}];
        reconstructConfig = [reconstructConfig {'n', json.tissueRefractiveIndex}];
    end
else
    reconstructConfig = [reconstructConfig {'dispersionQuadraticTerm', in.dispersionQuadraticTerm}];
    reconstructConfig = [reconstructConfig {'n', json.tissueRefractiveIndex}];
end

fp = cellfun(@(x)(awsModifyPathForCompetability([tiledScanInputFolder '\' x '\'])),json.octFolders,'UniformOutput',false);
focusPositionInImageZpix = in.focusPositionInImageZpix;
focusSigma = in.focusSigma;
OCTSystem = json.OCTSystem; %Provide OCT system to prevent unesscecary polling of file system
dimOneTile = ...
	yOCTLoadInterfFromFile([fp(1), reconstructConfig, {'OCTSystem',OCTSystem,'peakOnly',true}]);
tmp = zeros(size(dimOneTile.lambda.values(:)));
dimOneTileProcessed = yOCTInterfToScanCpx ([{tmp}, {dimOneTile},{'n'},{json.tissueRefractiveIndex}, reconstructConfig, {'peakOnly'},{true}]);
dimOneTile.z = dimOneTileProcessed.z; %Update only z, not lambda [lambda is changed because of equispacing]
dimOneTile.x.values = linspace(-0.5,0.5,length(dimOneTile.x.values))*json.xRange;
dimOneTile.y.values = linspace(-0.5,0.5,length(dimOneTile.y.values))*json.yRange;
dimOneTile.x.units = 'millimeters';
dimOneTile.y.units = 'millimeters';
dimOneTile = yOCTChangeDimensionsStructureUnits(dimOneTile,'mm');

%TODO(yonatan) move as well
%see if parameter exists
if (in.applyPathLengthCorrection && isfield(json.octProbe,'OpticalPathCorrectionPolynomial'))
    OP_p = json.octProbe.OpticalPathCorrectionPolynomial;
    OP_p = OP_p(:)';
else 
    %No correction
    OP_p = [0 0 0 0 0];
end

%% Create dimensions structure for the entire tiled volume
zDepths = json.zDepths;
xCenters = json.xCenters;
yCenters = json.yCenters;

% Dimensions of one tile (mm)
xOneTile = json.xOffset+json.xRange*linspace(-0.5,0.5,json.nXPixels+1); xOneTile(end) = [];
dx = diff(xOneTile(1:2));
yOneTile = json.yOffset+json.yRange*linspace(-0.5,0.5,json.nYPixels+1); yOneTile(end) = [];
dy = diff(yOneTile(1:2));
zOneTile = dimOneTile.z.values(:)'; %[mm]
dz = diff(zOneTile(1:2));

%Dimensions of the entire stack
xAll = (min(xCenters)+xOneTile(1)):dx:(max(xCenters)+xOneTile(end)+dx);xAll = xAll(:);
yAll = (min(yCenters)+yOneTile(1)):dy:(max(yCenters)+yOneTile(end)+dy);yAll = yAll(:);
zAll = (min(zDepths)+zOneTile(1)):dz:(max(zDepths)+zOneTile(end)+dz);zAll = zAll(:);

% Correct for the case of only one scan
if (length(xCenters) == 1)
    xAll = xAll(1:length(dimOneTileProcessed.x.values));
end
if (length(yCenters) == 1)
    yAll = yAll(1:length(dimOneTileProcessed.y.values));
end

% Remove Z positions that are way out of focus (if we are doing focus processing)
if(~isnan(focusPositionInImageZpix))
    zAll( ...
        ( zAll < min(zDepths) + zOneTile(round(max(focusPositionInImageZpix - cuttoffSigma*in.focusSigma,0))) ) ...
        | ...
        ( zAll > max(zDepths) + zOneTile(round(min(focusPositionInImageZpix + cuttoffSigma*in.focusSigma,length(zOneTile)))) ) ...
        ) = []; 
end

% Dimensions of everything
dimOutput.lambda = dimOneTile.lambda;
dimOutput.z = dimOneTile.z;
dimOutput.z.values = zAll(:)' - ...
    dz*(focusPositionInImageZpix-1)*zSetOriginAsFocusOfZDepth0;
if zSetOriginAsFocusOfZDepth0
    dimOutput.z.origin = 'z=0 is the focus positoin of OCT image when zDepths=0 scan was taken';
else
    dimOutput.z.origin = 'z=0 is the top of OCT image when zDepths=0 scan was taken';
end
dimOutput.x = dimOneTile.x;
dimOutput.x.origin = 'x=0 is OCT scanner origin when xCenters=0 scan was taken';
dimOutput.x.values = xAll(:)';
dimOutput.y = dimOneTile.x;
dimOutput.y.values = yAll(:)';
dimOutput.y.origin = 'y=0 is OCT scanner origin when yCenters=0 scan was taken';
dimOutput.aux = dimOneTile.aux;

%% Save some Y planes in a debug folder if needed
if ~isempty(in.yPlanesOutputFolder) && in.howManyYPlanes > 0
    isSaveSomeYPlanes = true;
    yPlanesOutputFolder = awsModifyPathForCompetability([in.yPlanesOutputFolder '/']);
    
    % Clear folder if it exists
    if awsExist(yPlanesOutputFolder)
        awsRmDir(yPlanesOutputFolder);
    end
    
    %Save some stacks, which ones?
    yToSaveI = round(linspace(1,length(dimOutput.y.values),in.howManyYPlanes));
else
    isSaveSomeYPlanes = false;
    yPlanesOutputFolder = '';
    yToSaveI = [];
end

%% Create indexing reference
%This specifies how to mesh together a tiled scan, each axis seperately

imOutSize = [...
    length(dimOutput.z.values) ...
    length(dimOutput.x.values) ...
    length(dimOutput.y.values)];

%For each Y, what files are being used
yGroupSF = zeros(length(json.yCenters),2); %[start yI, end yI]
yGroupFP =  cell(length(json.yCenters),1); %Which files to use for each group
for i=1:length(yGroupFP)
    yI = abs(dimOutput.y.values-json.yCenters(i)) <= json.yRange/2;
    yGroupSF(i,:) = [find(yI,1,'first') find(yI,1,'last')];
    yGroupFP(i) = {fp(json.gridYcc == json.yCenters(i))};
end

%% Main loop
printStatsEveryyI = max(floor(length(yAll)/20),1);
ticBytes(gcp);
if(v)
    fprintf('%s Stitching ...\n',datestr(datetime)); tt=tic();
end
yOCT2Tif([], outputPath, 'partialFileMode', 1); %Init
parfor yI=1:length(yAll) 
    try
        %Create a container for all data
        stack = zeros(imOutSize(1:2)); %#ok<PFBNS> %z,x,zStach
        totalWeights = zeros(imOutSize(1:2)); %z,x
        
        %Relevant OCT files for this y
        yGroup = find(yGroupSF(:,1) <= yI & yI <= yGroupSF(:,2),1,'first'); %#ok<PFBNS>
        fps = yGroupFP{yGroup}; %#ok<PFBNS>
        fileI = 1;
        
        %What is the y index in the file corresponding to this yI
        %Example, let us assume we took two tiles in the y direction each
        %tile has 4 scans in the y direction. Total number of scans is 8.
        %The first scan in the stitched file is tile #1 scan #1, then scan 2
        %and finally tile #1 scan #3. After this one, we need to move to 
        %tile #2 scan #1. We need to have a varible specifing which scan in
        %the file to grab. Altough we are processing the overall scan #5
        %thus yI=5, the scan we would like to grab is yInFile=1
        yIInFile = yI - yGroupSF(yGroup,1)+1;
        
        %Loop over all x stacks
        for xxI=1:length(xCenters)
            %Loop over depths stacks
            for zzI=1:length(zDepths)
                
                %Frame Name
                fpTxt = fps{fileI};
                fileI = fileI+1;
                
                %Load Frame
                [int1,dim1] = ...
                    yOCTLoadInterfFromFile([{fpTxt}, reconstructConfig, ...
                    {'dimensions', dimOneTile, 'YFramesToProcess', yIInFile}]);
                [scan1,~] = yOCTInterfToScanCpx ([{int1 dim1} reconstructConfig]);
                int1 = []; %#ok<NASGU> %Freeup some memory
                scan1 = abs(scan1);
                for i=length(size(scan1)):-1:3 %Average BScan Averages, A Scan etc
                    scan1 = squeeze(mean(scan1,i));
                end
                
                %Change dimensions to microns units
                dim1 = yOCTChangeDimensionsStructureUnits(dim1,'microns');
                
                %Lens abberation / optical path correction
                %TODO(yonatan) Add this to a seperate function as well so
                %we can call it from any Processing function
                correction = @(x,y)(x*OP_p(1)+y*OP_p(2)+x.^2*OP_p(3)+y.^2*OP_p(4)+x.*y*OP_p(5)); %x,y are in microns
                [xx,zz] = meshgrid(dim1.x.values,dim1.z.values); %um                
                scan1_min = min(scan1(:));
                scan1 = interp2(xx,zz,scan1,xx,zz+correction(xx,dim1.y.values),'nearest');
                scan1_nan = isnan(scan1);
                scan1(scan1<scan1_min) = scan1_min; %Dont let interpolation value go too low
                scan1(scan1_nan) = 0; %interpolated nan values should not contribute to image

                %Filter around the focus
                zI = 1:length(zOneTile); zI = zI(:);
                if ~isnan(focusPositionInImageZpix)
                    factorZ = exp(-(zI-focusPositionInImageZpix).^2/(2*focusSigma)^2) + ...
                        (zI>focusPositionInImageZpix)*exp(-3^2/2);%Under the focus, its possible to not reduce factor as much 
                    factor = repmat(factorZ, [1 size(scan1,2)]);
                else
                    factor = ones(length(zOneTile),length(xOneTile)); %No focus gating
                end
                factor(scan1_nan) = 0; %interpolated nan values should not contribute to image
            
                
                %Figure out what is the x,z position of each pixel in this file
                x = xOneTile+xCenters(xxI);
                z = zOneTile+zDepths(zzI);
                
                %Helps with interpolation problems
                x(1) = x(1) - 1e-10; 
                x(end) = x(end) + 1e-10; 
                z(1) = z(1) - 1e-10; 
                z(end) = z(end) + 1e-10; 
                
                %Add to stack
                [xxAll,zzAll] = meshgrid(xAll,zAll);
                stack = stack + interp2(x,z,scan1.*factor,xxAll,zzAll,'linear',0);
                totalWeights = totalWeights + interp2(x,z,factor,xxAll,zzAll,'linear',0);
                
                %Save Stack, some files for future (debug)
                if (isSaveSomeYPlanes && sum(yI == yToSaveI)>0)
                    
                    tn = [tempname '.tif'];
                    im = mag2db(scan1);
                    if ~isnan(focusPositionInImageZpix)
                        im(focusPositionInImageZpix,1:20:end) = min(im(:)); % Mark focus position on sample
                    end
                    yOCT2Tif(im,tn);
                    awsCopyFile_MW1(tn, ...
                        awsModifyPathForCompetability(sprintf('%s/y%04d_xtile%04d_ztile%04d.tif',yPlanesOutputFolder,yI,xxI,zzI)) ...
                        );
                    delete(tn);
                    
                    if (xxI == length(xCenters) && zzI==length(zDepths))
                        %Save the last weight
                        tn = [tempname '.tif'];
                        yOCT2Tif(totalWeights,tn);
                        awsCopyFile_MW1(tn, ...
                            awsModifyPathForCompetability(sprintf('%s/y%04d_totalWeights.mat',yPlanesOutputFolder,yI)) ...
                            );
                        delete(tn);
                    end
                end
            end
        end
                      
        %Dont allow factor to get too small, it creates an unstable solution
        minFactor1 = exp(-cuttoffSigma^2/2);
        totalWeights(totalWeights<minFactor1) = NaN; 
            
        %Normalization
        stackmean = stack./totalWeights;
        
        % Save
        whereAreMyFiles = yOCT2Tif(mag2db(stackmean), outputPath, ...
            'partialFileMode', 2, 'partialFileModeIndex', yI); 
            % Since we are in partialFileMode 1, whereAreMyFiles will
            % contain the folder that code is working on right now.
        
        % Is it time to print statistics?
        if mod(yI,printStatsEveryyI)==0 && v
            % Stats time!
            ds = fileDatastore(whereAreMyFiles,'ReadFcn',@(x)(x),'FileExtensions','.getmeout','IncludeSubfolders',true); %Count all artifacts
            isFile = cellfun(@(x)(contains(lower(x),'.json')),ds.Files);
            done = sum(isFile);
            fprintf('%s Completed yIs so far: %d/%d (%.1f%%)\n',datestr(datetime),done,length(yAll),100*done/length(yAll));
        end

    catch ME
        fprintf('Error happened in parfor, iteration %d:\n',yI); 
        disp(ME.message);
        for j=1:length(ME.stack) 
            ME.stack(j) 
        end  
        error('Error in parfor');
    end
end %parfor

if (v)
    fprintf('Done stitching, toatl time: %.0f[min]\n',toc(tt)/60);
    tocBytes(gcp)
end

%% Reorganizing files
% Move files outside of their folder
if (v)
    fprintf('%s Finalizing saving tif file ... ',datestr(datetime));
    tt=tic;
end

% Get the main data out
yOCT2Tif([], outputPath, 'metadata', dimOutput, 'partialFileMode', 3);

% Get saved y planes out
if isSaveSomeYPlanes
    if (v)
        fprintf('%s Reorg some y planes ... ',datestr(datetime));
    end
    awsCopyFile_MW2(yPlanesOutputFolder);
end
if (v)
    fprintf('Done! took %.1f[min]\n',toc(tt)/60);
end
