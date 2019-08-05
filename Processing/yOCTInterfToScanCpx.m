function [scanCpx,dimensions] = yOCTInterfToScanCpx (varargin)
%This function takes the interferogram loaded from yOCTLoadInterfFromFile
%and converts it to a complex scanCpx datastructure
%
%USAGE:
%       scanCpx = yOCTInterfToBScanCpx (interferogram,dimensions [,param1,value1,...])
%INPUTS
%   - interferogram - as loaded by yOCTLoadInterfFromFile. dimensions
%       should be (lambda,x,...)
%   - dimensions - Dimensions structure as loaded by yOCTLoadInterfFromFile.
%   - Optional Parameters
%       - 'dispersionParameterA', quadradic phase correction units of
%          [nm^2/rad]. Default value is 100
%       - 'band',[start end] - Use a Hann filter to filter out part of the
%          spectrum. Units are [nm]. Default is all spectrum
%       - 'interpMethod', see help yOCTEquispaceInterf for interpetation
%           methods
%		- 'n' - medium refractive index. default: 1.33
%		- 'peakOnly' - if set to true, only returns dimensions update. Default: false
%			dimensions = yOCTInterfToScanCpx (varargin)
%OUTPUT
%   BScanCpx - where lambda dimension is replaced by z
%	dimensions - updated dimesions, adding dimesions for z
%
%Author: Yonatan W (Dec 27, 2017)

%% Hendle Inputs
if (iscell(varargin{1}))
    %the first varible contains a cell with the rest of the varibles, open it
    varargin = varargin{1};
end 

interferogram = varargin{1};
dimensions = varargin{2};

%Optional Parameters
dispersionParameterA = 100; %Default Value
band = [];
interpMethod = []; %Default
n = 1.33;
peakOnly = false;
for i=3:2:length(varargin)
   eval([varargin{i} ' = varargin{i+1};']); %<-TBD - there should be a safer way
end

if (peakOnly)
	interferogram=interferogram(:,1,1,1,1,1); %In peak only mode only process one A Scan
end	

%% Check if interferogram is equispaced. If not, equispace it before processing
lambda = dimensions.lambda.values;
k = 2*pi./(lambda); %Get wave lumber in [1/nm]

if (abs((max(diff(k)) - min(diff(k)))/max(k)) > 1e-10)
    %Not equispaced, equispacing needed
    [interferogram,dimensions] = yOCTEquispaceInterf(interferogram,dimensions,interpMethod);
    
    lambda = dimensions.lambda.values;
    k = 2*pi./(lambda); %Get wave lumber in [1/nm]
end
s = size(interferogram);

%% Filter bands
filter = zeros(size(k));
if ~isempty(band)
    %Provide a warning if band is out of lambda range
    if (band(1) < min(dimensions.lambda.values) || band(2) > max(dimensions.lambda.values))
        warning('Requested band is outside data borders, shrinking band size');
    end
    
    %Band filter to select sub band
    fLambda = linspace(band(1),band(2),length(dimensions.lambda.values));
    fVal = hann(length(fLambda)); 
    filter = interp1(fLambda,fVal,dimensions.lambda.values,'linear',0); %Extrapolation is 0 for values outside the filter
else
    %No band filter, so apply Hann filter on the entire sample
    filter = hann(length(filter));
end

%Normalize filter
filter = filter(:);
filter = filter * (length(filter)/sum(filter)); %Normalization

%% Reshape interferogram for easy parallelization
interf = reshape(interferogram,s(1),[]);

%% Dispersion & Overall filter
if ~exist('dispersionParameterA','var') || isempty(dispersionParameterA)
    disp('Please Enter dispersionParameterA (quadratic correction), recomendations [nm^2/rad]:')
    disp('100 is a good value to start from');
    disp('You can also try running Demo_DispersionCorrection, to figure out the best Value for you');
    return;
end

%Quadratic term only omega
dispersionPhase = -dispersionParameterA .* (k(:)-k(1)).^2; %[rad]
%Technically dispersionPhase = -A*k^2. We added the term -A*(k-k0)^2
%because when doing the ifft, ifft assumes that the first term is DC. which
%in our case is not true. Thus by applying phase=-A*k^2 we introduce a
%multiplicative phase term: A*k0^2 which does not effect the final result
%however, if we run over A the phase term changes and in the fft world it
%translates to translation that move our image up & down. To Avoid it we
%subtract -A*(k-k0)^2
dispersionComp = exp(1i*dispersionPhase);

filterAll = repmat(dispersionComp.*filter,[1 size(interf,2)]);

%% Generate Cpx 
ft = ifft((interf.*filterAll));
scanCpx = ft(1:(size(interf,1)/2),:);

%% Reshape back
scanCpx = reshape(scanCpx,[size(scanCpx,1) s(2:end)]);

%% Update Dimensions
dimensions.z.order = 1;
lambda = mean(dimensions.lambda.values)/1000; %[um]
dlambda = diff(dimensions.lambda.values([1 end]))/1000;%[um]
N = length(dimensions.lambda.values);
zStepSizeAir = 1/2*lambda^2/dlambda; %1/2 factor is because light goes back and forth
zStepSizeMedium = zStepSizeAir/n;
dimensions.z.values = linspace(0,zStepSizeMedium*N/2,N/2); 
dimensions.z.units = 'microns [in medium]';

if (peakOnly)
	scanCpx = dimensions;
end	