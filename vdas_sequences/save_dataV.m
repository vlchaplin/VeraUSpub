function save_dataV(RcvData)
disp('Saving data..!.');
tic;
% use parameters in workspace
opath = evalin('base','opath');
oname = evalin('base','oname');
fs = evalin('base','fs');
Nacq = evalin('base','numacq');
PRF = evalin('base','PRF');
Resource = evalin('base','Resource');
Trans = evalin('base','Trans');
TW = evalin('base','TW');
TX = evalin('base','TX');
Receive = evalin('base','Receive');
numRcvSamples = evalin('base','rowsPerAcq');

% redefine transmit frequency in Hz
f0 = Trans.frequency*1e6;

% set up channel data matrix
%numRcvSamples = ceil((Receive(1).endDepth-Receive(1).startDepth)/f0*fs*2); % number of rows per acq
%numRcvSamples = 
tmpDat = RcvData;
tmpDat = tmpDat(1:Nacq*numRcvSamples,:,:); % cut out non-acq data
chandat = reshape(tmpDat,[numRcvSamples,Nacq,size(tmpDat,2),Resource.RcvBuffer(1).numFrames]);
clear tmpDat

% redefine dimensions and reshape data accordingly
%chandat = permute(chandat,[1 3 2 4]);
%chanSz = size(chandat);
%chandat = reshape(chandat,[chanSz(1), chanSz(2), prod(chanSz(3:4))]);

% store output parameters in structure
params.numRcvChannels = Resource.Parameters.numRcvChannels;
params.numRcvSamples = numRcvSamples;
params.c = Resource.Parameters.speedOfSound; % m/s
params.prf = PRF; % Hz
params.pitch = Trans.spacingMm/1000; % m
params.fs = fs; % Hz
params.f0 = f0; % Hz
params.numacq = Nacq;
params.numframes = Resource.RcvBuffer(1).numFrames;
lensCorr = Trans.lensCorrection/1000/params.c*f0; % convert to wavelengths
params.t0= round((Receive(1).startDepth+lensCorr.*2+TW(1).peak+max(TX(1).Delay))*Receive(1).samplesPerWave);

% save the data
if ~exist(opath,'dir'); mkdir(opath); end
save(fullfile(opath,[oname,'_params']),'params');
fid = fopen(fullfile(opath,[oname '.bin']),'w');
fwrite(fid,chandat,'int16');
fclose(fid);

% display when done saving
tsave = toc;
disp(['Data saved. Elapsed time is ' num2str(tsave)]);
close all;

end