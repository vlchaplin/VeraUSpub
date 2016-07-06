function save_bmode_rcvbuff(RcvData)
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

num_final_bmode_acs = evalin('base','num_final_bmode_acs');

tmpDat = RcvData;
tmpDat = tmpDat(1:num_final_bmode_acs*numRcvSamples,:,:); % cut out non-acq data
chandat = reshape(tmpDat,[numRcvSamples, num_final_bmode_acs, size(tmpDat,2), Resource.RcvBuffer(2).numFrames]);
clear tmpDat

% save bmode data
if ~exist(opath,'dir'); mkdir(opath); end
%save(fullfile(opath,[oname,'_paramss']),'params');
fid = fopen(fullfile(opath,[oname '.bin']),'a');
fwrite(fid,chandat,'int16');
fclose(fid);

% display when done saving
tsave = toc;
disp(['Data saved. Elapsed time is ' num2str(tsave)]);
close all;

end