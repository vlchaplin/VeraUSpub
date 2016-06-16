dir='C:\Users\vchaplin\SkyDrive\Notes,Maths\recon\screwdriver\';
dir='C:\Users\Vandiver\Data\Verasonics\sonalleve_0607\';
%dir='C:\Users\Vandiver\Data\Verasonics\h101_agar_0527\';

oname = '3foc-agar-1000av';


% dir='C:\Users\Vandiver\Data\Verasonics\deeper target\';
% oname = 'bubble_burst_center_400mV_5cyc_deep';
% 
% dir='C:\Users\Vandiver\Data\Verasonics\increased depth\';
% oname='bubble_burst_center_400mV_5cyc_more_depth';

fid = fopen([dir oname '.bin']);
chandat2 = fread(fid,Inf,'int16');
fclose(fid);

load([dir oname '_params']);

chandat2 = reshape(chandat2,[params.numRcvSamples params.numacq params.numRcvChannels params.numframes]);

% redefine dimensions and reshape data accordingly
%chandat2 = permute(chandat2,[1 3 2 4]);
%chanSz = size(chandat2);
%chandat2 = reshape(chandat2,[chanSz(1), chanSz(2), prod(chanSz(3:4))]);
% pre-allocate rf data
rf_data = zeros([params.numRcvSamples 128 params.numacq*params.numframes]);

for m=1:params.numframes
    m
    for k=1:params.numacq
        my_i = (m-1)*params.numacq+k;
        rf_data(:,:,my_i) = squeeze(chandat2(:,k,:,m));
    end
end;
