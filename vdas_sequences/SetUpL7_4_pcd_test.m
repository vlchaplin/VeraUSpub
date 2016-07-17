clear all
VSXName='C:\Users\verasonics\Desktop\Vantage-3.0.4/MatFiles/L7_4_pcd_test';

%% Define variables
opath = 'C:\Users\verasonics\Desktop\Vantage-3.0.4\Vandiver\'; % path to output files
oname = 'pcd_test';
numacq = 20; % number of acquisitions per frame... also the number of triggers recorded per frame
numFrames = 50; % number of frames
%PRF = 1000; % pulses per second
PRF=-1;
Nhalfcycles = 2; 

ext_trig = 1; %wait for external trigger

%% Trigger setup 
%!! NOT YET IMPLEMENTED... ONLY ONE ACQ PER TRIG NOW !!

%for ext_trig giving one image per trigger
acquisitions_per_trig = 1;

%time between acq in microseconds. If only one acq per exttrig, this time
% is determined solely by the external trigger rate
if acquisitions_per_trig == 1
    us_between_acq=0;
else
    us_between_acq = 1000;
end

%% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.simulateMode = 0;

Resource.VDAS.dmaTimeout=20e3; %increase the dma timeout to a large value. This is necesary because when waiting for an input trigger the timeout error was occuring.
%% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.

samplesPerWave=4;
wavesPerSample=1.0/samplesPerWave;

fs = samplesPerWave*Trans.frequency*1e6; % found in manual

%% Depth and scan size
wavelength_cm = 1e2*Resource.Parameters.speedOfSound/(Trans.frequency*1e6);

%make endDepth slightly larger than desired to account for longer path
%length of edge channels

desiredImageDepth_cm = 6.0;

startDepth_cm = 0;
endDepth_cm = sqrt(desiredImageDepth_cm^2 + (Trans.numelements*Trans.spacingMm*0.1)^2);

startDepth_wls = startDepth_cm / wavelength_cm;
endDepth_wls = endDepth_cm / wavelength_cm; %end depth in wavelengths

wlsPer128samples = 128/(samplesPerWave*2); % wavelengths in 128 samples for 4 samplesPerWave

%Adjust so buffer is a multiple of 128 samples.  Have to calculate buffer size via
%wavelengths for the scan format object.
%
% E.g., At 4 samples per wave, we want endDepth_wls-startDepth_wls to be a multiple of 16 wavelengths for optimal performance, as
%described in the Scan Fromat section of the Vantage manual. 
%Adjust endDepth upwards to get to this multiple.
qFact=(endDepth_wls-startDepth_wls)/wlsPer128samples;
multiplesOf128samples=ceil(qFact);

endDepth_wls = multiplesOf128samples*wlsPer128samples;
endDepth_cm = endDepth_wls*wavelength_cm;

%number of rows (samples) per acquisition. The above rounding will result
%in this being an integer number, fix() just converts to int data type
rowsPerAcq = fix( (endDepth_wls - startDepth_wls)*wlsPer128samples/2 );

%% Specify SFormat structure array.
SFormat.transducer = 'L7-4'; % 256 element linear array with 0.96 lambda spacing
SFormat.scanFormat = 'RLIN';     % rectangular linear array scan
SFormat.radius = 0;              % ROC for curved lin. or dist. to virt. apex
SFormat.theta = 0;
SFormat.numRays = 1;      % no. of Rays (1 for Flat Focus)
SFormat.FirstRayLoc = [0,0,0];   % x,y,z
SFormat.rayDelta = 128*Trans.spacing;  % spacing in radians(sector) or dist. between rays (wvlnghts)
SFormat.startDepth = startDepth_wls;   % Acquisition depth in wavelengths
SFormat.endDepth = endDepth_wls;

%% Specify PData structure array.
PData.sFormat = 1;      % use first SFormat structure.
PData.pdeltaX = Trans.spacing;
PData.pdeltaZ = 0.5;
PData.Size(1) = ceil((SFormat.endDepth-SFormat.startDepth)/PData.pdeltaZ); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.pdeltaX);
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,SFormat.startDepth]; % x,y,z of upper lft crnr.

%% Specify Resources.                                    
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = rowsPerAcq*numacq;%4096*33; % fits 50+ acquisions per frame
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numFrames;%300;%60;  % try 60 for 3 seconds worth of data % 50 frames used for RF cineloop.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).rowsPerFrame = rowsPerAcq;%1024; % this is for greatest depth
Resource.InterBuffer(1).colsPerFrame = PData.Size(2);
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).rowsPerFrame = rowsPerAcq;%1024;
Resource.ImageBuffer(1).colsPerFrame = PData.Size(2);
Resource.ImageBuffer(1).numFrames = numFrames;
Resource.DisplayWindow(1).Title = 'Flash';
Resource.DisplayWindow(1).pdelta = 0.45;
Resource.DisplayWindow(1).Position = [250,250, ...    % upper left corner position
    ceil(PData.Size(2)*PData.pdeltaX/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData.Size(1)*PData.pdeltaZ/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),PData.Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Colormap = gray(256);

%Make an additional recieve buffer resource to hold the initial b-mode acquisition
% Resource.RcvBuffer(2).datatype = 'int16';
% Resource.RcvBuffer(2).rowsPerFrame = rowsPerAcq*1;%4096*33; % fits 50+ acquisions per frame
% Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
% Resource.RcvBuffer(2).numFrames = 1;

%% Specify Transmit waveform structure. 
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,Nhalfcycles,1];   % A, B, C, D  for 8.1 MHz transmit frequency

%% Specify TX structure array.  
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', ones(1,Resource.Parameters.numTransmit), ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, numacq);
%TX(2).aperture = 1;  % Use the tx aperture that starts at element 65.
%TX(3).aperture = 1; % Use the tx aperture that starts at element 129.

%% Specify TGC Waveform structure.
TGC.CntrlPts = [234,368,514,609,750,856,1023,1023];
TGC.rangeMax = SFormat.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Receive structure arrays -

Receive = repmat(struct('Apod', zeros(1,128), ...
                        'startDepth', SFormat.startDepth, ...
                        'endDepth', SFormat.endDepth, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'samplesPerWave', 4, ...
                        'mode', 0, ...
                        'callMediaFunc', 1),1,numacq*Resource.RcvBuffer(1).numFrames + 1); 
                    
                    
% - Set event specific Receive attributes.
nn=0;

for i = 1:Resource.RcvBuffer(1).numFrames  % for each frame
    for j=1:numacq; % for each acquisition
        nn=nn+1;
        Receive(nn).Apod(1:128) = 1.0;
        Receive(nn).framenum = i;
        Receive(nn).acqNum = j; %mod(nn,20);

    end
end

%nn=nn+1;
%Recieve struct for the initial b-mode acquisition
% Receive(nn).bufnum=2;
% Receive(nn).Apod(1:128) = 1.0;
% Receive(nn).framenum = 1;
% Receive(nn).acqNum = 1;

bmode_save_buff=nn;

%%  Recon structure array
% - We need one Recon structure which will be used for each frame. 
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',1);

% Define ReconInfo structures.
% We need numacq ReconInfo structures for numacq transmits per frame.
ReconInfo = repmat(struct('mode', 0, ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, 1);

%% Specify Process structure array.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'norm',1,...        % normalization method(1 means fixed)
                         'pgain',5.0,...            % pgain is image processing gain
                         'persistMethod','none',...
                         'persistLevel',0,...
                         'interp',1,...      % method of interpolation (1=4pt interp)
                         'compression',0.5,...      % X^0.5 normalized to output word size
                         'reject',2,...
                         'mappingMode','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                     
% This one uses out external processing function
Process(2).classname = 'External';
Process(2).method = 'display_done'; 
Process(2).Parameters = {'srcbuffer','none',...
                         'srcbufnum',1,...
                         'srcframenum',0,...
                         'dstbuffer','none'};

Process(3).classname = 'External';
Process(3).method = 'save_dataV';
Process(3).Parameters = {'srcbuffer','receive',...
                         'srcbufnum',1,...
                         'srcframenum',0,...
                         'dstbuffer','none'};
                     
% Process(4).classname = 'External';
% Process(4).method = 'save_bmode_rcvbuff';
% Process(4).Parameters = {'srcbuffer','receive',...
%                          'srcbufnum',2,...
%                          'srcframenum',0,...
%                          'dstbuffer','none'};
                     
                     
%% Specify SeqControl structure arrays.
si=1;
SeqControl(si).command = 'jump'; % jump back to start.
SeqControl(si).argument = 1;
si=si+1;

ACQWAIT_seqidx=si; si=si+1;
SeqControl(ACQWAIT_seqidx).command = 'timeToNextAcq';  % time between real-time acquisitions
SeqControl(ACQWAIT_seqidx).argument = 1/50*1e6; % usec

RETURN_seqidx=si; si=si+1;
SeqControl(RETURN_seqidx).command = 'returnToMatlab';

NOOP_seqidx=si; si=si+1;
SeqControl(NOOP_seqidx).command = 'noop';
SeqControl(NOOP_seqidx).argument = 500000; % wait in usec.

SYNC_seqidx=si; si=si+1;
SeqControl(SYNC_seqidx).command = 'sync';
SeqControl(SYNC_seqidx).argument = 5e6; % wait in usec.

trigger=0;
if trigger
    TRIGOUT_seqidx=si;  si=si+1;
    SeqControl(TRIGOUT_seqidx).command = 'triggerOut';
    nsc = si; % nsc is count of SeqControl objects
elseif ext_trig
 
    TRIGIN_seqidx=si;  si=si+1;
    SeqControl(TRIGIN_seqidx).command = 'triggerIn';
    SeqControl(TRIGIN_seqidx).condition = 'Trigger_1_Rising'; % input BNC #1 trigger
    SeqControl(TRIGIN_seqidx).argument = 0; % 500 msec timeout delay
    % (Timeout range is 1:255 in 250 msec steps; 0 means timeout disabled)
    
    nsc = si;
else
    nsc = si;
end

%% Set up B-mode display loop at the beginning
n = 1; % n is count of Events

% real-time acquisition!
Event(n).info = 'Flash transmit';
Event(n).tx = 1;         % use 1st TX structure.
Event(n).rcv = 1;%numacq*i-(numacq-1);        % use 1st Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = [2,nsc]; % time between acqs.
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
n = n+1;

Event(n).info = 'reconstruct'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 1;      % reconstruction
Event(n).process = 0;    % process
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'process (Display B-Mode) and return to Matlab'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % reconstruction
Event(n).process = 1;    % process
Event(n).seqControl = RETURN_seqidx;
n = n+1;

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 1; % jump command
n = n+1;
lastBmodeEvent = n;

%% Set up acquisition events
%          
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:numacq 
        
        Event(n).info = 'Passive recieve';
        Event(n).tx = 1;%j;         % else no transmit
        Event(n).rcv = numacq*i-(numacq-j);   
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing

        if (i==1 && j==1)
            %For the initial triggered event, wait a long time for trigger, so HIFU
            %system has time to start

            
            
            Event(n).seqControl = [nsc, SYNC_seqidx];
            SeqControl(nsc).command = 'triggerIn';
            SeqControl(nsc).condition = 'Trigger_1_Rising'; % input BNC #1 trigger
            SeqControl(nsc).argument = 0;  % (Timeout range is 1:255 in 250 msec steps; 0 means timeout disabled)
            nsc=nsc+1;
            
%             SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
%             nsc = nsc+1;
            
        else
            %one image acq for each trigger
            Event(n).seqControl = [TRIGIN_seqidx, SYNC_seqidx];
%             SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
%             nsc = nsc+1;
        end
        
        n = n+1;
             
        
        %% testing 
%         Event(n).info = 'Transfer frame to host.';
%         Event(n).tx = 0;        % no TX
%         Event(n).rcv = 0;       % no Rcv
%         Event(n).recon = 0;     % no Recon
%         Event(n).process = 0; 
%         Event(n).seqControl = nsc; 
       
%         
%         Event(n).info = 'Wait for transfer to complete'; 
%         Event(n).tx = 0;         % no transmit
%         Event(n).rcv = 0;        % no rcv
%         Event(n).recon = 0;      % no reconstruction
%         Event(n).process = 0;    % no processing
%         Event(n).seqControl = nsc; % wait for transfer complete
%         SeqControl(nsc).command = 'waitForTransferComplete';
%         SeqControl(nsc).argument = nsc - 1;
%         nsc = nsc + 1;
%         n = n+1;
%                 
%         Event(n).info = 'Mark transfer processed'; 
%         Event(n).tx = 0;         % no transmit
%         Event(n).rcv = 0;        % no rcv
%         Event(n).recon = 0;      % no reconstruction
%         Event(n).process = 0;    % no processing
%         Event(n).seqControl = nsc; % wait for transfer complete
%            SeqControl(nsc).command = 'markTransferProcessed';
%            SeqControl(nsc).argument = nsc - 2;
%            nsc = nsc + 1;
%         n = n+1;
        
%         Event(n).info = 'reconstruct'; 
%         Event(n).tx = 0;         % no transmit
%         Event(n).rcv = 0;        % no rcv
%         Event(n).recon = 1;      % reconstruction
%         Event(n).process = 0;    % process
%         Event(n).seqControl = 0;
%         n = n+1;
% 
%         Event(n).info = 'process (Display B-Mode) and return to Matlab'; 
%         Event(n).tx = 0;         % no transmit
%         Event(n).rcv = 0;        % no rcv
%         Event(n).recon = 0;      % reconstruction
%         Event(n).process = 1;    % process
%         Event(n).seqControl = RETURN_seqidx;
%         n = n+1;
        
        %%
    end
%     
    Event(n).info = 'Transfer frame to host.';
    Event(n).tx = 0;        % no TX
    Event(n).rcv = 0;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0; 
    Event(n).seqControl = nsc; 
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'Wait for transfer complete & return to catch matlab events';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 0;    % process
    Event(n).seqControl = [nsc,nsc+1, RETURN_seqidx]; % wait for data to be transferred, mark processed
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'markTransferProcessed';
    SeqControl(nsc+1).argument = nsc-1;
    nsc = nsc+2;
    n = n+1;  


%         Event(n).info = 'process (Display B-Mode) and return to Matlab'; 
%         Event(n).tx = 0;         % no transmit
%         Event(n).rcv = 0;        % no rcv
%         Event(n).recon = 1;      % reconstruction
%         Event(n).process = 1;    % process
%         Event(n).seqControl = [SYNC_seqidx, RETURN_seqidx];
%         n = n+1;

end
% 
% 
% Event(n).info = 'Insert wait before flash b-mode';
% Event(n).tx = 0;%j;         %do T/R to get a B-mode image
% Event(n).rcv = 0;    
% Event(n).recon = 0;      % no reconstruction.
% Event(n).process = 0; 
% Event(n).seqControl = NOOP_seqidx;
% n=n+1;
% 
% %very last acquisition is a B-mode
% %for very first event, transmit plane wave, recieve and copy to
%             %host
% Event(n).info = 'Flash transmit';
% Event(n).tx = 1;%j;         %do T/R to get a B-mode image
% Event(n).rcv = bmode_save_buff;    
% Event(n).recon = 0;      % no reconstruction.
% Event(n).process = 0;    % no processing
% Event(n).seqControl = nsc;
% SeqControl(nsc).command = 'transferToHost';
% nsc=nsc+1;
% n=n+1;
% 
% Event(n).info = 'Wait for transfer complete';
% Event(n).tx = 0;         % no transmit
% Event(n).rcv = 0;        % no rcv
% Event(n).recon = 0;      % reconstruction
% Event(n).process = 0;    % process
% Event(n).seqControl = [nsc,nsc+1]; % wait for data to be transferred, mark processed
% SeqControl(nsc).command = 'waitForTransferComplete';
% SeqControl(nsc).argument = nsc-1;
% SeqControl(nsc+1).command = 'markTransferProcessed';
% SeqControl(nsc+1).argument = nsc-1;
% nsc = nsc+2;
% n = n+1;
% 
% 
Event(n).info = 'save data';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 3;
Event(n).seqControl = NOOP_seqidx;
n = n+1;
% 
% Event(n).info = 'save b-mode data';
% Event(n).tx = 0;         % no transmit
% Event(n).rcv = 0;        % no rcv
% Event(n).recon = 0;      % reconstruction
% Event(n).process = 4;    % process
% Event(n).seqControl = NOOP_seqidx;
% n = n+1;
% 
% Event(n).info = 'Display when sequence is done';
% Event(n).tx = 0;         % no transmit
% Event(n).rcv = 0;        % no rcv
% Event(n).recon = 0;      % no reconstruction
% Event(n).process = 2;    % process
% Event(n).seqControl = NOOP_seqidx; % NOOP
% n = n+1;

Event(n).info = 'Back to Matlab';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % no process
Event(n).seqControl =RETURN_seqidx; % Back to Matlab
n = n+1;
    
    
%% User specified UI Control Elements
% - Sensitivity Cutoff
sensx = 170;
sensy = 420;
UI(1).Control = {'Style','text',...        % popupmenu gives list of choices
                 'String','Sens. Cutoff',...
                 'Position',[sensx+10,sensy,100,20],... % position on UI
                 'FontName','Arial','FontWeight','bold','FontSize',12,...
                 'BackgroundColor',[0.8,0.8,0.8]};
UI(2).Control = {'Style','slider',...        % popupmenu gives list of choices
                 'Position',[sensx,sensy-30,120,30],... % position on UI
                 'Max',1.0,'Min',0,'Value',Recon(1).senscutoff,...
                 'SliderStep',[0.025 0.1],...
                 'Callback',{@sensCutoffCallback}};
UI(2).Callback = {'sensCutoffCallback.m',...
                 'function sensCutoffCallback(hObject,eventdata)',...
                 ' ',...
                 'sens = get(hObject,''Value'');',...
                 'ReconL = evalin(''base'', ''Recon'');',...
                 'for i = 1:size(ReconL,2)',...
                 '    ReconL(i).senscutoff = sens;',...
                 'end',...
                 'assignin(''base'',''Recon'',ReconL);',...
                 '% Set Control.Command to re-initialize Recon structure.',...
                 'Control = evalin(''base'',''Control'');',...
                 'Control.Command = ''update&Run'';',...
                 'Control.Parameters = {''Recon''};',...
                 'assignin(''base'',''Control'', Control);',...
                 '% Set the new cutoff value in the text field.',...
                 'h = findobj(''tag'',''sensCutoffValue'');',...
                 'set(h,''String'',num2str(sens,''%1.3f''));',...
                 'return'};
UI(3).Control = {'Style','edit','String',num2str(Recon(1).senscutoff,'%1.3f'), ...  % text
                 'Position',[sensx+20,sensy-40,60,22], ...   % position on UI
                 'tag','sensCutoffValue', ...
                 'BackgroundColor',[0.9,0.9,0.9]}; 
     
 % -- Enable DisplayWindow's WindowButtonDown callback function for switching acquisition loops.
UI(4).Statement = 'set(Resource.DisplayWindow(1).figureHandle,''WindowButtonDownFcn'',@wbdCallback);';
UI(4).Callback = {'wbdCallback.m',...
    'function wbdCallback(hObject,eventdata)',...
    ' ',...
    'persistent init wbFig wbAxesl',...
    'if isempty(init)',...
    '   wbFig = evalin(''base'',''Resource.DisplayWindow(1).figureHandle'');',...
    '   wbAxes = get(wbFig,''CurrentAxes'');',...
    '   init = 1;',...
    'end',...
    '% if left mouse button ...',...
    'if strcmp(get(hObject,''SelectionType''),''normal'')',...
    '   % set startEvent',...
    '   lastBmodeEvent = evalin(''base'',''lastBmodeEvent'')',...
    '   Control = evalin(''base'',''Control'');',...
    '   Control(1).Command = ''set&Run'';',...
    '   Control(1).Parameters = {''Parameters'',1,''startEvent'',lastBmodeEvent};',...
    '   evalin(''base'',''Resource.Parameters.startEvent = lastBmodeEvent;'');',...
    '   assignin(''base'',''Control'', Control);',...
    'end',...
    'return'};       

clear i j n sensx sensy

% Save all the structures to a .mat file.
display(['filename =''' VSXName ''';VSX'])
save(VSXName);
    
    