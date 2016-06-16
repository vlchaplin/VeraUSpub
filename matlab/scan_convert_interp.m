

%% setup variables
c=1520; % [m/s] assume tissue
Trans.spacing=params.pitch;
Trans.frequency=params.fs/4*1e-6;
lambda=c/(Trans.frequency*1e6);

[depth_m sensor_n, nBuff] = size(rf_data);

sampleZero=0; %should be zero if the physical image begins at sample 1

Fs = 4*Trans.frequency*1e6;
dt = 1/Fs;

%lambda=c/1.1e6;

dx = lambda/2;
dz = lambda/2;
%dz =c*dt;
%dx=Trans.spacing;
Nx = round(sensor_n*Trans.spacing/dx);
%Nz = 450;
Nz = round( (depth_m-sampleZero)/2.0*c*dt/dz);

%Nz = round( (6e-2)/dz);

xaxis = (-(Nx/2*dx):dx:(Nx/2*dx-dx));
zaxis = (0:dz:(Nz-1)*dz);
dsensor = Trans.spacing;

sensors_x = ((0:sensor_n-1) + 0.5 - sensor_n/2.0)*dsensor;


[ndZ, ndX, ndsensX] = ndgrid( zaxis, xaxis, sensors_x);

distances_ndgrid = sqrt( (ndX-ndsensX).^2 + (ndZ + sampleZero*c*dt).^2 );
delays_ndgrid = distances_ndgrid/c;

delayedsampling_ndgrid = round( delays_ndgrid/dt ) + 1;

inbounds3 =  delayedsampling_ndgrid <= depth_m;

ii = find(inbounds3);
[zii,xii,sii]=ind2sub([Nz Nx sensor_n],ii);
rf_sample_indices_flat = sub2ind([depth_m sensor_n], delayedsampling_ndgrid(inbounds3), sii); 


%% define a mask to implement aperture growth vs. depth.
% fnum is the F# (aka 'f/D', focal length vs aperture diameter)

fnum=0; %(set fnum=0 to turn off ap. growth).  2.0 is pretty standard for PE 
RxApertureMask = abs(ndZ./(2*(ndX-ndsensX))) >= (fnum);

%channel aperture
%RxApertureMask = abs(ndX-ndsensX) < 1*(0.5*Trans.spacing);
%RxApertureMask(:)=1;

% plot the aperture mask for a given recieve element:
% rx=1;
% imagesc(RxApertureMask(:,:,rx));

distances_flat = distances_ndgrid(inbounds3);
aperture_flat = RxApertureMask(inbounds3);


rf_coefficient =aperture_flat .* distances_flat;

%% FFT & filter

N = depth_m;
NFFT = 2^nextpow2(N); % Next power of 2 

fendidx=NFFT/2+1;
fMHz = Fs/2*linspace(0,1,fendidx)'*1e-6;

harmonicFilt = ones([NFFT 1]);

f0MHz = 1.2;
nharm=ceil(fMHz(end)/f0MHz);
%harmsToRemove=[10];
%nharm=length(harmsToRemove);
tw=200;
%harmonics
for n=nharm
    if n==1
        %dropout=n*f0MHz + [-0.04 0.04];
        dropout=[-0.5 f0MHz+0.1];
    else
        dropout=n*f0MHz + [-0.1 0.1];
    end
    harmonicFilt(1:fendidx) = harmonicFilt(1:fendidx).*(1.0 ./ (1.0 + exp(tw*(fMHz- dropout(1) )))  +  1.0 ./ (1.0 + exp(-tw*(fMHz- dropout(2) ))));
end
%ultraharmonics
for n=2:nharm
    dropout=(n+0.5)*f0MHz + [-0.03 0.03];
    harmonicFilt(1:fendidx) = harmonicFilt(1:fendidx).*(1.0 ./ (1.0 + exp(tw*(fMHz- dropout(1) )))  +  1.0 ./ (1.0 + exp(-tw*(fMHz- dropout(2) ))));
end

harmonicFilt = repmat(harmonicFilt,[1,size(rf_data,2)]);

%plot( fMHz(1:fendidx), harmonicFilt(1:fendidx,5) );

%%
figure(1);
clf;
hold on;
textlab=text( 0.05, 0.9, sprintf('[%d]', 0 ), 'FontSize',18,'FontWeight', 'bold', 'Color', [0.9 0 0], 'Units','normalized' );

imframe=imagesc(zeros([Nz Nx]));
colormap(gray);

set(gca, 'YDir', 'Reverse');
axis equal tight;

tic

framelist=1:nBuff;
framelist=[100:2:1000];
framelist=100:120
nframes=length(framelist);
I_r = zeros([Nz Nx nframes]);

delayed_recieve_stack = zeros(size(delayedsampling_ndgrid));

do_fft_filt=1;
%outer loop is over the frames 
for bi=1:nframes

    fi = framelist(bi);
    
    if do_fft_filt
        page_fft = fft(rf_data(:,:,fi),NFFT,1) ;
        rf_data_page = real(ifft( harmonicFilt.*page_fft ));
        rf_data_page = rf_data_page(1:depth_m,:);
    else
        rf_data_page = rf_data(:,:,fi);
    end

    %framelist(end+1)=bi;
    %--fastest method:
    
    %store the data with delays and aperture applied (and distance
    %multiplied)    
    delayed_recieve_stack(inbounds3) = rf_coefficient .* rf_data_page(rf_sample_indices_flat);
    %sum the delayed channel signals
    I_r(:,:,bi) = sum(delayed_recieve_stack,3);
    
    
%   frame_accum = zeros([Nz Nx]);
%     for rx = 1:sensor_n
%         
%         %mapping from rf sample number to image pixel index
%         rf_sample_indices = delayedsampling_ndgrid(:,:,rx);
%         inbounds =  rf_sample_indices <= depth_m;
% 
%         %this will make rf_sample_indices a flat array of indices
%         rf_sample_indices = rf_sample_indices(inbounds);
%         
%         recieve_aperature_mask = RxApertureMask(:,:,rx);
%         
%         distances = distances_ndgrid(:,:,rx);
%         %frame_accum(inbounds) = frame_accum(inbounds) + recieve_aperature_mask(inbounds) .* distances(inbounds).* rf_data_page(rf_sample_indices, rx);
%         frame_accum(inbounds) = frame_accum(inbounds) + recieve_aperature_mask(inbounds) .* distances(inbounds);
%     end

%I_r(:,:,bi) = frame_accum;
    
    %I_r(:,:,bi) = frame_accum;
    
    %the equivalent for-loop is below
    
%     for n_sens=1:sensor_n
%         u_loc_x = (n_sens-sensor_n/2)*dsensor;
%         u_loc_z = 0;
% 
%         for mz=1:Nz
%             for mx=1:Nx
%                 % compute distance from sensor to point in space |r-r_n|
%                 cur_loc_x = (mx-Nx/2)*dx;
%                 cur_loc_z = (mz)*dz;
% 
%                 cur_dist = sqrt((cur_loc_x-u_loc_x)^2+((cur_loc_z-u_loc_z))^2);
%                 delay =  cur_dist/c ;
%                 
%                 delay_i = round(delay/(dt))+1;
%                 
%                 if delay_i < depth_m
%                     %---> This is the regular beamformer
%                     %I_r(mz,mx,bi)=I_r(mz,mx,bi)+rf_data_page(delay_i,n_sens);
%                 
%                     %---> This is the passive power map 
%                     I_r(mz,mx,bi)=I_r(mz,mx,bi)+ cur_dist*(rf_data_page(delay_i,n_sens));
%                 end;
% 
%             end;
%         end;
%     end;
    toc


    %imagesc(log10(rf_data_page.^2+1));
    %imagesc(log10(I_r(:,:,bi).^2+1));
    %colormap(gray);
    imframe.CData =log10(I_r(:,:,bi).^2+1);
    set(textlab,'String', sprintf('[%d]', fi ) );
    
    drawnow
end


% I_hil = abs(hilbert(I_r(:,:,bi)));
% % normalized
% I_norm = I_hil./max(I_hil(:))+1;
% figure(1);
% imagesc(20*log10(I_norm),'XData',xaxis,'YData',zaxis,[0 6]);
% axis equal;
% axis tight

%If = I_r(:,:,framelist);
pcmap=((mean(I_r.^2,3)));
%%
figure(4);
%gray_r = flipud(gray);

colormap(gray);
subplot(121);
%imagesc((mean(If.^2,3)),'XData',xaxis,'YData',zaxis);
minpc=min(pcmap(:));
maxpc=max(pcmap(:));
imagesc(pcmap, 'XData',xaxis*100,'YData',zaxis*100, [minpc+0.0*(maxpc-minpc) 1.0*maxpc]);
colorbar();
axis equal;
axis tight


xlabel('cm');
ylabel('cm');

subplot(122);

imagesc(log10(pcmap+1),'XData',xaxis*100,'YData',zaxis*100);
colorbar();
xlabel('cm');
ylabel('cm');
axis equal;
axis tight

return
%%

figure(3);

%rf_series = rf_data(:,64,50);

%page_fft = fft( rf_data(:,:,50),NFFT,1)/depth_m;
filtpage = harmonicFilt.*page_fft;
%plot(rf_series);

hold on;

bbAmpVsChans = sum( abs(filtpage(1:fendidx,:)), 1);
bbPowVsChans = sum( abs(filtpage(1:fendidx,:)).^2, 1);

%sigFFT = fft(rf_series(:),NFFT)/depth_m;

plot( fMHz, 2*abs(page_fft(1:fendidx,64)) );
plot( fMHz, 2*abs(filtpage(1:fendidx,64)) );
%plot( fMHz, 2*abs(sigFFT(1:fendidx)) );
xlabel('MHz');


figure(5);
plot(bbAmpVsChans);
xlabel('channel');
ylabel('Broadband amplitude');
text(0.7,0.1,sprintf('sum: %0.2e',sum(bbAmpVsChans)),'Units','Normalized')

%% Spatial filter post-beamformed image

MzFFT = 2^nextpow2(Nz); % Next power of 2 
MxFFT = 2^nextpow2(Nx);

imgFFT = fftshift( fft2( pcmap,MzFFT,MxFFT) );

kZax = (1.0/dz)*linspace(-1,1,MzFFT)';
kXax = (1.0/dx)*linspace(-1,1,MxFFT)';

[vertsZ vertsX] = ndgrid(kZax, kXax);

radialF = sqrt(vertsZ.^2 + vertsX.^2);
lambda0 = 1515/(f0MHz*1e6);
kr0 = 1.0/lambda0;

spatialFilt = ones(size(imgFFT));
nrads = 10;
wr=.1;
for n=0.5:0.5:nrads
    dropout=[-0.1 .1]*kr0 + n*kr0;
    spatialFilt = spatialFilt.*(1.0 ./ (1.0 + exp(wr*(radialF- dropout(1) )))  +  1.0 ./ (1.0 + exp(-wr*(radialF- dropout(2) ))));
end


filtImg = ifft2(imgFFT.*spatialFilt);
%filtImg = ifftshift(ifft2(imgFFT));
%imagesc( kXax, kZax, abs(fftshift(imgFFT) ))
figure(4);
clf;
imagesc( xaxis*100, zaxis*100, (pcmap) );
%imagesc(log10(abs(imgFFT)))
colormap(gray);
axis equal tight;
figure(5);
clf;
imagesc( xaxis*100, zaxis*100, ( abs(filtImg(1:Nz,1:Nx) )) );
%imagesc(spatialFilt)
colormap(gray);
axis equal tight;



