

%% setup variables
c=1480; % [m/s] assume tissue
Trans.spacing=params.pitch;
Trans.frequency=params.fs/4*1e-6;
lambda=c/(Trans.frequency*1e6);

[depth_m sensor_n, nBuff] = size(rf_data);

sampleZero=0; %should be zero if the physical image begins at sample 1

Fs = 4*Trans.frequency*1e6;
dt = 1/Fs;

%lambda=c/1.1e6;

%dx = lambda/2;
dz = lambda/2;
%dz =c*dt;
dx=Trans.spacing;
Nx = round(sensor_n*Trans.spacing/dx);
%Nz = 450;
Nz = round( (depth_m-sampleZero)/2.0*c*dt/dz);

Nz = round( (6e-2)/dz);

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

fnum=0.0; %(set fnum=0 to turn off ap. growth).  2.0 is pretty standard for PE 
RxApertureMask = abs(ndZ./(2*(ndX-ndsensX))) >= (fnum);

% plot the aperture mask for a given recieve element:
%rx=64;
%imagesc(RxApertureMask(:,:,rx));

distances_flat = distances_ndgrid(inbounds3);
aperture_flat = RxApertureMask(inbounds3);


rf_coefficient =aperture_flat .* distances_flat;

%%
figure(1);
tic

framelist=1:nBuff;
framelist=[16];

nframes=length(framelist);
I_r = zeros([Nz Nx nframes]);

delayed_recieve_stack = zeros(size(delayedsampling_ndgrid));


%outer loop is over the frames 
for bi=1:nframes

    fi = framelist(bi);
    rf_data_page = rf_data(:,:,fi);
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


    imagesc(log10(rf_data_page.^2+1));
    %imagesc(log10(I_r(:,:,bi).^2+1));
    colormap(gray);
    drawnow
end

%%

% I_hil = abs(hilbert(I_r(:,:,bi)));
% % normalized
% I_norm = I_hil./max(I_hil(:))+1;
% figure(1);
% imagesc(20*log10(I_norm),'XData',xaxis,'YData',zaxis,[0 6]);
% axis equal;
% axis tight

%If = I_r(:,:,framelist);
pcmap=mean(I_r.^2,3);
figure(2);
colormap(gray);
subplot(121);
%imagesc((mean(If.^2,3)),'XData',xaxis,'YData',zaxis);
minpc=min(pcmap(:));
maxpc=max(pcmap(:));
imagesc(pcmap, 'XData',xaxis,'YData',zaxis, [minpc+0.0*(maxpc-minpc) 0.5*maxpc]);
colorbar();
axis equal;
axis tight

subplot(122);

imagesc(log10(pcmap+1),'XData',xaxis,'YData',zaxis);
colorbar();

axis equal;
axis tight

return
%%

figure(3);

rf_series = rf_data(:,62,10);

%plot(rf_series);

hold on;
N = length(rf_series);
NFFT = 2^nextpow2(N); % Next power of 2 
sigFFT = fft(rf_series,NFFT)/N;

fendidx=NFFT/2+1;
fMHz = Fs/2*linspace(0,1,fendidx)*1e-6;

plot( fMHz, 2*abs(sigFFT(1:fendidx)) );











