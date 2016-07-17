%% repackage saved data

% pre-allocate rf data
rf_data = zeros([1792 128 50*100]);

for m=1:100
    m
    for k=1:50
        my_i = (m-1)*50+k;
        rf_data(:,:,my_i) = squeeze(chandat2(:,k,:,m));
    end
end;

%%
for k=1:1000
    imagesc(log10(1+abs(hilbert(rf_data(:,:,k)))));colormap(gray);
    pause(0.01);
end;
