%% Compute drawdown using U-Net 
for i=101:151
figure(i);
imagesc(reshape(Data(:,:,end,i),[128 128])'); colorbar(); caxis([-2 1])
end