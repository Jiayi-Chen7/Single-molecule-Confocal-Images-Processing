function [MultiLayerNucMask] = MultiLayerNucIdentify(DapiRawMask,Nr,RawData)
%modified from spmask.m from Xu, Heng, et al.
% "Combining protein and mRNA quantification to decipher transcriptional regulation." 
% Nature methods 12.8 (2015): 739-742.

N_layer = size(DapiRawMask,3);
max3 = DapiRawMask;   
max2 = zeros(size(DapiRawMask)); 
for I_layer = 1:N_layer
    max2(:,:,I_layer) = imdilate(DapiRawMask(:,:,I_layer),strel('square',3));
%     max3(:,:,I_layer) = imdilate(RawMaskstack(:,:,I_layer),strel('square',3));
end

Nr = ceil(Nr);
Ar = -Nr:Nr;
zmin = max(1,(1+Ar));
zmax = min(N_layer,(N_layer+Ar));
for Ir = 1:length(Ar)
    max3(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
end

mask_out = max3;
MultiLayerNucMask = zeros(size(mask_out));
for j = 1:N_layer
    MultiLayerNucMask(:,:,j) = imreconstruct(mask_out(:,:,j),DapiRawMask(:,:,j));
figure();
[B,L] = bwboundaries(MultiLayerNucMask(:,:,j));
imshow(RawData(:,:,j),[]);
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
   hold on
end
     hold off      
     saveas(gcf, sprintf('RefineNucleiWithMask%d.tif',j));

end

end
