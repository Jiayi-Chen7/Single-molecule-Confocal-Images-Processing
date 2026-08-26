function [MultiLayerMask,mask_out] = MultiLayerSpotIdentifyTwoLayer(RawMaskstack,Nr)
%modified from spmask.m from Xu, Heng, et al.
% "Combining protein and mRNA quantification to decipher transcriptional regulation." 
% Nature methods 12.8 (2015): 739-742.
% max3 = imregionalmax(RawMaskstack);   %%% 3D local maxima
% max2 = false(size(max3));   %%% 2D local maxima
% max20 = max2;
% N_layer = size(RawMaskstack,3);
% for I_layer = 1:N_layer
%     
% %     max20(:,:,I_layer) = imregionalmax(RawMaskstack(:,:,I_layer));
% %     max2(:,:,I_layer) = imdilate(max20(:,:,I_layer),strel('square',3));
% end
N_layer = size(RawMaskstack,3);
max3 = RawMaskstack;  
max3_temp = RawMaskstack; 
max3_temp2 = RawMaskstack; 
max3_temp3 = RawMaskstack; 
max2 = zeros(size(RawMaskstack)); 
for I_layer = 1:N_layer
    max2(:,:,I_layer) = imdilate(RawMaskstack(:,:,I_layer),strel('square',3));
%     max3(:,:,I_layer) = imdilate(RawMaskstack(:,:,I_layer),strel('square',3));
end

Nr = ceil(Nr);
Ar = -Nr:Nr;
zmin = max(1,(1+Ar));
zmax = min(N_layer,(N_layer+Ar));
% for Ir = 1:length(Ar)
Ir = 1;
    max3_temp(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
% Ir = 2;
%     max3_temp2(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
Ir = 3;
    max3_temp3(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
% end
max3 = max3_temp| max3_temp3;

mask_out = max3;
% class(mask_out)
% class(RawMaskstack)
RawMaskstack = logical(RawMaskstack);
% MultiLayerMask = zeros(size(mask_out));
% for j = 1:N_layer
%     MultiLayerMask(:,:,j) = imreconstruct(mask_out(:,:,j),RawMaskstack(:,:,j));
% end
MultiLayerMask = imreconstruct(mask_out,RawMaskstack);
end
