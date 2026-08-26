function [CytoGFPMask] = SinglePlaneGFPMask_cjy(SingleGFP,SinglePlaneNuc,ExRNACore,IxRNACore,GFPStackIndex)
Ex = fspecial('gaussian',21,ExRNACore);%10
Ix = fspecial('gaussian',21,IxRNACore);%20
outE = imfilter(single(SingleGFP),Ex,'replicate'); 
outI = imfilter(single(SingleGFP),Ix,'replicate'); 
outims = outE - outI;  

% 
% figure();
% imshow(outims,[]);



LowHighGFP = stretchlim(uint16(outims));
AdjOutims = imadjust(uint16(outims),LowHighGFP,[]);



NumOfCanGFPMask = zeros(1,200);
ThVal = linspace(0,1,200);
for i = 1:200
    CanRawGFPMask = imbinarize(AdjOutims,ThVal(i));
    BWCanRawGFPMask = bwlabel(CanRawGFPMask);
    StsCanRawGFPMask = regionprops(BWCanRawGFPMask,'Area');
    CanGFPArea = [StsCanRawGFPMask.Area];
    NumOfCanGFPMask(i) = length(CanGFPArea);
end
% figure();
% axis square
% plot(NumOfCanRNAMask);
ElbowPointIndexRelative = GetElbowIndex(NumOfCanGFPMask)
RawGFPMask = imbinarize(AdjOutims,ThVal(ElbowPointIndexRelative));
RawGFPMaskErode = imerode(RawGFPMask,strel("disk",2));
NucSpotMarker = (RawGFPMaskErode & SinglePlaneNuc);
RemoveGFP = imreconstruct(NucSpotMarker,RawGFPMask);
CytoGFPMask = RawGFPMask.*(~RemoveGFP);
% bw2 = ~bwareaopen(~RawRNAMask,3);
% D = bwdist(~bw2);
% D = -D;
% %local minimum mark
% LocalMinMask = imextendedmin(D,2);
% %varified local minima location 
% % imshowpair(D,LocalMinMask,'blend');
% % Modify the distance transform so it only has minima at the desired locations,
% D2 = imimposemin(D,LocalMinMask);
% L = watershed(D2);
% %the logical ture of the watershed segmentation located in the region of
% %the raw RNA mask
% ZerosCross = RawRNAMask & (L==0);
% %assign the position find above to 0 in the mask image 
% RawRNAMask(ZerosCross == 1) = 0;


% [B,L] = bwboundaries(RawRNAMask);
% Bdr_fig=zeros(2048);
% for k = 1:length(B)
%    boundary = B{k};
%    for j=1:length(boundary)
%     Bdr_fig(boundary(j,1),boundary(j,2)) =1;
%    end
% end
%  
% figure();
% Com = imoverlay(uint16(SinglemRNA),Bdr_fig,'red');
% imshow(Com);



%remove the single mRNA in the nuclei region and transcription foci
% % % % 
% % % % RawRNAMaskErode = imerode(RawRNAMask,strel("disk",1));
% % % % NucSpotMarker = (RawRNAMaskErode & SinglePlaneNuc) | (RawRNAMaskErode & mRNAFociRawMaskPlane);
% % % % RemoveRNA = imreconstruct(NucSpotMarker,RawRNAMask);
% % % % RemNuc = RawRNAMask.*(SinglePlaneNuc);
% % % % TotRemRNA = RemoveRNA|RemNuc;
% % % % CytoRNAMask = RawRNAMask.*(~TotRemRNA);
% % % % % class(CytoRNAMask)
% % % % % class(AdjOutims)
% % % % OutimsAfterMask = double(AdjOutims).*CytoRNAMask;
% % % % 
% % % % RawRNAMask2 = imregionalmax(OutimsAfterMask);
% % % % CytoSingleRNAMask = RawRNAMask2;
% % % % 
% % % % RawRNAMaskDilate2 = imdilate(RawRNAMask2,strel('disk',1));
% % % % figure();
% % % % imshow(RawRNAMaskDilate2);
% NucSpotMarker2 = (RawRNAMask2 & SinglePlaneNuc) | (RawRNAMask2 & mRNAFociRawMaskPlane);
% RemoveRNA2 = imreconstruct(NucSpotMarker2,RawRNAMask2);
% CytoRNAMask2 = RawRNAMask2.*(~RemoveRNA2);
% CytoRNAMask2erode = imdilate(CytoRNAMask2,strel("disk",2));


%plot
% figure();
% imshow(SinglemRNA(1025:2048,1025:2048),[]);
% hold on
% [B,L] = bwboundaries(CytoRNAMask(1025:2048,1025:2048));
% [BN,LN] = bwboundaries(SinglePlaneNuc(1025:2048,1025:2048));
% [BR,LR] = bwboundaries(mRNAFociRawMaskPlane(1025:2048,1025:2048));
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
%    hold on
% end
% for k = 1:length(BN)
%    boundary = BN{k};
%    plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 0.5);
%    hold on
% end
% for k = 1:length(BR)
%    boundary = BR{k};
%    plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 0.7);
%    hold on
% end
% hold off     
% saveas(gcf, sprintf('RawCytoRNAMask%d.tif',IndexOfZstack));

% (1025:2048,1025:2048)
% figure();
% imshow(Adj,[]);
% hold on
% [B,L] = bwboundaries(RawGFPMask);
% [BN,LN] = bwboundaries(SinglePlaneNuc);
% 
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
%    hold on
% end
% for k = 1:length(BN)
%    boundary = BN{k};
%    plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 0.5);
%    hold on
% end
% hold off  
% LLCytoSingleRNAMask = bwlabel(CytoSingleRNAMask);
% CytoSingleRNAMaskStats = regionprops(LLCytoSingleRNAMask,SinglemRNA,'MeanIntensity');
% CytoSingleRNAMeanIntensity = [CytoSingleRNAMaskStats.MeanIntensity];
% figure();
% h = histogram(CytoSingleRNAMeanIntensity);
% [HisVal, HisEdge] = histcounts(CytoSingleRNAMeanIntensity,'Normalization','probability');
% figure();
% plot(HisVal);
[B,L] = bwboundaries(CytoGFPMask);
Bdr_fig=zeros(2048);
for k = 1:length(B)
   boundary = B{k};
   for j=1:length(boundary)
    Bdr_fig(boundary(j,1),boundary(j,2)) =1;
   end
end
figure();%"Visible","off"
Com = imoverlay(uint16(SingleGFP),Bdr_fig,'red');
imshow(Com);
saveas(gcf, sprintf('RawGFPMaskZstack%d.tif',GFPStackIndex));

end
