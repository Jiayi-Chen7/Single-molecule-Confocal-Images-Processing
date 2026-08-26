function [CytoSingleGFPMask,GauGFPFil] = SinglePlaneGFPMaskAdjTh_cjy0(SingleGFP,SinglePlaneNuc,ExGFPCore,IxGFPCore,IndexOfZstack,GauKer)
%%
% filcore = fspecial('average',3);
% filSinglemRNA = imfilter(SinglemRNA,filcore,'replicate');
Ex = fspecial('gaussian',21,ExGFPCore);%10
Ix = fspecial('gaussian',21,IxGFPCore);%20
outE = imfilter(single(SingleGFP),Ex,'replicate'); 
outI = imfilter(single(SingleGFP),Ix,'replicate'); 
outims = outE - outI;  

% 
% figure();
% imshow(outims,[]);

 if isempty(GauKer)
      GauKer = 1.5;
 end

Gauoutims = imgaussfilt(single(outims),GauKer);
% figure();
% imshow(outims(1536:2048,1536:2048),[]);
LowHighRNA = stretchlim(uint16(Gauoutims));
AdjOutims = imadjust(uint16(Gauoutims),LowHighRNA,[]);
GauGFPFil = imgaussfilt(single(SingleGFP),1.5);
% figure();
% imshow(AdjOutims(1536:2048,1536:2048));

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
% plot(NumOfCanGFPMask,'LineWidth',2);
ElbowPointIndexRelative = GetElbowIndex(NumOfCanGFPMask)
RawGFPMask = imbinarize(AdjOutims,ThVal(ElbowPointIndexRelative));

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

RawGFPMaskErode = imerode(RawGFPMask,strel("disk",1));
NucSpotMarker = (RawGFPMaskErode & SinglePlaneNuc);
RemoveGFP = imreconstruct(NucSpotMarker,RawGFPMask);
RemNuc = RemoveGFP.*(SinglePlaneNuc);
TotRemGFP = RemoveGFP|RemNuc;
CytoGFPMask = RawGFPMask.*(~TotRemGFP);
% class(CytoRNAMask)
% class(AdjOutims)
% OutimsAfterMask = double(AdjOutims).*CytoRNAMask;
% 
% RawRNAMask2 = imregionalmax(OutimsAfterMask);

CytoSingleGFPMask = CytoGFPMask;
% 
% RawRNAMaskDilate2 = imdilate(RawRNAMask2,strel('disk',1));
% figure();
% imshow(RawRNAMaskDilate2);
% NucSpotMarker2 = (RawRNAMask2 & SinglePlaneNuc) | (RawRNAMask2 & mRNAFociRawMaskPlane);
% RemoveRNA2 = imreconstruct(NucSpotMarker2,RawRNAMask2);
% CytoRNAMask2 = RawRNAMask2.*(~RemoveRNA2);
% CytoRNAMask2erode = imdilate(CytoRNAMask2,strel("disk",2));

% CytoSingleRNAMask = CytoRNAMask;
%plot
% % % % % LowHighRNA2 = stretchlim(uint16(SingleGFP));
% % % % % AdjOutims2 = imadjust(uint16(SingleGFP),LowHighRNA2,[]);
% % % % % figure();
% % % % % imshow(AdjOutims2(1536:2048,1536:2048),[]);
% % % % % hold on
% % % % % [B,L] = bwboundaries(CytoSingleGFPMask(1536:2048,1536:2048));
% % % % % 
% % % % % for k = 1:length(B)
% % % % %    boundary = B{k};
% % % % %    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
% % % % %    hold on
% % % % % end
% % % % % 
% % % % % hold off    
% saveas(gcf, sprintf('RawCytoRNAMask%d.tif',IndexOfZstack));

% 
% figure();
% imshow(SingleGFP(1:1024,1:1024),[]);
% hold on
% [B,L] = bwboundaries(CytoSingleGFPMask(1:1024,1:1024));
% [BN,LN] = bwboundaries(SinglePlaneNuc(1:1024,1:1024));
% % [BR,LR] = bwboundaries(mRNAFociRawMaskPlane(1025:2048,1025:2048));
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
% 
% hold off  
% LLCytoSingleRNAMask = bwlabel(CytoSingleRNAMask);
% CytoSingleRNAMaskStats = regionprops(LLCytoSingleRNAMask,SinglemRNA,'MeanIntensity');
% CytoSingleRNAMeanIntensity = [CytoSingleRNAMaskStats.MeanIntensity];
% figure();
% h = histogram(CytoSingleRNAMeanIntensity);
% [HisVal, HisEdge] = histcounts(CytoSingleRNAMeanIntensity,'Normalization','probability');
% figure();
% plot(HisVal);
% [B,L] = bwboundaries(CytoSingleRNAMask);
% Bdr_fig=zeros(2048);
% for k = 1:length(B)
%    boundary = B{k};
%    for j=1:length(boundary)
%     Bdr_fig(boundary(j,1),boundary(j,2)) =1;
%    end
% end
% figure();%"Visible","off"
% Com = imoverlay(uint16(SinglemRNA),Bdr_fig,'red');
% imshow(Com);
% saveas(gcf, sprintf('RawRNAMaskZstack%d.tif',IndexOfZstack));

end
