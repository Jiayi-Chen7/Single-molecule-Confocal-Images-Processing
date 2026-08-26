function [CytoSingleRNAMask,GauFil] = SinglePlaneRNAMaskRaw_cjy(SinglemRNA,SinglePlaneNuc,ExRNACore,IxRNACore,IndexOfZstack,mRNAFociRawMaskPlane)
%%
% filcore = fspecial('average',3);
% filSinglemRNA = imfilter(SinglemRNA,filcore,'replicate');
Ex = fspecial('gaussian',21,ExRNACore);%10
Ix = fspecial('gaussian',21,IxRNACore);%20
outE = imfilter(single(SinglemRNA),Ex,'replicate'); 
outI = imfilter(single(SinglemRNA),Ix,'replicate'); 
outims = outE - outI;  

% 
% figure();
% imshow(outims,[]);


Gauoutims = imgaussfilt(single(outims),1.5);
LowHighRNA = stretchlim(uint16(Gauoutims));
AdjOutims = imadjust(uint16(Gauoutims),LowHighRNA,[]);
GauFil = imgaussfilt(single(SinglemRNA),1.5);


NumOfCanRNAMask = zeros(1,200);
ThVal = linspace(0,1,200);
for i = 1:200
    CanRawRNAMask = imbinarize(AdjOutims,ThVal(i));
    BWCanRawRNAMask = bwlabel(CanRawRNAMask);
    StsCanRawRNAMask = regionprops(BWCanRawRNAMask,'Area');
    CanRNAArea = [StsCanRawRNAMask.Area];
    NumOfCanRNAMask(i) = length(CanRNAArea);
end
% figure();
% axis square
% plot(NumOfCanRNAMask);
ElbowPointIndexRelative = GetElbowIndex(NumOfCanRNAMask);
RawRNAMask = imbinarize(AdjOutims,ThVal(ElbowPointIndexRelative));

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

RawRNAMaskErode = imerode(RawRNAMask,strel("disk",1));
NucSpotMarker = (RawRNAMaskErode & SinglePlaneNuc) | (RawRNAMaskErode & mRNAFociRawMaskPlane);
RemoveRNA = imreconstruct(NucSpotMarker,RawRNAMask);
RemNuc = RawRNAMask.*(SinglePlaneNuc);
TotRemRNA = RemoveRNA|RemNuc;
CytoRNAMask = RawRNAMask.*(~TotRemRNA);
% class(CytoRNAMask)
% class(AdjOutims)
% % % OutimsAfterMask = double(AdjOutims).*CytoRNAMask;
% % % 
% % % RawRNAMask2 = imregionalmax(OutimsAfterMask);

CytoSingleRNAMask = CytoRNAMask;
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

% 
% figure();
% imshow(SinglemRNA(1025:2048,1025:2048),[]);
% hold on
% [B,L] = bwboundaries(RawRNAMaskDilate2(1025:2048,1025:2048));
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
