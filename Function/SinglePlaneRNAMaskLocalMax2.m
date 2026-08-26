function [CytoSingleRNAMask,GauFil] = SinglePlaneRNAMaskLocalMax2(SinglemRNA,SinglePlaneNuc,ExRNACore,IxRNACore,IndexOfZstack,mRNAFociRawMaskPlane)
%%
% filcore = fspecial('average',3);
% filSinglemRNA = imfilter(SinglemRNA,filcore,'replicate');
RemoveRegion = SinglePlaneNuc| mRNAFociRawMaskPlane;
SinglemRNA = double(SinglemRNA).*(~RemoveRegion);

Ex = fspecial('gaussian',21,ExRNACore);%10
Ix = fspecial('gaussian',21,IxRNACore);%20
outE = imfilter(single(SinglemRNA),Ex,'replicate'); 
outI = imfilter(single(SinglemRNA),Ix,'replicate'); 
outims = outE - outI;  

LowHighRNARaw1 = stretchlim(uint16(SinglemRNA));
AdjRNARaw1 = imadjust(uint16(SinglemRNA),LowHighRNARaw1,[]);
% figure();
% imshow(AdjRNARaw1(1:300,1:300));

% % 
% figure();
% imshow(outims,[]);


Gauoutims = imgaussfilt(single(outims),1.5);
% figure();
% imshow(Gauoutims,[]);
LowHighRNA = stretchlim(uint16(Gauoutims));

AdjOutims = imadjust(uint16(Gauoutims),LowHighRNA,[]);
%20240713
%GauFil = imgaussfilt(single(SinglemRNA),1.5);
GauFil = AdjOutims;

% figure();
% imshow(AdjOutims);

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
% plot(NumOfCanRNAMask,'LineWidth',2);
% xlabel("Threshold index");
% set(gca,'fontname','Arial');
% ylabel("The number of masks")
% set(gca,'fontname','Arial');
% set (gca,'linewidth',2,'fontsize',12);
% set(gcf,'position',[100 100 300 300]);



ElbowPointIndexRelative = GetElbowIndex(NumOfCanRNAMask(1:end-1));
% ThVal(22)
RawRNAMask = imbinarize(AdjOutims,ThVal(ElbowPointIndexRelative));



    
% figure();
% imshow(RawRNAMask);
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

RawRNAMaskDilate = imdilate(RawRNAMask,strel("disk",1));
% RawRNAMaskErode = RawRNAMask;
NucSpotMarker = (RawRNAMaskDilate & SinglePlaneNuc) | (RawRNAMaskDilate & mRNAFociRawMaskPlane);
RemoveRNA = imreconstruct(NucSpotMarker,RawRNAMask);
RemNuc = RawRNAMask.*(SinglePlaneNuc);
TotRemRNA = RemoveRNA|RemNuc;
CytoRNAMask = RawRNAMask.*(~TotRemRNA);
% class(CytoRNAMask)
% class(AdjOutims)
%%
%plot
% % RawRNAAfterMask = double(SinglemRNA).*CytoRNAMask;
% % LowHighRNARaw = stretchlim(uint16(RawRNAAfterMask));
% % AdjRNARaw = imadjust(uint16(RawRNAAfterMask),LowHighRNARaw,[]);
% % figure();
% % imshow(AdjRNARaw(1:300,1:300));

%%
OutimsAfterMask = double(SinglemRNA).*CytoRNAMask;

RawRNAMask2 = imregionalmax(OutimsAfterMask);

CytoSingleRNAMask = RawRNAMask2;
%%
CytoRNAImg = double(Gauoutims).*CytoSingleRNAMask;
CytoRNAImg_Hmax = imhmax(CytoRNAImg,400);
CytoSingleRNAMask2 = imregionalmax(CytoRNAImg_Hmax,26);
%%
% LowHighRNARaw1 = stretchlim(uint16(SinglemRNA));
% AdjRNARaw1 = imadjust(uint16(SinglemRNA),LowHighRNARaw1,[]);
% CytoRNALL = bwlabel(CytoSingleRNAMask(1:300,1:300));
% CytoRNALLsts = regionprops3(CytoRNALL,SinglemRNA(1:300,1:300),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
% size(CytoRNACen)
% % figure();
% % imshow(AdjRNARaw1(1:300,1:300));
% % hold on
% % scatter(CytoRNACen(:,1),CytoRNACen(:,2),4,'r','filled');
% % hold off
% figure();
% LowHighRNA2 = stretchlim(uint16(SinglemRNA));
% AdjOutims2 = imadjust(uint16(SinglemRNA),LowHighRNA2,[]);
% imshow(AdjOutims2);
% hold on
% GoodRNAMask2D =RawRNAMask2;
% GoodRNAMask2DLL = bwlabel(GoodRNAMask2D);
% GoodRNAMask2Dsts = regionprops3(GoodRNAMask2DLL,SinglemRNA,'WeightedCentroid');
% GoodRNAMask2DWC = [GoodRNAMask2Dsts.WeightedCentroid];
% scatter(GoodRNAMask2DWC(:,1),GoodRNAMask2DWC(:,2),14,'r','filled');
% hold off

% figure();
% imshow(CytoSingleRNAMask);


% 
% RawRNAMaskDilate2 = imdilate(RawRNAMask2,strel('disk',1));
% figure();
% imshow(RawRNAMaskDilate2);
% NucSpotMarker2 = (RawRNAMask2 & SinglePlaneNuc) | (RawRNAMask2 & mRNAFociRawMaskPlane);
% RemoveRNA2 = imreconstruct(NucSpotMarker2,RawRNAMask2);
% CytoRNAMask2 = RawRNAMask2.*(~RemoveRNA2);
% CytoRNAMask2erode = imdilate(CytoRNAMask2,strel("disk",2));

% CytoSingleRNAMask = CytoRNAMask;
% plot
% % % % LowHighRNA2 = stretchlim(uint16(SinglemRNA));
% % % % AdjOutims2 = imadjust(uint16(SinglemRNA),LowHighRNA2,[]);
% % % % figure();
% % % % Img = uint16(SinglemRNA);
% % % % imshow(Img(512:1024,512:1024),[]);
% % % % 
% % % % figure();
% % % % imshow(AdjOutims2(512:1024,512:1024),[]);
% % % % hold on
% % % % [B,L] = bwboundaries(CytoRNAMask(512:1024,512:1024));
% % % % [BN,LN] = bwboundaries(SinglePlaneNuc(512:1024,512:1024));
% % % % [BR,LR] = bwboundaries(mRNAFociRawMaskPlane(512:1024,512:1024));
% % % % for k = 1:length(B)
% % % %    boundary = B{k};
% % % %    plot(boundary(:,2), boundary(:,1), 'color','r', 'LineWidth', 0.3);
% % % %    hold on
% % % % end
% % % % for k = 1:length(BN)
% % % %    boundary = BN{k};
% % % %    plot(boundary(:,2), boundary(:,1), 'color','c', 'LineWidth', 2);
% % % %    hold on
% % % % end
% % % % for k = 1:length(BR)
% % % %    boundary = BR{k};
% % % %    plot(boundary(:,2), boundary(:,1),'color','y', 'LineWidth', 1.5);
% % % %    hold on
% % % % end
% % % % hold off     
% % % % % % saveas(gcf, sprintf('RawCytoRNAMask%d.tif',IndexOfZstack));
% % % % % 
% % % % % % 
% % % % % % figure();
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
