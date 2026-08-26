function [bw_nuclei] = Nuclei_seg_cjy_eve(EmSingleOri, ExVal, IxVal,low_th,High_th,cir_th,i,NCEst,bwth,GauKernel,fig,Smooth,grayth,graythAdj,bw_inten_th,WaterShed)
%parameters initiation
if isempty(ExVal)
       ExVal = 3;
end

if isempty(IxVal)
    if NCEst == 14
       IxVal = 30;%30
    else
       IxVal = 40;
    end
end
   
if isempty(low_th)
    if NCEst == 14
       low_th = 500;
    else
       low_th = 1000;
    end
end


if isempty(High_th)
    if NCEst == 14
       High_th = 4000;
    else
       High_th = 8000;
    end
end


 if isempty(cir_th)
       cir_th = 0.65; 
 end

 if isempty(bwth)
      bwth = 3000;
 end
 if isempty(GauKernel)
     GauKernel = 2;
 end


 if isempty(Smooth)
     Smooth = 3;
 end


H = fspecial('disk',Smooth);%3
EmSingle = imfilter(EmSingleOri,H,'replicate'); 
% EmSingle = EmSingleOri;
Ex = fspecial('gaussian',100,ExVal);%10
Ix = fspecial('gaussian',100,IxVal);%20

outE = imfilter(single(EmSingle),Ex,'replicate'); 
outI = imfilter(single(EmSingle),Ix,'replicate'); 
outims = outE - outI; 

% LowHighRNA = stretchlim(uint16(outims));
% AdjOutims = imadjust(uint16(outims),LowHighRNA,[]);
% Gauoutims = imgaussfilt(AdjOutims,GauKernel);
if grayth == 1
    [T EM] = graythresh(outims);
    
    T = T+ graythAdj;
    
    RawNucMask = imbinarize(outims,T);
 
else
    RawNucMask = imbinarize(Gauoutims,bwth);
end

figure();
imshow(RawNucMask);
% [T EM] = graythresh(outims);
% NumOfCanNucMask = zeros(1,200);
% ThVal = linspace(0,1,200);
% for IThVal = 1:200
%     CanRawNucMask = imbinarize(outims,ThVal(IThVal));
%     BWCanRawNucMask = bwlabel(CanRawNucMask);
%     StsCanRawNucMask = regionprops(BWCanRawNucMask,'Area');
%     CanNucArea = [StsCanRawNucMask.Area];
%     NumOfCanNucMask(IThVal) = length(CanNucArea);
% end
% figure();
% plot(NumOfCanNucMask);
% ElbowPointIndexRelative = GetElbowIndex(NumOfCanNucMask);
% ThVal(ElbowPointIndexRelative);



% se = strel('disk',5);
% RawNucMask = imopen(RawNucMask,se);
if (NCEst == 14 | WaterShed == 1)

   %watershed-bassed segemntation algorithm with the local minimum mark
   bw2 = ~bwareaopen(~RawNucMask, 5);%10 15
   D = bwdist(~bw2);
   D = -D;
   %local minimum mark
   LocalMinMask = imextendedmin(D,2);%2
   %varified local minima location 
   % imshowpair(D,LocalMinMask,'blend');
   % Modify the distance transform so it only has minima at the desired locations,
   D2 = imimposemin(D,LocalMinMask);
   L = watershed(D2);
   % imshow(label2rgb(L));
   %ovelap image
   % figure();
   % imshow(labeloverlay(AdjEmSingleOri,L,'Transparency',0.6));
   %the logical ture of the watershed segmentation located in the region of
   %the raw nuclei mask
   ZerosCross = RawNucMask & (L==0);
   %assign the position find above to 0 in the mask image 
   RawNucMask(ZerosCross == 1) = 0;
%    RefineNucMask = RawNucMask;
end
% figure();
% imshow(RawNucMask);

RawNucMask = imfill(RawNucMask,'holes');

se = strel('disk',10);

CloseNucMask=imclose(RawNucMask,se);
% se2 = strel('disk',5);
% CloseNucMask=imclose(OpenNucMask,se2);
CloseNucMaskOpen=imopen(CloseNucMask,se);
% figure();
% imshow(CloseNucMaskOpen);

bw_prop = regionprops(CloseNucMaskOpen,EmSingleOri,'Area','Perimeter','MeanIntensity');
bw_area = [bw_prop.Area];
bw_perim = [bw_prop.Perimeter];
bw_meanintensity = [bw_prop.MeanIntensity];
figure();
histogram(bw_area);

%
if isempty(bw_inten_th)
    ind_true = find(bw_area >= low_th & bw_area < High_th & 4*pi*bw_area./bw_perim.^2 >= cir_th);
else
ind_true = find(bw_area >= low_th & 4*pi*bw_area./bw_perim.^2 >= cir_th & bw_meanintensity >= bw_inten_th); 
end

bw_true = ismember(bwlabel(CloseNucMaskOpen),ind_true);
 
%     bw_true = OpenNucMask;

bw_true2 = imfill(bw_true,'holes');


figure();
imshow(bw_true2);
%  saveas(gcf,sprintf('bw_true2%d.tif',i));
%  
if fig ==1

figure();
[B,L] = bwboundaries(bw_true2);
imshow(EmSingleOri,[]);
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
   hold on
end
     hold off      
%      saveas(gcf, sprintf('NucleiWithMask%d.tif',i));

end
bw_nuclei = bw_true2;
end