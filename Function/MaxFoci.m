function [MaxTranscriptionFoci,MaxTranscriptionFociMask] = MaxFoci(MIP,ExFociCore,InFociCore,seval)
% %maximum intensity projection foci
% figure();
% imshow(MIP(1024:2048,1024:2048),[]);
MaxFilFocicore = fspecial('average',3);
MaxFilFoci = imfilter(MIP,MaxFilFocicore,'replicate');
ExFoci = fspecial('gaussian',40,ExFociCore);%10
IxFoci = fspecial('gaussian',40,InFociCore);%20

outEFociMax = imfilter(single(MaxFilFoci),ExFoci,'replicate'); 
outIFociMax = imfilter(single(MaxFilFoci),IxFoci,'replicate'); 
outimsFociMax = outEFociMax - outIFociMax;  
% figure();
% imshow(outimsFociMax(1024:2048,1024:2048),[]);
% GauoutimsFoci = imgaussfilt(MIP,10);
LowHighFociMax = stretchlim(uint16(outimsFociMax));
AdjOutimsFociMax = imadjust(uint16(outimsFociMax),LowHighFociMax,[]);
% figure();
% imshow(AdjOutimsFociMax(1024:2048,1024:2048));
[T EM] = graythresh(AdjOutimsFociMax);
FociRawMaxMask = imbinarize(AdjOutimsFociMax,T);
se = strel('disk',seval);
OpenFociRawMaxMaskRaw = imopen(FociRawMaxMask,se);
% OpenFociRawMaxMaskMarker = OpenFociRawMaxMaskRaw & MaxNuc_true;
% OpenFociRawMaxMask = imreconstruct(OpenFociRawMaxMaskMarker,OpenFociRawMaxMaskRaw);
% 
MaxTranscriptionFociMask = OpenFociRawMaxMaskRaw;
MaxTranscriptionFoci = MIP.*OpenFociRawMaxMaskRaw;
end
