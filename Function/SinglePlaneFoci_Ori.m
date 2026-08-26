function [RawFociMask] = SinglePlaneFoci(SinglemRNA,SinglePlaneNuc,ExFociCore,InFociCore,img_ind,MaxTranscriptionFociMask)
%% identify transcription foci 
% FilFocicore = fspecial('average',5);
% FilFoci = imfilter(SinglemRNA,FilFocicore,'replicate');
Blank = zeros(size(SinglemRNA));
ExFoci = fspecial('gaussian',10,ExFociCore);%10 20
IxFoci = fspecial('gaussian',10,InFociCore);%20 20

outEFoci = imfilter(single(SinglemRNA),ExFoci,'replicate'); 
outIFoci = imfilter(single(SinglemRNA),IxFoci,'replicate'); 
outimsFoci = outEFoci - outIFoci;  
% % figure();
% % imshow(outimsFoci);
 GauoutimsFoci = imgaussfilt(outimsFoci,5);

LowHighFoci = stretchlim(uint16(outimsFoci));
AdjOutimsFoci = imadjust(uint16(outimsFoci),LowHighFoci,[]);


[T EM] = graythresh(AdjOutimsFoci);
FociRawMask = imbinarize(AdjOutimsFoci,T);

se = strel('disk',4);
OpenFociRawMask = imopen(FociRawMask,se);
% figure();
% imshow(OpenFociRawMask,[]);

%identify true transcription foci in single z stack plane with the maximum
%intensity projection result
if sum(sum(SinglePlaneNuc)) == 0
    RawFociMask = Blank;
else
UniSet = MaxTranscriptionFociMask & OpenFociRawMask;

TrueFociSinglePlane = imreconstruct(UniSet,OpenFociRawMask);
RawFociMask = TrueFociSinglePlane;
end
