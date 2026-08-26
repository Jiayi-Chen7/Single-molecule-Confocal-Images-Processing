function [RawFociMask] = SinglePlaneFoci(SinglemRNA,SinglePlaneNuc,ExFociCore,InFociCore,img_ind,SeindInput,MaxTranscriptionFociMask)
%% identify transcription foci 
% FilFocicore = fspecial('average',5);
% FilFoci = imfilter(SinglemRNA,FilFocicore,'replicate');
Blank = zeros(size(SinglemRNA));
ExFoci = fspecial('gaussian',40,ExFociCore);%10%40
IxFoci = fspecial('gaussian',40,InFociCore);%20%40

outEFoci = imfilter(single(SinglemRNA),ExFoci,'replicate'); 
outIFoci = imfilter(single(SinglemRNA),IxFoci,'replicate'); 
outimsFoci = outEFoci - outIFoci;  
% % figure();
% % imshow(outimsFoci);
% GauoutimsFoci = imgaussfilt(outimsFoci,5);

LowHighFoci = stretchlim(uint16(outimsFoci));
AdjOutimsFoci = imadjust(uint16(outimsFoci),LowHighFoci,[]);


[T EM] = graythresh(AdjOutimsFoci);

FociRawMask = imbinarize(AdjOutimsFoci,T);
seInd = SeindInput;
se = strel('disk',seInd);
OpenFociRawMask = imopen(FociRawMask,se);
% figure();
% imshow(OpenFociRawMask,[]);

%identify true transcription foci in single z stack plane with the maximum
%intensity projection result
if sum(sum(SinglePlaneNuc)) == 0
    RawFociMask = Blank;
else
UniSet = MaxTranscriptionFociMask & OpenFociRawMask & SinglePlaneNuc;
seUni = strel('disk',1);
OpenUniSet = imopen(UniSet,seUni);
 TrueFociSinglePlane = imreconstruct(OpenUniSet,OpenFociRawMask);
% TrueFociSinglePlane = imreconstruct(UniSet,MaxTranscriptionFociMask);

RawFociMask = TrueFociSinglePlane;
end
