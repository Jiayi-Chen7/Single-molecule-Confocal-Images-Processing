function[x_top,y_top,x_back,y_back]=get_Axis(I_mask)

% INPUT
%
% I_embryo - matrix containing the image of the embryo
% I_mask - matrix containg the mask of the embryo
%
% OUTPUT
%
% x_top,y_top,x_back,y_back - positions of the top and the back of the
% AP axis
%
% REVISION HISTORY
%
% 05/24/10 Written by Julien Dubuis

L=bwlabel(I_mask);
Props=regionprops(L,'Extrema');

ext=Props.Extrema;

x_top=round(ext(8,1));
y_top=round(0.5*(ext(7,2)+ext(8,2)));

y_back=y_top;

L2=bwlabel(I_mask(y_top,:));
Props2=regionprops(L2,'Extrema');

ext2=Props2.Extrema;

x_back=round(ext2(3,1));
save AP2

% figure, imshow(imadjust(I_embryo)), title('Please double-click on the head of the embryo and press Enter')
% 
% pause
% 
% [x y]=getpts;
% 
% x_top_bis=round(x);
% y_top_bis=round(y);
% 
% close
% 
% figure, imshow(imadjust(I_embryo)), title('Please double-click on the tail of the embryo and press Enter')
% 
% pause
% 
% [x y]=getpts;
% 
% x_back_bis=round(x);
% y_back_bis=round(y);
% 
% close
end
