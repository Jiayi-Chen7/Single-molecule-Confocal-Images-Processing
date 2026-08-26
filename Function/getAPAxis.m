function [AP_x, AP_y] = getAPAxis(embryoMask, old_AP_x, old_AP_y,flipAP)
% find AP axis by searching two points  in the perimeter
% of the embryo mask with the largest distance. be caution, sperm extrusion
% tips can mess up the search
[ycoords, xcoords] = find(bwperim(embryoMask));
len = length(xcoords);
if length(xcoords) > 2500, xcoords = xcoords(round(1:len/2500:end)); ycoords = ycoords(round(1:len/2500:end)); len = length(xcoords); end
xcoords = single(xcoords); ycoords = single(ycoords);
len = length(xcoords);
xdist = xcoords * ones(1,len,'single') - ones(len,1,'single') * xcoords';
ydist = ycoords * ones(1,len,'single') - ones(len,1,'single') * ycoords';
distances = xdist.^2 + ydist.^2;
topDistances = distances > 0.995*max(distances(:));

%stats = regionprops(embryoMask,'MajorAxisLength');
%APlength=stats.MajorAxisLength;
%topDistances = distances > .995*APlength;
%topDistances = distances < 1.005*APlength;
topDistancesIdx = find(topDistances);
[distancesSorted, idx] = sort(distances(topDistances));
m = zeros(length(idx),1); n = zeros(length(idx),1);
for i=1:length(idx), [m(i),n(i)] = ind2sub(size(distances),topDistancesIdx(idx(i))); end
mn = [m; n];
X = [xcoords(mn) ycoords(mn)];
kIdx = kmeans(X, 2, 'emptyaction', 'singleton', 'display', 'off');
coordIdx1 = mn(kIdx==1); coordIdx2 = mn(kIdx==2);
AP_x = double([mean(xcoords(coordIdx1)), mean(xcoords(coordIdx2))]);
AP_y = double([mean(ycoords(coordIdx1)), mean(ycoords(coordIdx2))]);
APMask = embryoMask & line2mask(AP_x, AP_y, size(embryoMask,1), size(embryoMask,2));
cols = sum(APMask, 1);
if flipAP == 1
% AP_x(1) = find(cols, 1, 'first');
% AP_x(2) = find(cols, 1, 'last');
    AP_x(2) = find(cols, 1, 'first');
    AP_x(1) = find(cols, 1, 'last');
else
    AP_x(1) = find(cols, 1, 'first');
   AP_x(2) = find(cols, 1, 'last');
end
rows1 = APMask(:, AP_x(1));
rows2 = APMask(:, AP_x(2));
% rows1 = APMask(:, AP_x(2));
% rows2 = APMask(:, AP_x(1));
if find(rows1, 1, 'first') <= find(rows2, 1, 'first')
    AP_y(1) = find(rows1, 1, 'first');
    AP_y(2) = find(rows2, 1, 'last');
else
    AP_y(1) = find(rows1, 1, 'last');
    AP_y(2) = find(rows2, 1, 'first');
end
if exist('old_AP_y') && ~isempty(old_AP_y)
    dist = (AP_x(1) - old_AP_x(1))^2 + (AP_y(1) - old_AP_y(1))^2;
    distFlipped = (AP_x(1) - old_AP_x(2))^2 + (AP_y(1) - old_AP_y(2))^2;
    if distFlipped < dist, AP_x = fliplr(AP_x); AP_y = fliplr(AP_y); end
end
end