function [ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn] = convertXY2AP(intens, nx, ny, emmask, fig,flipAP) %(NucleiIntensity, NucleiPosition(:,1), NucleiPosition(:,2), embryomask, fig)
    %% calculates what is ventral, dorsal from the coordinates
    % takes intensities of nuclei at (nx,ny); and embryo mask
    % col changes the color, 
    % orientation=1 will set the DV axis so the left side is
    % Dorsal(flat),anterior ple face down
    % ventralData and dorsalData in format [AP, DV, intensity, X, Y]
    % last modified July 13 2009
    
    dataind=1:length(intens);
    dataind=dataind';
    
    % test the first orientation of AP axis
    [AP_x, AP_y] = getAPAxis(emmask,[],[],flipAP);
    AP_n = myXY2AP(nx, ny, AP_x, AP_y);

    
    dorsalPts = AP_n(:,2)<0;

%     dorsalData = sortrows([AP_n(dorsalPts, 1), intens(dorsalPts)]);

    % flip the orientation of the AP if Anteripole face down
% %     if dorsalData(1,2)<dorsalData(size(dorsalData,1),2),
% %         AP_x = fliplr(AP_x); AP_y = fliplr(AP_y);
% %         AP_n = myXY2AP(nx, ny, AP_x, AP_y);
% %         dorsalPts = AP_n(:,2)<0;
% %     end
% %     
     
    ventralPts = AP_n(:,2)>=0;

    [vp dp] = deal(AP_n(ventralPts,2), AP_n(dorsalPts,2));
    % flip the orientation of the DV axis depending on the var orientation,
    % nuclei on the dorsal side are close to AP axis. However, this may not
    % work for dorsal view when the shape of the embryo is actrully
    % symetric
    if max(abs(vp)) < max(abs(dp))
        tmp = dorsalPts;
        dorsalPts = ventralPts;
        ventralPts = tmp;
    end
    
    
    
    dorsalData = sortrows([AP_n(dorsalPts, :) intens(dorsalPts) nx(dorsalPts) ny(dorsalPts) dataind(dorsalPts)]);
    ventralData = sortrows([AP_n(ventralPts, :) intens(ventralPts) nx(ventralPts) ny(ventralPts) dataind(ventralPts)]);
    
    % scale the data to [0, 1] range along the AP axis
    v = [diff(AP_x), diff(AP_y)];
    APlength=norm(v);
    %APlength=max(max(dorsalData(:,1)),max(ventralData(:,1)));
    dorsalData(:,1) = dorsalData(:,1) / APlength;
    ventralData(:,1) = ventralData(:,1) / APlength;
    NormAPn = AP_n(:,1)./APlength;
    
    %dorsalData(:,1) = dorsalData(:,1) / max(dorsalData(:,1));
   % ventralData(:,1) = ventralData(:,1) / max(ventralData(:,1));

    % some visualization of which side we are taking as the ventral and
    % anteripole
    if fig
        hold on;
        plot(nx(ventralPts), ny(ventralPts), 'b.');
        line(AP_x, AP_y, 'Color',[.8 .8 .0]);

        v = [diff(AP_x), diff(AP_y)];
        n = [-v(2) v(1)]; n = n/norm(n);
        s = 10;
        p = repmat(0:s, 2, 1)' .* repmat(v/s, s+1, 1) + repmat([AP_x(1) AP_y(1)], s+1, 1);
        p1 = p + repmat(n, s+1, 1) * 170;
        p2 = p - repmat(n, s+1, 1) * 170;

        ii = line([p1(:,1) p2(:,1)]', [p1(:,2) p2(:,2)]', 'lineStyle','--', 'Color',[.8 .8 0]);
        
        for i=1:s+1
            percent = int2str(100/s*(s+1-i));
            text(p2(i,1),p2(i,2),['\color{magenta}', percent, '% EL'], 'horizontalAlignment', 'right');
        end
 
        hold off;
       saveas(gcf, sprintf('AP aixs validation.fig'));
    end
     
end

