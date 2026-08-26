function centroids = findCentroids3D(labeledImage)
    % Find unique labels excluding background (assumed to be 0)
    labels = unique(labeledImage);
    labels(labels == 0) = [];
    
    % Initialize an array to hold centroids
    centroids = zeros(length(labels), 4); % [label, x, y, z]
    
    % Loop through each label
    for i = 1:length(labels)
        
        label = labels(i);
        label
        % Find the indices of the current label
        [x, y, z] = ind2sub(size(labeledImage), find(labeledImage == label));
        
        % Calculate the centroid for the current label
        centroid_x = mean(x);
        centroid_y = mean(y);
        centroid_z = mean(z);
        
        % Store the label and its centroid
        centroids(i, :) = [label, centroid_x, centroid_y, centroid_z];
    end
end
