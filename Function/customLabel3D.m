function labeledImage = customLabel3D(binaryImage)
    % Check if the input is a 3D binary image
    if ndims(binaryImage) ~= 3 || ~islogical(binaryImage)
        error('Input must be a 3D binary image.');
    end

    % Get the size of the input image
    [rows, cols, slices] = size(binaryImage);

    % Initialize the labeled image
    labeledImage = zeros(size(binaryImage));
    
    % Initialize the label counter
    currentLabel = 0;
    
    % Define the 26-connectivity neighborhood
    [dx, dy, dz] = ndgrid(-1:1, -1:1, -1:1);
    neighborhoodOffsets = [dx(:), dy(:), dz(:)];
    neighborhoodOffsets(14, :) = []; % remove the (0,0,0) offset

    % Function to check if a voxel is within the image bounds
    function inBounds = isInBounds(x, y, z)
        inBounds = (x > 0) && (x <= rows) && ...
                   (y > 0) && (y <= cols) && ...
                   (z > 0) && (z <= slices);
    end

    % DFS function to label the connected components
    function dfs(x, y, z)
        stack = [x, y, z];
        while ~isempty(stack)
            [cx, cy, cz] = deal(stack(end, 1), stack(end, 2), stack(end, 3));
            stack(end, :) = [];
            for k = 1:size(neighborhoodOffsets, 1)
                nx = cx + neighborhoodOffsets(k, 1);
                ny = cy + neighborhoodOffsets(k, 2);
                nz = cz + neighborhoodOffsets(k, 3);
                if isInBounds(nx, ny, nz) && binaryImage(nx, ny, nz) && ~labeledImage(nx, ny, nz)
                    labeledImage(nx, ny, nz) = currentLabel;
                    stack = [stack; nx, ny, nz];
                end
            end
        end
    end

    % Loop through each voxel in the image
    for i = 1:rows
        for j = 1:cols
            for k = 1:slices
                if binaryImage(i, j, k) && ~labeledImage(i, j, k)
                    currentLabel = currentLabel + 1;
                    labeledImage(i, j, k) = currentLabel;
                    dfs(i, j, k);
                end
            end
        end
    end
end