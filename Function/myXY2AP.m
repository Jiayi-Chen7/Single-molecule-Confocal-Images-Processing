function AP_n = myXY2AP(nx, ny, AP_x, AP_y)
    % transform point from XY coordinates to AP coordinates
	AP_vec = [AP_x(1)-AP_x(2); AP_y(1)-AP_y(2)];
    AP_transform = planerot(AP_vec)';
    xy = [nx-AP_x(2), ny-AP_y(2)];
    AP_n = xy * AP_transform;
end