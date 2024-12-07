function rotatedPoints = rotate2D(points, theta)

    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    rotatedPoints = zeros(2,size(points , 2)) ; 
    for ii = 1:size(points , 2)
                rotatedPoints(:, ii) = R*points(:, ii); 
    end

    
end

