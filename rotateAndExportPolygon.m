function rotateAndExportPolygon(init, x, y, filename) 

    % Validate input
    if length(x) ~= length(y)
        error('x and y must be the same length.');
    end
    if nargin < 3
        filename = 'output.stl';
    end

    % Generate rotation symmetry about x-axis
    theta = linspace(0, 2*pi, 100); % Angles for full rotation
    [Theta, R] = meshgrid(theta, y); % Create mesh grid for rotation

    % Compute coordinates for 3D surface
    Z = R .* sin(Theta);          % Z-coordinates (rotated y values)
    Y = R .* cos(Theta);          % Y-coordinates (rotated y values)
    X = repmat(x(:), 1, length(theta)); % Repeat x-coordinates along rows

    % Surface plot
    figure;
    s = surf(X, Y, Z, 'EdgeColor', 'none'); % Create surface
    direction = [0 1 0];  
    rotate(s,direction,-rad2deg(init.theta0))
    colormap('turbo'); % Color map for visualization
    lighting gouraud; % Enhance lighting
    title('3D Surface of Rotated Polygon');
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    axis equal; % Equal scaling for axes
    grid on; % Enable grid

    

    % % Convert surface to triangulation
    % fv = surf2patch(X, Y, Z, 'triangles'); % Convert to face-vertex format
    % 
    % 
    % % Ensure triangulation compatibility
    % vertices = fv.vertices; % Vertex list
    % faces = fv.faces;       % Triangular faces
    % 
    % % Write to STL
    % fprintf('Writing to STL file: %s\n', filename);
    % % writeSTL(filename, vertices, faces);
end
