
function init = Von_Karman(init)
    % L = init.T_PLV{"Hfairing","Stage2"} ;
    L = 6; 
    Rb = init.r_fuselage;

    init.nR = 4; init.xc = Rb/2;

    C = 0; %Von Karman Ogive
    nPoints = 1000;
    resolution = 3;
    x_ogive = round( linspace(0, L, nPoints) ,resolution) ; 
    theta = acos(1-2/L.*x_ogive) ; 
    y_ogive = round((Rb.*sqrt(theta - sin(2.*theta)/2 + C.*sin(theta).^3) ) / sqrt(pi) ,resolution) ; 
    C_N = 2;
    
    rho = (Rb^2 + L^2) / (2*Rb) ;
    A = pi*Rb^2 ;
    V = pi * (L*rho^2 - L^3 /3 - ((rho-Rb)*rho^2)*asin(L/rho)) ; 

    nR = init.nR; 
    xc = init.xc; 

    theta = linspace(0, pi, nPoints);
    x_sphere = round(xc - Rb/nR * cos(theta) ,resolution);  % x-coordinates of the spherical section
    y_sphere = round(Rb/nR * sin(theta) ,resolution);       % y-coordinates of the spherical section

    switch init.plot_nosecone
        case 1
            figure
            plot(x_ogive, y_ogive, 'r', x_sphere, y_sphere, 'b', x_ogive, -y_ogive, 'r', x_sphere, -y_sphere, 'b')
            axis equal
            title('Tangent Ogive Nose Cone Cross Section')
    end

    [x0,y0,iout,jout] = intersections(x_sphere, y_sphere, x_ogive, y_ogive, 1) ; 

    pos_integers_sphere = (round(iout(:)) == (iout(:))) ; 
    
    pos_integers_ogive = (round(jout(:)) == (jout(:))) ; 

    pos_integers_intersection = (pos_integers_sphere == pos_integers_ogive & pos_integers_sphere == 1 & pos_integers_ogive == 1) ;

    x_intersect = x0(pos_integers_intersection) ;
    y_intersect = y0(pos_integers_intersection) ;


    % for i = 1:numel(x_intersect)-1
    %     tangent_point = i; 
    %     blunt_ogive_x =x_ogive(x_ogive >= x_intersect(end-tangent_point)) ;
    %     blunt_ogive_y = y_ogive(x_ogive >= x_intersect(end-tangent_point));
    %     blunt_sphere_x = x_sphere(x_sphere <= x_intersect(end-tangent_point)); 
    %     blunt_sphere_y = y_sphere(x_sphere <= x_intersect(end-tangent_point)) ;
    % 
    %     figure
    %     plot(-blunt_ogive_x, blunt_ogive_y, 'r', -blunt_sphere_x, blunt_sphere_y, 'b')
    %     hold on 
    %     plot(-blunt_ogive_x, -blunt_ogive_y, 'r', -blunt_sphere_x, -blunt_sphere_y, 'b')
    %     hold off
    %     axis equal
    % end

    blunt_ogive_x = -x_ogive(x_ogive >= x_intersect(1)) ;
    blunt_ogive_y = y_ogive(x_ogive >= x_intersect(1));
    blunt_sphere_x = -x_sphere(x_sphere <= x_intersect(1)); 
    blunt_sphere_y = y_sphere(x_sphere <= x_intersect(1)) ;

    shiftx = -min(blunt_ogive_x);
    
    init.cone.C_N = C_N; 
    init.cone.X_n_low = 0.5*L; 
    init.cone.X_n_high = V/A ;
    

    init.cone.x = [blunt_ogive_x, blunt_sphere_x] + shiftx;
    init.cone.y = [blunt_ogive_y, blunt_sphere_y] ;

    [init.cone.x, sort_idx] = sort(init.cone.x);
    init.cone.y = init.cone.y(sort_idx);
    
    conex = [init.cone.x, init.cone.x]; 
    coney = [init.cone.y, -init.cone.y]; 
    [coney, sort_idx] = sort(coney);
    conex = conex(sort_idx);
            
    points = [conex; coney]; 

    init.cone.points = points; 

    init.cone.L = max(points(1,:)); 


    switch init.plot_nosecone
        case 1
            figure 
            plot(points(1,:), points(2,:))
            axis equal
            title('Blunt Tangent Ogive Nose Cone @Rb')
    end

    switch init.plot_nosecone
        case 1
            theta = 90; 
            
            rotatedPoints = rotate2D(points, theta) ;
            
            figure
            plot(rotatedPoints(1,:), rotatedPoints(2,:))
            title('Nose Cone: Upright Orientation')
            axis equal

            init.cone.rotatedPoints = rotatedPoints;

    end


    %A blunt tip reduces the complexity of manufacturing the nose cone. 

    %https://semanticscholar.org/paper/A-Chief-Engineer's-View-of-the-NASA-X-43A-Scramjet-Marshall-Corpening/754f226edb2e38875b0f0eba8ec890e00b1d5ac5/figure/2


end

