function init = plot_rocket_body(init)
    HStage1 = init.T_PLV{"Hstage","Stage1"} ;
    HStage2 = init.T_PLV{"Hstage","Stage2"} - init.T_PLV{"Hfairing","Stage2"} ;

    xStage1 = linspace(0, HStage1, 100);
    xStage2 = linspace(HStage1, HStage1+HStage2, 100); 

    init = Von_Karman(init); 
    cone = init.cone.points;
    cone(1,:) = cone(1,:) + max(xStage2); 

    Stages = [xStage1, xStage2(2:end)] ; 

    Stages = [Stages; ones(1,numel(Stages)) * init.r_fuselage]; 

    Stage1 = [xStage1; ones(1,numel(xStage1)) * init.r_fuselage]; 

    Stage2 = [xStage2; ones(1,numel(xStage2)) * init.r_fuselage]; 

    tail = [zeros(1,100); linspace(-init.r_fuselage, init.r_fuselage,100)];

    InterstageLine = [ones(1,100) * max(xStage1) ; linspace(-init.r_fuselage, init.r_fuselage,100)] ; 
    PayloadLine = [ones(1,100)* max(xStage2); linspace(-init.r_fuselage, init.r_fuselage,100)] ;
    

    rocket = [tail, Stages, [Stages(1,:); -Stages(2,:)],cone];
    [rocket(2,:),sort_idx_o] = sort(rocket(2,:)) ; 
    rocket(1,:) = rocket(1,sort_idx_o); 

    coordinates = rocket'; 
    C = coordinates; 
    pgon = polyshape(C(:,1), C(:,2)); 
    vertices = pgon.Vertices;
    X = vertices;
    X(isnan(X(:,1)),:) = []; 

    k = boundary(C(:,1), C(:,2)); 
    init.body.full_boundary = [C(k,1)' ; C(k,2)']; 


    coordinates_stage1 = [tail,  Stage1,  [Stage1(1,:); -Stage1(2,:)],  InterstageLine];
    coordinates_stage2 =  [InterstageLine,  Stage2, [Stage2(1,:); -Stage2(2,:)],  PayloadLine];
    coordinates_stage2wfairing = [InterstageLine,  Stage2, [Stage2(1,:); -Stage2(2,:)], cone];
    
    %
    [coordinates_stage1(2,:),sort_idx_o] = sort(coordinates_stage1(2,:)) ; 
    coordinates_stage1(1,:) = coordinates_stage1(1,sort_idx_o); 

    coordinates1 = coordinates_stage1'; 
    C1 = coordinates1; 
    pgon = polyshape(C1(:,1), C1(:,2)); 
    vertices = pgon.Vertices;
    X = vertices;
    X(isnan(X(:,1)),:) = []; 

    k = boundary(C1(:,1), C1(:,2)); 
    init.body.stage1_boundary = [C1(k,1)' ; C1(k,2)']; 
    %
    [coordinates_stage2wfairing(2,:),sort_idx_o] = sort(coordinates_stage2wfairing(2,:)) ; 
    coordinates_stage2wfairing(1,:) = coordinates_stage2wfairing(1,sort_idx_o); 

    coordinates2 = coordinates_stage2wfairing'; 
    C2 = coordinates2; 
    pgon = polyshape(C2(:,1), C2(:,2)); 
    vertices = pgon.Vertices;
    X = vertices;
    X(isnan(X(:,1)),:) = []; 

    k = boundary(C2(:,1), C2(:,2)); 
    init.body.stage2wfairing_boundary = [C2(k,1)' ; C2(k,2)']; 
    
    reducedX = reducepoly(X);
    rocket =reducedX';

    [rocket(2,:),sort_idx_o] = sort(rocket(2,:)) ; 
    rocket(1,:) = rocket(1,sort_idx_o); 

    init.body.stage1_fuelline = InterstageLine; 
    init.body.stage2_fuelline = PayloadLine; 
    
    init.body.rocket2Dfilled = rocket;

    rocketupperstage = [InterstageLine, [xStage2(2:end) ;ones(1,numel(xStage2(2:end))) * init.r_fuselage], [xStage2(2:end); -ones(1,numel(xStage2(2:end))) * init.r_fuselage],cone];
    [rocketupperstage(2,:),sort_idx_o] = sort(rocketupperstage(2,:)) ; 
    rocketupperstage(1,:) = rocketupperstage(1,sort_idx_o); 
    init.body.rocketupperstage = rocketupperstage; 

    rocket90 = rotate2D(rocket, 90) ; 
    InterstageLine90 = rotate2D(InterstageLine  , 90) ;
    PayloadLine90 = rotate2D(PayloadLine, 90) ;


    rocketTheta0 = rotate2D(rocket, rad2deg(init.theta0)) ; 
    InterstageLineTheta0 = rotate2D(InterstageLine  , rad2deg(init.theta0)) ;
    PayloadLineTheta0 = rotate2D(PayloadLine, rad2deg(init.theta0)) ;



    switch init.plot_rocketbody
        case 1

            figure 
            plot(rocket(1,:) ,rocket(2,:) )
            hold on 
            plot(InterstageLine(1,:) ,InterstageLine(2,:) , 'g')
            plot(PayloadLine(1,:) ,PayloadLine(2,:) , 'g')
            hold off
            axis equal
            xlim([-1 20])
            title("Pegasus Launch Vehicle 2D Rigid Body")
            subtitle("CFD Orientation")
        
            figure 
            plot(rocket90(1,:) ,rocket90(2,:) )
            hold on 
            plot(InterstageLine90(1,:) ,InterstageLine90(2,:) , 'g')
            plot(PayloadLine90(1,:) ,PayloadLine90(2,:) , 'g')
            hold off
            axis equal
            ylim([0, 20])
            title("Pegasus Launch Vehicle 2D Rigid Body")

            figure 
            plot(rocketTheta0(1,:) ,rocketTheta0(2,:) )
            hold on 
            plot(InterstageLineTheta0(1,:) ,InterstageLineTheta0(2,:) , 'g')
            plot(PayloadLineTheta0(1,:) ,PayloadLineTheta0(2,:) , 'g')
            hold off
            axis equal
            title("Launch Rail Orientation")
            ylim([-1, 20])

            figure
            plot(C(k,1), C(k,2)); 
            axis equal
            title('PLV Boundary')

            
            rotateAndExportPolygon(init, C(k,1), C(k,2), 'testpolygon.stl')
            
            coordinates = [init.body.full_boundary', zeros(size(C(k,1), 1), 1)] ; size(coordinates);
            
            varnames = ["x" , "y", "z"]; 
            numaxes = 1:3;
            polygonT = array2table(coordinates(:,numaxes), "VariableNames",varnames(numaxes)) ;
            filename = "C:\Users\franc\Downloads\plvPolygon.txt"; 
            writetable(polygonT, filename, 'Delimiter',' ','WriteVariableNames', false)

            save('plvPolygon.mat', 'coordinates')

            
            % figure
        
            % txt = strcat('Total Height:', num2str(max(cone_x)), 'm');
            % txt2 = strcat('Maximum Diameter:', num2str(init.r_fuselage*2), 'm');
            % txt3 = strcat('Stage 1 Length:', num2str(HStage1), 'm');
            % txt4 = strcat('Stage 2 Length:', num2str(HStage2 ), 'm');
            % txt5 = strcat('Fairing Length:', num2str(init.cone.L), 'm');
            % title('Pegasus Launch Vehicle 2D-Rigid Body')
            % 
            % pos_txt_x = max(cone_x)/ 2;
            % 
            % t = text(pos_txt_x ,6,txt,"HorizontalAlignment","center");
            % t = text(pos_txt_x , 5,txt2,"HorizontalAlignment","center");
            % t = text(pos_txt_x, -5,txt3,"HorizontalAlignment","right");
            % t = text(pos_txt_x , -4,txt4,"HorizontalAlignment","center");
            % t = text(pos_txt_x , -3,txt5,"HorizontalAlignment","left");
    end
    
    
    
end
