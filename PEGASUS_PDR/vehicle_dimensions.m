function [T, init,init2] = vehicle_dimensions(init, init2)
    go = 9.81;
    T = array2table(zeros(5, 2), "VariableNames", {'Stage1', 'Stage2'}, "RowNames", {'Isp', 'Ml', 'Ms', 'Mp', 'Mo'}); 
    T.Stage2 = [init2.variable_Isp, init2.Ml, init2.Ms, init2.Mp, 0]' ; T{"Mo", "Stage2"} = sum(T.Stage2(2:end)) ;
    T.Stage1 = [init.variable_Isp, T{"Mo", "Stage2"}, init.Ms, init.Mp, 0]' ; T{"Mo", "Stage1"} = sum(T.Stage1(2:end)) ;

    T{"MLambda", :} = round(T{"Ml",:} ./ (T{"Mp",:} + T{"Ms",:}) ...
    ,4) ;

    T{"MEpsilon", :} = round(T{"Ms",:} ./ (T{"Mp",:} + T{"Ms",:}) ...
        ,4) ;
    
    T{"MR", :} = round(  (1 + T{"MLambda",:}) ./ (T{"MEpsilon",:} + T{"MLambda",:}) ...
        ,4) ;
    
    T{"DeltaV", :} = round( (T{"Isp",:} *go) .* log(T{"MR", :} ) ...
        ,4) ;

    h_engine = 4.374; 
    h_fairing = 5.625; 
    init.r_fuselage = 1.5; 
    A_cs_fuselage = round(pi*init.r_fuselage^2 ,4) ; 

    rho_methane = 446.0; %@ 100 bar
    rho_LOX = 1141; 

    T{"O/F", :} = [3.5 3.5] ;

    T{"Mfuel", :}   = 1 ./ (T{"O/F", :}  +1) .* T{"Mp", :} ;
    
    T{"Mox", :} = T{"Mp", :} - T{"Mfuel", :} ;
    
    T{"Vol_fuel", :}   = T{"Mfuel", :} / rho_methane;
    
    T{"Vol_ox", :} = T{"Mox", :}  / rho_LOX ;
    
    T{"Hfairing", :} = [0 h_fairing] ;
    
    T{"Hfueltank", :} = round( T{"Vol_fuel", :} ./ A_cs_fuselage ...
        , 4) ;
    T{"Hoxtank", :} = round(T{"Vol_ox", :} ./ A_cs_fuselage ...
        , 4) ;
    
    T{"Hengine", :} = h_engine ; 
    T{"Hstage", :} = [T{"Hfueltank", :} + T{"Hoxtank", :} + h_engine] ;
    T{"Hstage", "Stage2"} = T{"Hstage", "Stage2"} + h_fairing ; 

    T{"Hstage_o", :} = [T{"Hfueltank", :} + T{"Hoxtank", :} + h_engine];
    
    T{"Area_stage", :} = A_cs_fuselage ; 
    
    total_vehicle_height_SI = sum(T{"Hstage", :}) 

    total_vehicle_height_USC = total_vehicle_height_SI *3.28084 ;

    T.LV = zeros([size(T,1) 1] ) ; 

    T{"Total_LV_Height", :} = [0 0 total_vehicle_height_SI] ;

    

    init.T_PLV = T;
    init2.T_PLV = T;


end

