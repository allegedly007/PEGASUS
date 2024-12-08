function [mdot_throttle, Ft] = thrust_model(params, thrustmodelCases_pos, Patm)
    g0 = 9.81;  % Gravity at Earth's surface (m/s^2)
    Isp = params.variable_Isp; 
    Mo =params.Mo;  
    T_W = params.T_W; 

    variable_mdot = params.variable_mdot; 
    variable_Pe = params.variable_Pe ;
    variable_Ae = params.variable_Ae ;

    % Constant Thrust 
    W = Mo*g0; 
    T = W*T_W; 
    mdot = T/(Isp*g0); % Mass flow rate (kg/s

 

    switch thrustmodelCases_pos
        case 1   
            mdot_throttle = mdot; 
            % Constant Thrust 
            Ft = Isp * mdot_throttle * g0;    % Thrust force (N)
            
        case 2 
            % Variable Thrust 
            mdot_throttle = variable_mdot; 
            Ft = Isp *g0 * mdot_throttle + variable_Ae*(variable_Pe - Patm);

        case 3 
            mdot_throttle =  mdot; 
            Ft = Isp*g0 * mdot_throttle + variable_Ae*(variable_Pe - Patm);

    end

end
