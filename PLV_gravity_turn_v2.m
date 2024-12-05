%% Author: Liam Hunte
%% I, Liam Hunte, affirm that I am the original author of this MATLAB program and that it is completely of my own design.
%% Coursework: AE441 Rocket Preliminary Design

clear; close all; clc; format long g
%% Independent Variables: Mass & Isp
Isp_stage1 = 338.4867;   % Specific impulse (s)
Isp_stage2 = 338.4867;   % Specific impulse (s)

Mp_stage1 = 11000;            % Propellant mass (kg)
Ms_stage1 = 5000;             % Structural mass (kg)
Mp_stage2 = 7500;             % Propellant mass (kg)
Ms_stage2 = 3000;             % Structural mass (kg)
Ml = 1500;                    % Payload mass (kg)

Mo_stage2 = Mp_stage2 + Ms_stage2 + Ml;  % Initial mass of stage 2 (kg)
Mo = Mp_stage1 + Ms_stage1 + Mo_stage2;  % Total initial mass (kg)

% Output total rocket mass
fprintf('Total Rocket Mass: %f kg\n\n\n', Mo)

%% Toggle for desired output (init - stage 1, init2 - stage 2)
init = struct('select_drag_model', 'og', ...
    'plot_nosecone', false, 'plot_rocketbody', false, ...
    'create_gif', false, 'gif_name', 'stage1', 'plot_trajectory', false, ...
    'print_telemetry', true, 'export_data', false, 'seperationCase', 'burnout', 'select_thrust_model', 'hybrid', ...
    'variable_mdot', 583.2219, 'variable_Pe', 4.7184e3, 'variable_Isp', 338.4867, 'variable_Ae', 0.8683^2 *pi, ... 
    'end_at_peak', true);

init2 = struct('select_drag_model', init.select_drag_model, ...
    'plot_nosecone', false, 'plot_rocketbody', false, ...
    'create_gif', false, 'gif_name', 'stage2', 'plot_trajectory', false, ...
    'print_telemetry', true, 'export_data', false, 'seperationCase', 'peak', 'select_thrust_model', 'hybrid', ...
    'variable_mdot', 0, 'variable_Pe', 4.7184e3, 'variable_Isp', 338.4867, 'variable_Ae', 0.8683^2 *pi, ...
    'end_at_peak', true);

% Set stage parameters
init.Mo = Mo; init.Mp = Mp_stage1; init.Isp = Isp_stage1; init.Ms = Ms_stage1;
init.Mo2 = Mo_stage2; init2.Ml = Ml; init2.Mo = Mo_stage2; init2.Mp = Mp_stage2; 
init2.Isp = Isp_stage2; init2.Ms = Ms_stage2;


%% Initial Launch Angle and T/W Ratio
init.theta0 = 62.8* pi / 180;  % Launch angle in radians
init.T_W = 5.5; init2.T_W = 2.25; %Dynamic Engine throttling  

[T, init, init2] = vehicle_dimensions(init, init2);
display(T);
init = plot_rocket_body(init);
init2.body.rocket = init.body.rocketupperstage;


%% Initial Conditions
init.L_rail = 100; init.v0 = 0; init.x0 = 0; init.y0 = 0;
init.r_cs = 1.5; init.t_total = 1000; init.thrust_off_t = NaN; 
init.thrust_off_h = NaN; init.jettison_stage = true;
init.Me_jettison = Ms_stage1; init.Ml = Mo_stage2; init.jettison_peak = false;
init.Me_peak_jettison = NaN; init.theta_peak_margin = NaN; init.Stage_Name = "Stage-1";

[burnout_T1, t, x, y, theta1, init] = simulate_gravity_turn(init);
display(burnout_T1);

%% Hot Stage - Stage 2
init2.theta0 = init.seperation_theta; init2.L_rail = 0; 
init2.v0 = init.seperation_v; init2.x0 = init.seperation_x; 
init2.y0 = init.seperation_h; init2.r_cs = init.r_cs;
init2.t_total = 1000; init2.thrust_off_t = NaN; init2.thrust_off_h = NaN; 
init2.jettison_stage = false; init2.Me_jettison = NaN;
init2.jettison_peak = true; init2.Me_peak_jettison = Ms_stage2;
init2.theta_peak_margin = NaN; init2.Stage_Name = "Stage-2";

[burnout_T2, t2, x2, y2, theta2, init2] = simulate_gravity_turn(init2);
display(burnout_T2);



function [state_T, t, x, y, theta, init] = simulate_gravity_turn(init)
    % Use a structure for parameters to improve readability and maintenance
    params = struct();

    % Initialize parameters from the `init` structure
    init_vars = {'Isp', 'Mo', 'Mp', 'theta0', 'L_rail', 'v0', 'x0', 'y0', 'T_W', 'r_cs', 't_total', ...
                 'seperationCase', 'select_drag_model', 'Ml', 'thrust_off_t', 'thrust_off_h', 'select_thrust_model', ... 
                 'variable_mdot', 'variable_Pe', 'variable_Isp', 'variable_Ae'};
    
    % Assign values to the params structure
    for i = 1:length(init_vars)
        params.(init_vars{i}) = init.(init_vars{i});
    end


    % Extract values from params structure
    Isp = params.Isp;
    Mo = params.Mo;
    Mp = params.Mp;
    theta0 = params.theta0;
    L_rail = params.L_rail;
    v0 = params.v0;
    x0 = params.x0;
    y0 = params.y0;
    T_W = params.T_W;
    r_cs = params.r_cs;
    t_total = params.t_total;
    Ml = params.Ml;

    g0 = 9.81;  % Gravity at Earth's surface (m/s^2)
    

    % Simulation loop
    on_rail = true; % Flag to indicate if the rocket is still on rail 
    seperated = false;  % Flag to indicate if the rocket's upper stage /payload has seperated
    burnout = false; % Flag to indicate if the rocket's engine has reached burnout for the particular stage
    
    seperationCases = ["burnout", "peak", "rand"];
    seperationCase_pos = find(strcmp(seperationCases, params.seperationCase));

    dragmodelCases = ["og", "cfd", "g6"];
    dragmodelCases_pos = find(strcmp(dragmodelCases, params.select_drag_model));

    thrustmodelCases = ["constant", "variable", "hybrid"];
    thrustmodelCases_pos = find(strcmp(thrustmodelCases, params.select_thrust_model)); 


    [mdot, ~] = thrust_model(params, thrustmodelCases_pos, 101325); 
    init.mdot = mdot; 


    A = pi * r_cs^2;         % Cross-sectional area of rocket (m^2)
    R = 6378000;             % Radius of Earth (m)

    % Time settings
    t_step = 0.001;
    burn_time = Mp / mdot;  % Time until fuel depletion (s)


    % Preallocate memory for trajectory data arrays
    t = 0:t_step:t_total;
    O = zeros(size(t));  % Preallocate a single array with the same size as the time vector
    [x, y, vx, vy, theta, mass, Machs, Drags, axs, ays, v_norm, rhos, Thrusts, ReNums] = deal(O);

    % Initial Conditions
    x(1) = x0; y(1) = y0; vx(1) = v0 * cos(theta0); vy(1) = v0 * sin(theta0); 
    theta(1) = theta0; mass(1) = Mo;
    
    

    for i = 2:length(t)
        % Calculate altitude-dependent parameters
        h = y(i-1);                    % Altitude (m)
        [Tatm, Patm, rho] = NASA_Atmospheric_Model(h);
        rhos(i) = rho; 

        [mu, ~] = SutherlandsLaw(Tatm);
        g = g0 * (R / (R + h))^2;  % Gravity with altitude (m/s^2)

        [~, Ft] = thrust_model(params, thrustmodelCases_pos, Patm); 

        % Update mass until burnout
        switch seperationCase_pos
            case 1   % "burnout"
                if t(i) <= burn_time 
                    mass(i) = mass(i-1) - mdot * t_step;
                    Ft_current = Ft;
                else
                    if ~burnout
                        burnout = true; 
                        init.burnout_n = i; 
                    end

                    if ~seperated
                        mass(i) = mass(i-1) - Ml;
                        Ft_current = 0; 
                        seperated = true; 
                        init.seperation_theta = theta(i-1);
                        init.seperation_v = v_norm(i-1);
                        init.seperation_x = x(i-1);
                        init.seperation_h = y(i-1);
                        init.seperation_n = i;
                        fprintf('Payload Launched at Burnout Altitude %.4f km at %.4f deg\n', y(i-1)/1000, rad2deg(theta(i-1)))
                    else
                        mass(i) = mass(i-1);
                        Ft_current = 0; 
                    end
                end
            case 2   % "peak"
                if t(i) <= burn_time 
                    mass(i) = mass(i-1) - mdot * t_step;
                    Ft_current = Ft;
                else
                    if ~burnout
                        burnout = true; 
                        init.burnout_n = i; 
                    end
                    
                    if and(rad2deg(theta(i-1)) <= 0.0005, ~seperated)
                        mass(i) = mass(i-1) - Ml;
                        Ft_current = 0; 
                        seperated = true; 
                        init.seperation_theta = theta(i-1);
                        init.seperation_v = v_norm(i-1);
                        init.seperation_x = x(i-1);
                        init.seperation_h = y(i-1);
                        init.seperation_n = i; 
                        fprintf('Payload Launched at Peak Altitude %.4f km at %.4f deg\n', y(i-1)/1000, rad2deg(theta(i-1)))
                    else
                        mass(i) = mass(i-1);
                        Ft_current = 0; 
                    end
                end
            case 3  % "rand"
                % Implement the random delay case (if needed)
        end

        Thrusts(i) = Ft_current;

        % Calculate velocity magnitude and Mach number
        v = sqrt(vx(i-1)^2 + vy(i-1)^2);  % Total velocity (magnitude)
        a = 340.29 * sqrt(Tatm / 288.15);  % Speed of sound (m/s) at altitude
        Mach = v / a;  % Mach number
        Machs(i) = Mach;

        ReNum = ReVal(rho, v , 16, mu); ReNums(i) = ReNum;
        Cd = drag_model(Mach, dragmodelCases_pos);
        Fd = 0.5 * Cd * rho * A * v^2;  % Drag force (N)
        Drags(i) = Fd;

        % Compute acceleration
        if on_rail
            ax = (Ft_current / mass(i)) * cos(theta0) - (Fd / mass(i)) * cos(theta0);
            ay = (Ft_current / mass(i)) * sin(theta0) - g - (Fd / mass(i)) * sin(theta0);
        else
            ax = (Ft_current / mass(i)) * cos(theta(i-1)) - (Fd / mass(i)) * (vx(i-1) / v);
            ay = (Ft_current / mass(i)) * sin(theta(i-1)) - g - (Fd / mass(i)) * (vy(i-1) / v);
        end

        axs(i) = ax;
        ays(i) = ay;

        % Update velocity and position
        vx(i) = vx(i-1) + ax * t_step;
        vy(i) = vy(i-1) + ay * t_step;
        x(i) = x(i-1) + vx(i) * t_step;
        y(i) = y(i-1) + vy(i) * t_step;

        v_norm(i) = sqrt(vx(i)^2 + vy(i)^2);

        % Check if the rocket has left the rail
        rail_distance = sqrt(x(i)^2 + y(i)^2);  % Total distance traveled along the rail
        if rail_distance >= L_rail
            on_rail = false;  % Rocket has left the rail
        end

        % Update flight angle
        if ~on_rail
            theta(i) = atan2(vy(i), vx(i));
        else
            theta(i) = theta0;
        end

        % Stop if the rocket hits the ground
        if y(i) < 0
            idx = 1:i;
            [x, y, vx, vy, v_norm, axs, ays, theta, t, Machs, Drags, Thrusts, ReNums, rhos, mass] = ...
                deal(x(idx), y(idx), vx(idx), vy(idx), v_norm(idx), axs(idx), ays(idx), theta(idx), t(idx), Machs(idx), Drags(idx), Thrusts(idx), ReNums(idx), rhos(idx), mass(idx));     
            init.impact_time = t(i); 
            break;

        elseif theta(i) < 0 && init.end_at_peak
            idx = 1:i;
            [x, y, vx, vy, v_norm, axs, ays, theta, t, Machs, Drags, Thrusts, ReNums, rhos, mass] = ...
                deal(x(idx), y(idx), vx(idx), vy(idx), v_norm(idx), axs(idx), ays(idx), theta(idx), t(idx), Machs(idx), Drags(idx), Thrusts(idx), ReNums(idx), rhos(idx), mass(idx));     
            init.impact_time = t(i); 
            break;

        end
    end

    % Save results to init structure
    [init.x, init.vx, init.axs, init.y, init.vy, init.ays, ...
     init.v, init.v_norm, init.theta, init.t, init.Machs, init.Drags, ...
     init.Thrusts, init.rhos, init.mass, init.ReNums] = deal(x, vx, axs, y, vy, ays, ...
        sqrt(vx.^2 + vy.^2), v_norm, theta, t, Machs, Drags, Thrusts, ...
        rhos, mass, ReNums);
 
    % [init, state_T] = telemetry_results(init);  % Final telemetry and state


    [init,state_T] = reporting(init) ; 
end



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


function [T, P, rho] = NASA_Atmospheric_Model(h)
    if h < 11000
        T = (15.04 - 0.00649 * h) + 273.15; % Temperature (Kelvin)
        P = 101.29 * (T / 288.08)^5.256; % Pressure (KPa)
    elseif h >= 11000 && h <= 25000
        T = -56.46 + 273.15; % Temperature (Kelvin)
        P = 22.65 * exp(1.73 - 0.000157 * h); % Pressure (KPa)
    else
        T = (-131.21 + 0.002999 * h) + 273.15; % Temperature (Kelvin)
        P = 2.4888 * (T / 216.6)^(-11.388); % Pressure (KPa)
    end

    rho = P / (0.2869 * T);
    P = P * 1000; 
end


function [mu, k] = SutherlandsLaw(Tatm)
    T = Tatm; 
    mu_o = 1.716e-5; %Air viscosity 
    To = 273;
    S_mu = 111;
    mu = mu_o * (T/To) ^(3/2) *(To+S_mu) / (T+S_mu) ; 

    ko = 0.0241; 
    S_k = 194; 
    k = ko * (T/To) ^(3/2) *(To+S_k) / (T+S_k) ; 
end


function Re = ReVal(rho, u , L, mu)
    Re = (rho*u*L) / mu;
end

function Cd = drag_model(Mach, dragmodelCases_pos)
    if dragmodelCases_pos == 1
        if Mach < 0.8
            Cd = 0.1; 
        elseif Mach < 1.2
            Cd = 0.3;
        else
            Cd = 0.2; 
        end
    elseif dragmodelCases_pos == 2
        % Add logic for case 2 if necessary
    elseif dragmodelCases_pos == 3
        if Mach <= 1
            % Handle case when Mach <= 1
        else
            % Handle case when Mach > 1 if needed
        end
    end
end


function [init,state_T] = telemetry_results(init)

    [~,n] = max(init.v);
    init.n1 = n; 

    init.theta_at_burnout = rad2deg(init.theta(1, n)) ; 
        
    init.time_at_burnout = init.t(n) ; 
    
    init.altitude_at_burnout = init.y(n) /1000 ;
    
    init.distance_at_burnout = init.x(n) /1000 ; 
        
    init.velocity_at_burnout = init.v(n) ; 

    init.Mach_at_burnout = init.Machs(n) ; 
    
    [~, n2] = max(init.y/1000) ;
    init.n2 = n2; 

    init.altitude_at_peak = init.y(n2)/1000 ; 
    init.velocity_at_peak = init.v(n2) ; 
    init.Mach_at_peak = init.Machs(n2) ;
    init.ThetaRadians_at_peak = init.theta(1, n2); 
    init.ThetaDegrees_at_peak = rad2deg(init.ThetaRadians_at_peak) ; 
    init.altitude_pct_err_at_peak = abs(init.y(n2)/1000 -35)/35 *100;
    init.Mach_pct_err_at_peak = abs(init.Machs(n2) -10)/10 *100 ;

    state = [init.theta_at_burnout, init.velocity_at_burnout, init.distance_at_burnout*1000, init.altitude_at_burnout*1000, init.time_at_burnout, init.Mach_at_burnout]; 
    state_T = array2table(state); 

    state_T.Properties.VariableNames = {'theta_b', 'v_b', 'x_b', 'h_b', 't_b', 'M_b'};

    if init.print_telemetry

        fprintf('Peak Altitude: %f km\n', init.altitude_at_peak)
        fprintf('Velocity at Peak: %f m.s\n', init.velocity_at_peak)
        fprintf('Mach at Peak: %f\n', init.Mach_at_peak)
        fprintf('Theta at Peak: %f deg\n', init.ThetaDegrees_at_peak) % Sanity Check
        fprintf('Pct Error from 35 Km: %f %%\n', init.altitude_pct_err_at_peak )
        fprintf('Pct Error from Mach 10: %f %%\n', init.Mach_pct_err_at_peak)

    end
end



function [T, init,init2] = vehicle_dimensions(init, init2)
    go = 9.81;
    T = array2table(zeros(5, 2), "VariableNames", {'Stage1', 'Stage2'}, "RowNames", {'Isp', 'Ml', 'Ms', 'Mp', 'Mo'}); 
    T.Stage2 = [init2.Isp, init2.Ml, init2.Ms, init2.Mp, 0]' ; T{"Mo", "Stage2"} = sum(T.Stage2(2:end)) ;
    T.Stage1 = [init.Isp, T{"Mo", "Stage2"}, init.Ms, init.Mp, 0]' ; T{"Mo", "Stage1"} = sum(T.Stage1(2:end)) ;

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


function [init,state_T] = reporting(init)
    dataframe = [init.t', init.x', init.y', init.vx', init.vy', init.axs', init.ays', init.theta', init.Machs' , init.Drags' , init.Thrusts', init.mass', init.rhos', init.ReNums'] ;
    
    plotTrajectory(init);

    [init,state_T] = telemetry_results(init) ; 

    init.dataframe = array2table(dataframe) ; 
    init.dataframe.Properties.VariableNames = {'t', 'x', 'y', 'vx' , 'vy', 'axs', 'ays', 'theta', 'Mach', 'Drag', 'Thrust', 'mass', 'rho', 'Re'};
    
    filename_precursor1 = "AE441_Pegasus_LaunchTOBurnout_" + init.Stage_Name + "_Trajectory_SI_" ; 
    exportData(init.dataframe(1:init.n1, :), filename_precursor1, init.export_data)

    filename_precursor2 = "AE441_Pegasus_BurnoutTOPeak_" + init.Stage_Name + "_Trajectory_SI_" ; 
    exportData(init.dataframe(init.n1:init.n2, :), filename_precursor2, init.export_data)

    filename_precursor3 = "AE441_Pegasus_PeakTOImpact_" + init.Stage_Name +"_Trajectory_SI_" ; 
    exportData(init.dataframe(init.n2:end, :), filename_precursor3, init.export_data)

    filename_precursor4 = "Ae441_Pegasus_Parameters" ; 
    exportData(init.T_PLV, filename_precursor4, init.export_data)


end



function plotTrajectory(init)
    if init.plot_trajectory
        %Plotting the trajectory
        figure
        plot(init.x/1000, init.y/1000);
        hold on 
        plot(init.x(init.seperation_n)/1000, init.y(init.seperation_n)/1000, "diamond",'MarkerFaceColor','green');
        plot(init.x(init.burnout_n)/1000, init.y(init.burnout_n)/1000, '-p','MarkerFaceColor','magenta');
        hold off
        xlabel('Horizontal Distance (km)');
        ylabel('Vertical Distance (km)');
        title('Gravity Turn Trajectory');
        subtitle(init.Stage_Name)
        grid on;
        legend(['Flight Dynamics', "Stage Seperation", "Burnout"], Location="south")

        figure 
        plot(init.t, init.y/1000)
        hold on 
        plot(init.t(init.seperation_n), init.y(init.seperation_n)/1000, "diamond",'MarkerFaceColor','green');
        plot(init.t(init.burnout_n), init.y(init.burnout_n)/1000, '-p','MarkerFaceColor','magenta');
        hold off
        xlabel('Time (s)');
        ylabel('Altitude (km)');
        title('Altitude v. Time');
        subtitle(init.Stage_Name)
        grid on;
        legend(['Flight Dynamics', "Stage Seperation", "Burnout"], Location="south")


        figure 
        plot(init.t, init.v)
        hold on 
        plot(init.t(init.seperation_n), init.v(init.seperation_n), "diamond",'MarkerFaceColor','green');
        plot(init.t(init.burnout_n), init.v(init.burnout_n), '-p','MarkerFaceColor','magenta');
        hold off
        xlabel('Time (s)');
        ylabel('Velocity (m.s)');
        title('Velocity v. Time');
        subtitle(init.Stage_Name)
        grid on;
        legend(['Flight Dynamics', "Stage Seperation", "Burnout"], Location="south")

        figure 
        plot(init.t, rad2deg(init.theta))
        hold on 
        plot(init.t(init.seperation_n), rad2deg(init.theta(init.seperation_n)), "diamond",'MarkerFaceColor','green');
        plot(init.t(init.burnout_n), rad2deg(init.theta(init.burnout_n)), '-p','MarkerFaceColor','magenta');
        hold off
        xlabel('Time (s)');
        ylabel('Theta (deg)');
        title('Theta v. Time');
        subtitle(init.Stage_Name)
        grid on;
        legend(['Flight Dynamics', "Stage Seperation", "Burnout"], Location="south")


        % figure 
        % plot(init.t, init.ays)
        % hold on 
        % plot(init.t(init.seperation_n), init.ays(init.seperation_n), "diamond",'MarkerFaceColor','green');
        % plot(init.t(init.burnout_n), init.ays(init.burnout_n), '-p','MarkerFaceColor','magenta');
        % hold off
        % xlabel('Time (s)');
        % ylabel('Ay (m/s^2)');
        % title('Ay v. Time');
        % subtitle(init.Stage_Name)
        % grid on;
        % legend(['Flight Dynamics', "Stage Seperation", "Burnout"], Location="south")


        % figure
        % plot(init.Drags, init.ReNums)
        % hold on 
        % plot(init.Drags(init.seperation_n), init.ReNums(init.seperation_n), "diamond",'MarkerFaceColor','green');
        % plot(init.Drags(init.burnout_n), init.ReNums(init.burnout_n), '-p','MarkerFaceColor','magenta');
        % hold off
        % xlabel('Re');
        % ylabel('Drag (N)');
        % title('Drag v. Re');
        % subtitle(init.Stage_Name)
        % grid on;
        % legend(['Flight Dynamics', "Stage Seperation", "Burnout"], Location="south")



        

    end
end


function exportData(data, filename_precursor, state)
    if state
        date_time = string(datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm-ss')) ; 
        
        filename = filename_precursor + date_time + ".xlsx"  ;  
        fprintf('%s\n', filename)
        
        writetable(data, filename,'Sheet',1, 'WriteRowNames',true)
    end
end


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

    Stage1 = []; 

    max(xStage1)

    Stage2 = [] ; 

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

    
    reducedX = reducepoly(X);
    rocket =reducedX';

    [rocket(2,:),sort_idx_o] = sort(rocket(2,:)) ; 
    rocket(1,:) = rocket(1,sort_idx_o); 
    

    init.body.rocket = rocket;

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

            k = boundary(C(:,1), C(:,2)); 
            val = [C(k,1), C(k,2)];
            figure
            plot(C(k,1), C(k,2)); 
            axis equal
            title('PLV Boundary')
            
            rotateAndExportPolygon(init, val(:,1), val(:,2), 'testpolygon.stl')
            
            coordinates = [val, zeros(size(val, 1), 1)] ; size(coordinates);
            
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



function rotatedPoints = rotate2D(points, theta)

    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    rotatedPoints = zeros(2,size(points , 2)) ; 
    for ii = 1:size(points , 2)
                rotatedPoints(:, ii) = R*points(:, ii); 
    end

    
end



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