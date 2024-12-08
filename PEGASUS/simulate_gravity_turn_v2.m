function init = simulate_gravity_turn_v2(init)
    % Use a structure for parameters to improve readability and maintenance
    params = struct();

    % Initialize parameters from the `init` structure
    init_vars = {'target_theta','Mo', 'Mp', 'theta0', 'L_rail', 'v0', 'x0', 'y0', 'T_W', 'r_cs', 't_total', ...
                 'seperationCase', 'select_drag_model', 'Ml', 'thrust_off_t', 'thrust_off_h', 'select_thrust_model', ... 
                 'variable_mdot', 'variable_Pe', 'variable_Isp', 'variable_Ae'};
    
    % Assign values to the params structure
    for i = 1:length(init_vars)
        params.(init_vars{i}) = init.(init_vars{i});
    end


    % Extract values from params structure
    Mo = params.Mo;
    Mp = params.Mp;
    theta0 = params.theta0;
    L_rail = params.L_rail;
    v0 = params.v0;
    x0 = params.x0;
    y0 = params.y0;
    r_cs = params.r_cs;
    t_total = params.t_total;
    Ml = params.Ml;
    target_theta = params.target_theta; 

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
                        init.burnout_n = i-1;
                        init.burnout_theta = theta(i-1);
                        init.burnout_v = v_norm(i-1);
                        init.burnout_x = x(i-1);
                        init.burnout_h = y(i-1);
                        init.burnout_Mach = Machs(i-1); 
                        init.burnout_time = t(i-1); 
                        fprintf(strcat(init.Stage_Name, ' Engine Burnout:\n\n'))
                        
                        fprintf('Time: %f s\n', init.burnout_time)
                        fprintf('Altitude: %f km\n', init.burnout_h/1000)
                        fprintf('Velocity: %f m.s\n', init.burnout_v)
                        fprintf('Mach: %f\n', init.burnout_Mach)
                        fprintf('Theta: %f deg\n', rad2deg(init.burnout_theta)) % Sanity Check
                        fprintf('%%Error from 35 Km: %f %%\n',  abs(init.burnout_h/1000 -35)/35 *100)
                        fprintf('%%Error from Mach 10: %f %%\n\n\n\n',  abs(init.burnout_Mach -10)/10 *100)
                    end

                    if ~seperated
                        mass(i) = mass(i-1) - Ml;
                        Ft_current = 0; 
                        seperated = true; 
                        init.seperation_n = i-1; 
                        init.seperation_theta = theta(i-1);
                        init.seperation_v = v_norm(i-1);
                        init.seperation_x = x(i-1);
                        init.seperation_h = y(i-1);
                        init.seperation_Mach = Machs(i-1); 
                        init.seperation_time = t(i-1); 
                        fprintf(strcat(init.Stage_Name,' Payload Launched at Burnout Altitude:\n\n'))
                        
                        fprintf('Time: %f s\n', init.seperation_time)
                        fprintf('Altitude: %f km\n', init.seperation_h/1000)
                        fprintf('Velocity: %f m.s\n', init.seperation_v)
                        fprintf('Mach: %f\n', init.seperation_Mach)
                        fprintf('Theta: %f deg\n', rad2deg(init.seperation_theta)) % Sanity Check
                        fprintf('%%Error from 35 Km: %f %%\n',  abs(init.seperation_h/1000 -35)/35 *100)
                        fprintf('%%Error from Mach 10: %f %%\n\n\n\n',  abs(init.seperation_Mach -10)/10 *100)
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
                        init.burnout_n = i-1;
                        init.burnout_theta = theta(i-1);
                        init.burnout_v = v_norm(i-1);
                        init.burnout_x = x(i-1);
                        init.burnout_h = y(i-1);
                        init.burnout_Mach = Machs(i-1); 
                        init.burnout_time = t(i-1); 
                        fprintf(strcat(init.Stage_Name,' Engine Burnout:\n\n'))
                        
                        fprintf('Time: %f s\n', init.burnout_time)
                        fprintf('Altitude: %f km\n', init.burnout_h/1000)
                        fprintf('Velocity: %f m.s\n', init.burnout_v)
                        fprintf('Mach: %f\n', init.burnout_Mach)
                        fprintf('Theta: %f deg\n', rad2deg(init.burnout_theta)) % Sanity Check
                        fprintf('%%Error from 35 Km: %f %%\n',  abs(init.burnout_h/1000 -35)/35 *100)
                        fprintf('%%Error from Mach 10: %f %%\n\n\n\n',  abs(init.burnout_Mach -10)/10 *100)

                        
                    end
                    
                    if rad2deg(theta(i-1)) <= target_theta && ~seperated 
                        mass(i) = mass(i-1) - Ml;
                        Ft_current = 0; 
                        seperated = true; 
                        init.seperation_n = i-1; 
                        init.seperation_theta = theta(i-1);
                        init.seperation_v = v_norm(i-1);
                        init.seperation_x = x(i-1);
                        init.seperation_h = y(i-1);
                        init.seperation_Mach = Machs(i-1); 
                        init.seperation_time = t(i-1); 
                        fprintf(strcat(init.Stage_Name,' Payload Launched Near Peak Altitude:\n\n'))
                        
                        fprintf('Time: %f s\n', init.seperation_time)
                        fprintf('Altitude: %f km\n', init.seperation_h/1000)
                        fprintf('Velocity: %f m.s\n', init.seperation_v)
                        fprintf('Mach: %f\n', init.seperation_Mach)
                        fprintf('Theta: %f deg\n', rad2deg(init.seperation_theta)) % Sanity Check
                        fprintf('%%Error from 35 Km: %f %%\n',  abs(init.seperation_h/1000 -35)/35 *100)
                        fprintf('%%Error from Mach 10: %f %%\n\n\n\n',  abs(init.seperation_Mach -10)/10 *100)
                        
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
     init.v, init.v_norm, init.a, init.theta, init.t, init.Machs, init.Drags, ...
     init.Thrusts, init.rhos, init.mass, init.ReNums] = deal(x, vx, axs, y, vy, ays, ...
        sqrt(vx.^2 + vy.^2), v_norm, sqrt(axs.^2 + ays.^2), theta, t, Machs, Drags, Thrusts, ...
        rhos, mass, ReNums);
 
    % [init, state_T] = telemetry_results(init);  % Final telemetry and state


    init = reporting(init) ; 
end

