function init = telemetry_results(init)

    [~, n2] = max(init.y/1000) ;
    init.n2 = n2; 

    init.altitude_at_peak = init.y(n2)/1000 ; 
    init.velocity_at_peak = init.v(n2) ; 
    init.Mach_at_peak = init.Machs(n2) ;
    init.ThetaRadians_at_peak = init.theta(1, n2); 
    init.ThetaDegrees_at_peak = rad2deg(init.ThetaRadians_at_peak) ; 
    init.altitude_pct_err_at_peak = abs(init.y(n2)/1000 -35)/35 *100;
    init.Mach_pct_err_at_peak = abs(init.Machs(n2) -10)/10 *100 ;

    if init.print_results
        fprintf('Peak Altitude: %f km\n', init.altitude_at_peak)
        fprintf('Velocity at Peak: %f m.s\n', init.velocity_at_peak)
        fprintf('Mach at Peak: %f\n', init.Mach_at_peak)
        fprintf('Theta at Peak: %f deg\n', init.ThetaDegrees_at_peak) % Sanity Check
        fprintf('Pct Error from 35 Km: %f %%\n', init.altitude_pct_err_at_peak )
        fprintf('Pct Error from Mach 10: %f %%\n', init.Mach_pct_err_at_peak)
        fprintf('Minimum Acceleration: %fg\n', min(init.a(init.a(1:init.burnout_n) > 0) ) /9.81)
        fprintf('Maximum Acceleration: %fg\n', max(init.a(1:init.burnout_n)) / 9.81) 
        fprintf('Minimum Thrust: %f N\n', min(init.Thrusts(init.Thrusts > 0) ))
        fprintf('Maximum Thrust: %f N\n\n\n\n', max(init.Thrusts)) 

    end
end

