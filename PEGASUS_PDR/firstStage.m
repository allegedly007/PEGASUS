function init = firstStage(init)
init.L_rail = 100; init.v0 = 0; init.x0 = 0; init.y0 = 0;
init.r_cs = 1.5; init.t_total = 1000; init.thrust_off_t = NaN; 
init.thrust_off_h = NaN; init.jettison_stage = true;
init.jettison_peak = false;
init.Me_peak_jettison = NaN; init.theta_peak_margin = NaN; init.Stage_Name = "Stage-1";
end
