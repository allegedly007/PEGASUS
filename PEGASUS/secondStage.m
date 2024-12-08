function init2 = secondStage(init, init2)
init2.theta0 = init.seperation_theta; init2.L_rail = 0; 
init2.v0 = init.seperation_v; init2.x0 = init.seperation_x; 
init2.y0 = init.seperation_h; init2.r_cs = init.r_cs;
init2.t_total = 1000; init2.thrust_off_t = NaN; init2.thrust_off_h = NaN; 
init2.jettison_stage = false; 
init2.jettison_peak = true;
init2.theta_peak_margin = NaN; init2.Stage_Name = "Stage-2";
end

