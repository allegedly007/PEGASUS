%% Author: Liam Hunte
%% I, Liam Hunte, affirm that I am the original author of this MATLAB program and that it is completely of my own design.
%% Coursework: AE441 Rocket Preliminary Design
clear; close all; clc; format long g ; 
load init1.mat
load init2.mat

init1.end_at_peak = false;init2.end_at_peak = false; 
init1.create_gif = true; init2.create_gif = true; 

[~, init1,init2] = vehicle_dimensions(init1, init2) ; 
init1 = plot_rocket_body(init1) ; 
init2.body.stored = init1.body.stage2wfairing_boundary; 

init1 = firstStage(init1); 
init1 = simulate_gravity_turn_v2(init1);

init2 = secondStage(init1, init2) ; 
init2 = simulate_gravity_turn_v2(init2) ;

init1.body.stored = init1.body.full_boundary; 

init1.gif_name = "Stage1Trajectory" ; 
init2.gif_name = "Stage2Trajectory" ; 

init1.body.gif.fuelline = init1.body.stage1_fuelline ; 
init2.body.gif.fuelline = init1.body.stage2_fuelline ;



[~, ~] = gif_generator_2D(init1, 1, "Stage1");

[~, ~] = gif_generator_2D(init2, 0, "Stage2");

