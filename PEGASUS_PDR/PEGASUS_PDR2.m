%% Author: Liam Hunte
%% I, Liam Hunte, affirm that I am the original author of this MATLAB program and that it is completely of my own design.
%% Coursework: AE441 Rocket Preliminary Design

clear; close all; clc; format long g ; 
load init1.mat
load init2.mat


init1.end_at_peak=false; init2.end_at_peak=false;
theta_o= init1.theta0; 
unc_pos = ones(numel(theta_o),3) .* [1.001 1.002 1.003] ; unc_pos = theta_o .*unc_pos ; 
unc_neg = ones(numel(theta_o),3) .* [0.999 0.998 0.997] ; unc_neg = theta_o .*unc_neg ; 
m_c_i = [unc_neg, theta_o, unc_pos]; %Theta confidence inteval

n = size(m_c_i, 2);
O = zeros(size(0:0.001:1000,1),n); 
Ox1 = O; Oy1 = O; Ox2 = O; Oy2 = O;

for ii = 1:n
    init1.theta0 = m_c_i(ii);

    init1 = firstStage(init1); 
    init1 = simulate_gravity_turn_v2(init1);

    len1= numel(init1.x);
    % Ox1 = Ox1(1:len1,:); Oy1 = Oy1(1:len1,:);
    Ox1(1:len1, ii) = init1.x;
    Oy1(1:len1, ii) = init1.y;

    init2 = secondStage(init1, init2) ; 
    init2 = simulate_gravity_turn_v2(init2);

    len2= numel(init2.x); 
    % Ox2 = Ox2(1:len1,:); Oy2 = Oy2(1:len1,:);
    Ox2(1:len2, ii) = init2.x;
    Oy2(1:len2, ii) = init2.y;

end

Ox1 = Ox1([boolean(1) ;Ox1(2:end, 1) ~= 0], :) ;
Oy1 = Oy1([boolean(1) ;Oy1(2:end, 1) ~= 0], :) ; 

Ox2 = Ox2([boolean(1) ;Ox2(2:end, 1) ~= 0], :) ;
Oy2 = Oy2([boolean(1) ;Oy2(2:end, 1) ~= 0], :) ; 

figure 
plot(Ox1(:,4) , Oy1(:,4), LineStyle="-",Color= 'k', LineWidth=2)
hold on
plot(Ox2(:,4) , Oy2(:,4) , LineStyle="-",Color= 'red', LineWidth=2)
hold off
title('Nominal Gravity Turn Trajectory')
subtitle(['Launch Angle: ' '62.8' char(176)])
legend(["Stage 1",  "Stage 2"], Location="south")
axis equal
grid on

figure
plot(Ox1(:,4) , Oy1(:,4), LineStyle="-", Color= 'k', LineWidth=2)
hold on 
plot(Ox1(:,1) , Oy1(:,1), LineStyle="--")
plot(Ox1(:,2) , Oy1(:,2), LineStyle="--")
plot(Ox1(:,3) , Oy1(:,3), LineStyle="--")
plot(Ox1(:,5) , Oy1(:,5), LineStyle="--")
plot(Ox1(:,6) , Oy1(:,6), LineStyle="--")
plot(Ox1(:,7) , Oy1(:,7), LineStyle="--")
hold off
title('1st Stage Normal Trajectories')
subtitle('Launch Angle Uncertainty')
legend(["+/-0%",  "-0.1%", "-0.2%", "-0.3%" ,"+0.1%", "+0.2%", "+0.3%"], Location="northwest")
axis equal
grid on

figure
plot(Ox2(:,4) , Oy2(:,4), LineStyle="-", Color= 'k', LineWidth=2)
hold on 
plot(Ox2(:,1) , Oy2(:,1), LineStyle="--")
plot(Ox2(:,2) , Oy2(:,2), LineStyle="--")
plot(Ox2(:,3) , Oy2(:,3), LineStyle="--")
plot(Ox2(:,5) , Oy2(:,5), LineStyle="--")
plot(Ox2(:,6) , Oy2(:,6), LineStyle="--")
plot(Ox2(:,7) , Oy2(:,7), LineStyle="--")
hold off
title('2nd Stage Normal Trajectories')
subtitle('Launch Angle Uncertainty')
legend(["+/-0%",  "-0.1%", "-0.2%", "-0.3%" ,"+0.1%", "+0.2%", "+0.3%"], Location="northwest")
axis equal
grid on
