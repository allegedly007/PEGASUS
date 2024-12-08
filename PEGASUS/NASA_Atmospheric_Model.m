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
