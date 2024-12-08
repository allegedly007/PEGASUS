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
