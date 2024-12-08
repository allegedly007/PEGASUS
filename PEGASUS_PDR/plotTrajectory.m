
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

