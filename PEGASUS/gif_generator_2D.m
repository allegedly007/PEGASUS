function [bb, fb] = gif_generator_2D(init, sep, stage)
    switch init.create_gif 

        case 1
            figure
            p = plot(nan,nan);
            axis equal
            xlim([-2 20])
            pos_txt_x = 0;
            gifnm = strcat(init.gif_name, '.gif'); 
            set(gca,'XTickLabel',[]);
            set(gca,'YTickLabel',[]);
            thetaList = rad2deg(init.dataframe{:, "theta"}); 
            yList = init.dataframe{:, "y"};
            MachList = init.dataframe{:, "Mach"};
            n = numel(thetaList) ; 

            tankVol = init.T_PLV{"Vol_fuel", stage} + init.T_PLV{"Vol_ox", stage} ; 
            areaStage = init.T_PLV{"Area_stage", stage} ; 
            fuelline_o = init.body.gif.fuelline; 

            h_tank = init.T_PLV{"Hfueltank", stage} + init.T_PLV{"Hoxtank", stage} ; 
            %% Apply Rocket Telemetry Text
            % txt = strcat('Altitude:', num2str(yList(1)), 'm');
            % t1 = text(pos_txt_x ,0,txt,"HorizontalAlignment","center");
            % txt = strcat('Mach:', num2str(MachList(1)));
            % t2 = text(pos_txt_x ,-1,txt,"HorizontalAlignment","center");
            % txt = strcat('Theta:', num2str(thetaList(1)));
            % t3 = text(pos_txt_x ,-2,txt,"HorizontalAlignment","center");
            for i = 1:1000:n 
                theta = thetaList(i); 
                y = yList(i) /1000;
                Mach  = MachList(i);
                
                %% Apply Rocket Telemetry Text
                % delete(t1); 
                % delete(t2); 
                % delete(t3);
                % delete(findall(gcf,'type','annotation'))
                % txt = strcat('Altitude(km):', num2str(y));
                % t1 = text(pos_txt_x, 0 - i/50000, txt,"HorizontalAlignment","left",VerticalAlignment="baseline", Color='r');
                % txt = strcat('Mach:', num2str(Mach));
                % t2 = text(pos_txt_x, -1 - i/50000, txt,"HorizontalAlignment","left",VerticalAlignment="baseline", Color='r');
                % txt = strcat('Theta:', num2str(theta));
                % t3 = text(pos_txt_x, -2 - i/50000, txt,"HorizontalAlignment","left",VerticalAlignment="baseline", Color='r');
               
                if i < init.burnout_n
                     Mp_now = init.Mp - init.variable_mdot * init.t(i);
                     h_fuel = (Mp_now / tankVol) / areaStage ;
                     h_diff = h_tank - h_fuel ;
 
                end 
                
                fuelline_o(1,:) = fuelline_o(1,:) - h_diff ; 
                new_fuelline = fuelline_o; 
               
 

                body = init.body.stored; 

                if sep
                    if i > init.burnout_n 
                        body = init.body.stage1_boundary;             
                    end
                end
                
                bb = body; 
                fb = new_fuelline; 

                break
                
                
                % rocketTheta = rotate2D(body, theta); 
                % fuellineTheta = rotate2D(new_fuelline, theta);
                p.XData = rocketTheta(1,:);  p.YData = rocketTheta(2,:) ;
                
                
                exportgraphics(gcf,gifnm,'Append',true);

            end
    end
end
