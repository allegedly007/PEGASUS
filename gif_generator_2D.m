function gif_generator_2D(init)
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

                rocketTheta = rotate2D(init.body.rocket2Dfilled, theta); 
                p.XData = rocketTheta(1,:);  p.YData = rocketTheta(2,:) ; 
                exportgraphics(gcf,gifnm,'Append',true);

            end
    end
end
