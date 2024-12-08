function init = reporting(init)
    dataframe = [init.t', init.x', init.y', init.vx', init.vy', init.axs', init.ays', init.theta', init.Machs' , init.Drags' , init.Thrusts', init.mass', init.rhos', init.ReNums'] ;
    
    plotTrajectory(init);

    init = telemetry_results(init) ; 

    init.dataframe = array2table(dataframe) ; 
    init.dataframe.Properties.VariableNames = {'t', 'x', 'y', 'vx' , 'vy', 'axs', 'ays', 'theta', 'Mach', 'Drag', 'Thrust', 'mass', 'rho', 'Re'};
    
    filename_precursor1 = "AE441_Pegasus_LaunchTOBurnout_" + init.Stage_Name + "_Trajectory_SI_" ; 
    exportData(init.dataframe(1:init.burnout_n, :), filename_precursor1, init.export_data)

    filename_precursor2 = "AE441_Pegasus_BurnoutTOPeak_" + init.Stage_Name + "_Trajectory_SI_" ; 
    exportData(init.dataframe(init.burnout_n:init.n2, :), filename_precursor2, init.export_data)

    filename_precursor3 = "AE441_Pegasus_PeakTOImpact_" + init.Stage_Name +"_Trajectory_SI_" ; 
    exportData(init.dataframe(init.n2:end, :), filename_precursor3, init.export_data)

    filename_precursor4 = "Ae441_Pegasus_Parameters" ; 
    exportData(init.T_PLV, filename_precursor4, init.export_data)


end

