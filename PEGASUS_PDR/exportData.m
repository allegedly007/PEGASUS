function exportData(data, filename_precursor, state)
    if state
        date_time = string(datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm-ss')) ; 
        
        filename = filename_precursor + date_time + ".xlsx"  ;  
        fprintf('%s\n', filename)
        
        writetable(data, filename,'Sheet',1, 'WriteRowNames',true)
    end
end

