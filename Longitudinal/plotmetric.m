function metplot = plotmetric(ECoG_conn, patient, nettype, gmetric, band)
    
    if(strcmp(nettype, 'w')) networktype = 'Thresholded';
    elseif(strcmp(nettype, 'u')) networktype = 'Binzarized';
    else error('Invalid network type. Choose ''w'' for weighted or ''u'' for unweighted');
    end
    
    plotdata = [];

    conditions = fieldnames(ECoG_conn.WPLI.(patient));
    offset = 0;
    for iCond = 1:length(conditions)
        thisCond = conditions{iCond};
        
        cond_dat = ECoG_conn.WPLI.(patient).(thisCond).NetworkStats.(networktype).(gmetric).(band);
        metplot = plot(cond_dat(:,1)+offset,cond_dat(:,2), 'o');
        hold on
        offset = offset + max(cond_dat(:,1));
       
    end



end