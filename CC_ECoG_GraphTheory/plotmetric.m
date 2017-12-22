function metplot = plotmetric(ECoG_conn, patient, nettype, gmetric, band)
    %PLOTMETRIC Plot the progression of a network metric over time.
    %   P = PLOTMETRIC(ECoG_conn, patient, nettype, gmetric, band)
    %       nettype may be 'w' for weighted or 'u' for unweighted
    %       gmetric may be 'AvgDeg', 'Transitivity', 'CharPathLen', 'Assortivity', 'ClusteringCoef', 'Modularity', or 'SWP'
    
    if(strcmp(nettype, 'w')) networktype = 'Thresholded';
    elseif(strcmp(nettype, 'u')) networktype = 'Binzarized';
    else error('Invalid network type. Choose ''w'' for weighted or ''u'' for unweighted');
    end

    conditions = fieldnames(ECoG_conn.WPLI.(patient));
    offset = 0;
    for iCond = 1:length(conditions)
        thisCond = conditions{iCond};
        
        cond_dat = ECoG_conn.WPLI.(patient).(thisCond).NetworkStats.(networktype).(gmetric).(band);
        metplot = plot(cond_dat(:,1)+offset,cond_dat(:,2), 'o');
        offset = offset + max(cond_dat(:,1));
        xlabel('Data Segment');
        ylabel(gmetric);
        hold on;
    end
    legend(conditions);
    hold off;
end