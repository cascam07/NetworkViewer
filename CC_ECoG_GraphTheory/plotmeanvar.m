function varplot = plotmeanvar(VarTrack, patient, band, ptype)
    %PLOTMEANVAR Plot the variablility in the mean of a network wPLI as a function of time or segment length.
    %   P = PLOTMEANVAR(VarTrack, patient, band, ptype)
    %       ptype may be 'time' or 'seglen'
    
    plotdata = [];
    
    cond = fieldnames(VarTrack.(patient));
    cond = cond{1};
    segments = fieldnames(VarTrack.(patient).(cond).(band));


    if(strcmp(ptype, 'time'))
        for iSeg = 1:length(segments)
            thisSeg = segments{iSeg};

            seg_dat = VarTrack.(patient).(cond).(band).(thisSeg);
            varplot = plot(seg_dat.Time, seg_dat.Value);
            xlabel('Time');
            ylabel('wPLI');
            hold on;
        end
        legend(strrep(segments,'_',' '));
        hold off;
    
    elseif(strcmp(ptype, 'seglen'))
        for iSeg = 1:length(segments)
            thisSeg = segments{iSeg};
            seg_len = str2double(strrep(thisSeg,'segl_',''));
            seg_dat = VarTrack.(patient).(cond).(band).(thisSeg);
            varplot = scatter(repmat(seg_len,length(seg_dat.Value),1), seg_dat.Value);
            xlabel('Segment Length');
            ylabel('wPLI');
            hold on;
        end
        legend(strrep(segments,'_',' '));
        hold off;
    
    
    else error('Invalid network type. Choose ''w'' for weighted or ''u'' for unweighted');
    end

end