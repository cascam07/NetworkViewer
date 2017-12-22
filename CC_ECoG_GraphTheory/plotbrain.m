function bp = plotbrain(net, atlaspath, params, patient, alpha)
    %PLOTBRAIN Plot ECoG electrodes on a model brain and show network connectivity
    %   P = PLOTBRAIN(net, atlaspath, params, patient)
    %       net is an adjacency matrix representation of the network to be plotted
    %       atlaspath is the file path for the FieldTrip brain atlas being used e.g. 
    %                       '/MATLAB/ToolBoxes/fieldtrip-20171031/template/anatomy/surface_pial_both.mat'
    %       params is a batchParams object containing the electrode positions for the patient
    %       patient is the ID for the patient who will be referenced in the params object
    %       alpha is the transparency of the brain with 0 completely transparent and 1 completely opaque
    
    %Grab first alertness condition for patient to get electrode data from later
    cond = fieldnames(params.(patient));
    cond = cond{1};

    %Create Brain Mesh
    mesh = ft_read_atlas(atlaspath);
    ft_plot_mesh(mesh, 'facecolor', [0.781, 0.762, 0.664], 'facealpha', alpha, 'EdgeColor', 'none');
    view([-90 25])
    lighting gouraud
    material shiny
    camlight

    %Overlay Electrodes
    elecpos = [params.(patient).(cond).ECoGchannels.xCoord];
    elecpos = [elecpos; params.(patient).(cond).ECoGchannels.yCoord];
    elecpos = [elecpos; params.(patient).(cond).ECoGchannels.zCoord]';
    elecnames = {params.(patient).(cond).ECoGchannels.AnatReg};
    elecnames = [elecnames; cellfun(@num2str, {params.(patient).(cond).ECoGchannels.chanNum}, 'UniformOutput', false)]';
    elecnames = strcat(elecnames(:,1),'_',elecnames(:,2));
    elec = [];
    elec.label = elecnames; 
    elec.chanpos = elecpos;
    elec.elecpos = elecpos;
    hs = ft_plot_sens(elec);

    %Draw Edges Between Electrodes
    edges = adj2edge(net);
    for edj = 1:length(edges)
        a = edges(edj,1);
        b = edges(edj,2);
        wei = edges(edj,3);

        coords = elecpos(a,:);
        coords = [coords ; elecpos(b,:)];

        line(coords(:,1),coords(:,2),coords(:,3),'LineWidth',wei.*2)
    end

end