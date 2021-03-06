%% Setup Environment

%Change to fit your work environment
computerSpecPrefix = '/Users/cpcasey3/Documents/ECoG/';

if ~exist('defaultPath','var')
    dataPath = [computerSpecPrefix 'Data and analysis' filesep 'ECoG data' filesep];
end
if ~exist('electrodeFilePath','var')
    electrodeFilePath = [computerSpecPrefix 'Data and analysis' filesep ...
        'ECoG data' filesep 'Electrode MNI coordinates and ROI' filesep];
end
if ~exist('outPath','var')
    outPath = [computerSpecPrefix 'Data and analysis' filesep 'Longitudinal' filesep];
end
% Ensure output directory exists
if ~exist(outPath, 'dir')
    mkdir(outPath);
end
disp(['Data will be saved to `' outPath '`']);

[batchParams_gen,setName] = getBatchParams(electrodeFilePath);%,defaultPath);
patientIDs_gen = fieldnames(batchParams_gen);

if exist([outPath 'ECoG_conn_' setName '.mat'],'file')
    load([outPath 'ECoG_conn_' setName '.mat'],'ECoG_conn','batchParams');
    patientIDs_load = fieldnames(ECoG_conn.coh);
    %Figure out if there are new patients in getBatchParams
    genIndx = ~ismember(patientIDs_gen,patientIDs_load);
    patientsToAdd = patientIDs_gen(genIndx);
    if ~isempty(patientsToAdd)
        for iPatient = 1:length(patientsToAdd)
            batchParams.(patientsToAdd{iPatient}) = batchParams_gen.(patientsToAdd{iPatient});
        end
    end
    %Now get indices for added patients, in case just want to analyze those
    patientIDs = fieldnames(batchParams);
    patientIndx = ~ismember(patientIDs,patientIDs_load);
else
    ECoG_conn = struct;
    batchParams = batchParams_gen;
    patientIDs = fieldnames(batchParams);
    patientIndx = 1:length(patientIDs);
end

% Freq analysis params
tinyCriterion = 0.01;
freqsOfInterest = [2,6,10,20,40,80]; %delta, theta, alpha, beta, gamma, high gamma
freqSmoothing = [2,2,2,8,10,20]; %delta, theta, alpha, beta, gamma, high gamma
bands = {'delta','theta','alpha','beta','gamma'};%,'highGamma'};
dsFs = 250;
justRunNew = 1;

segment_l = 10; %Length of segments within a 5 minute recording period
trial_l = 2;    %Length of trials within a segment of length segment_l


%% wPLI Analysis
patientIDs = fieldnames(batchParams)';
patientSet = 1:length(patientIDs);

for iPatient = patientSet
    thisName = patientIDs{iPatient};
    params = batchParams.(patientIDs{iPatient});
    conditions = fieldnames(params)';
    for iCond = 1:length(conditions)
        thisCond = conditions{iCond};
        disp('------------------------');
        disp(['Patient ' thisName ' - Condition: ' thisCond]);
        disp('------------------------');
        [loadedData,params.(thisCond)] = loadECoGData(params.(thisCond),[dataPath filesep thisName filesep thisCond filesep]);
        % loadedData is matrix of nChan x nSamples
        [params.(thisCond).chanSide, params.(thisCond).chanROI, params.(thisCond).chanNum,...
            params.(thisCond).sortOrder] = getChannelInfo(params.(thisCond));
        [data_ECoG] = convertToFTFormat(loadedData,params.(thisCond));
        % patient 399R has artifact for first ~10 secs in ctrl condition
        if strcmp(thisName,'patient399R') && iCond == 1
            data_ECoG.trial{1} = data_ECoG.trial{1}(:,1.5e4:end);
            data_ECoG.sampleinfo = [1,size(data_ECoG.trial{1},2)];
            data_ECoG.time{1} = data_ECoG.time{1}(1:data_ECoG.sampleinfo(2));
        end
        
        if params.(thisCond).dT < 0.0025
            % downsample the data to speed up component analysis
            cfg = [];
            cfg.resamplefs = dsFs;
            cfg.detrend    = 'yes';
            % the following line avoids numeric round off issues in the time axes upon resampling
            data_ECoGds = ft_resampledata(cfg, data_ECoG);
            data_ECoGds.sampleinfo = [1, size(data_ECoGds.trial{1,1},2)];
        else
            data_ECoGds = data_ECoG;
        end
        clear data_ECoG
        
        %Redefine data into segments of length segment_l, each of which will be further subdivided into trials to calculate a wPLI value
        %for each segment
        
        cfg         = [];
        cfg.length  = segment_l; %trial length in ms
        cfg.overlap = params.(thisCond).trialOverlap; %1/4 segment overlap
        cfg.feedback   = 'no';
        data_ECoGds   = ft_redefinetrial(cfg, data_ECoGds);
        
        %For each segment, divide into trials of length trial_l and calculate wPLI for each segment
        counter = 0;
        for iSegment = 1:length(data_ECoGds.trial)
            counter = counter + 1;
            
            seg_dat = [];
            seg_dat.trial = data_ECoGds.trial(iSegment);
            seg_dat.time = data_ECoGds.time(iSegment);
            seg_dat.fsample = data_ECoGds.fsample;
            seg_dat.label = data_ECoGds.label;
            seg_dat.sampleinfo = data_ECoGds.sampleinfo(iSegment,:);
            
            cfg         = [];
            cfg.length  = trial_l; %trial length in ms
            cfg.overlap = params.(thisCond).trialOverlap; %1/4 segment overlap
            seg_dat   = ft_redefinetrial(cfg, seg_dat);
            
            for iBand = 1:length(bands)
                thisBand = bands{iBand};
                % compute Fourier spectra, to be used for cross-spectral density computation.
                cfg            = [];
                cfg.method     = 'mtmfft';
                cfg.output     = 'fourier';
                cfg.taper     = 'dpss';
                cfg.foi        = freqsOfInterest(iBand);
                cfg.tapsmofrq  = freqSmoothing(iBand) ;
                cfg.pad        = 'nextpow2';
                cfg.keeptrials = 'yes';
                cfg.feedback   = 'no';
                FFT_ECoG       = ft_freqanalysis(cfg, seg_dat);

                % compute connectivity as wpli
                cfg         = [];
                cfg.method  ='wpli_debiased';
                cfg.jackknife = 'yes';
                cfg.feedback   = 'no';
                ECoG_conn.WPLI.(thisName).(thisCond).(['seg_' num2str(iSegment)]).(thisBand) = ft_connectivityanalysis(cfg,FFT_ECoG);
                               
            end            
        end   
                
        
    end %Loop over conditions
    batchParams.(patientIDs{iPatient}) = params;
end
save([outPath 'ECoG_wPLI_long_' num2str(segment_l) '.mat'],'-v7.3','ECoG_conn','batchParams');


%% Add Networks to ECoG_conn object based on adjacency matrices
adjmatname = 'wpli_debiasedspctrm';
patientIDs = fieldnames(ECoG_conn.WPLI);
for iPatient = 1:length(patientIDs)
    thisName = patientIDs{iPatient};
    conditions = fieldnames(ECoG_conn.WPLI.(thisName));
    for iCond = 1:length(conditions) 
        thisCond = conditions{iCond};
       
        segs = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond));
        for iSegment = 1:length(segs)
            thisSeg = segs{iSegment};
            if(strcmp(thisSeg, 'NetworkStats')) continue; end %Skip NetworkStats if these already exist
            bands = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg));
            for iBand = 1:length(bands)
               thisBand = bands{iBand}; 
               adjmat = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).(adjmatname);
               
               if(any(any(adjmat < 0))) %If network has positive and negative values, make negative values zero
                   adjmat(adjmat < 0) = 0; 
               end
               
               %Calculate absolute threshold from subject's most alert state
               if iCond == 1
                   num_pos_edges = (length(adjmat).^2)/2; %possible number of edges
                   perc_thr = PercolationThr(adjmat, 'cutoff');
                   thr_ind = round(num_pos_edges * perc_thr);
                   edges = triu(adjmat);
                   edges(isnan(edges)) = [];
                   edges = sort(edges(:), 'descend');
                   abs_thr = edges(thr_ind);
                   ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).AbsThr = abs_thr;
               else
                   most_alert = fieldnames(ECoG_conn.WPLI.(thisName));
                   most_alert = most_alert(1);
                   if isfield(ECoG_conn.WPLI.(thisName).(most_alert{1}),thisSeg)
                       abs_thr = ECoG_conn.WPLI.(thisName).(most_alert{1}).(thisSeg).(thisBand).AbsThr;
                   else
                       maxseg = fieldnames(ECoG_conn.WPLI.(thisName).(most_alert{1}));
                       maxseg = maxseg{end};
                       abs_thr = ECoG_conn.WPLI.(thisName).(most_alert{1}).(maxseg).(thisBand).AbsThr;
                   end
                   ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).AbsThr = abs_thr;
               end

               ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Networks.Ref_Thresholded = threshold_absolute(adjmat, abs_thr);
               ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Networks.Ref_Binarized = weight_conversion(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Networks.Ref_Thresholded, 'binarize');
            end
        end
    end
end        

%% Calculate Network Stats (Avg. Degree, Transitivity, Char. Path, Assortitivity Coef., Clustering Coef., Modularity, Small World Propensity)
adjmatname = 'wpli_debiasedspctrm';
patientIDs = fieldnames(ECoG_conn.WPLI);
for iPatient = 1:length(patientIDs)
    thisName = patientIDs{iPatient};
    conditions = fieldnames(ECoG_conn.WPLI.(thisName));
    for iCond = 1:length(conditions) 
        thisCond = conditions{iCond};
        
        %Initialize Storage Variables
        avgdeg_wei = [];
        trans_wei = [];
        characterpath_wei = [];
        assort_wei = [];
        clustco_wei = [];
        modular_wei = [];
        swp_wei = [];
        
        avgdeg_bin = [];
        trans_bin = [];
        characterpath_bin = [];
        assort_bin = [];
        clustco_bin = [];
        modular_bin = [];
        swp_bin = [];
                      
        segs = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond));
        for iSegment = 1:length(segs)
            thisSeg = segs{iSegment};
            if(strcmp(thisSeg, 'NetworkStats')) continue; end %Skip NetworkStats if these already exist
            bands = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg));
            for iBand = 1:length(bands)
               thisBand = bands{iBand};
               thisThr = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).AbsThr;

               %Compute Stats on Weighted Networks 
               if(~isfield(avgdeg_wei,thisBand))        avgdeg_wei.(thisBand) = [];          end              
               if(~isfield(trans_wei,thisBand))         trans_wei.(thisBand) = [];           end
               if(~isfield(characterpath_wei,thisBand)) characterpath_wei.(thisBand) = [];   end
               if(~isfield(assort_wei,thisBand))        assort_wei.(thisBand) = [];          end
               if(~isfield(clustco_wei,thisBand))       clustco_wei.(thisBand) = [];         end
               if(~isfield(modular_wei,thisBand))       modular_wei.(thisBand) = [];         end
               if(~isfield(swp_wei,thisBand))           swp_wei.(thisBand) = [];             end               
                                           
               refNet = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Networks.Ref_Thresholded;
 
               %Average Degree
               avgdeg_wei.(thisBand) = [avgdeg_wei.(thisBand) ; iSegment, mean(degrees_und(refNet))];
               %Transitivity
               trans_wei.(thisBand) = [trans_wei.(thisBand) ; iSegment, transitivity_wu(refNet)];
               %Characteristic Path Length
               characterpath_wei.(thisBand) = [characterpath_wei.(thisBand) ; iSegment, charpath(distance_wei(refNet), 0, 0)];                   
               %Assortitivity Coefficient                   
               assort_wei.(thisBand) = [assort_wei.(thisBand) ; iSegment, assortativity_wei(refNet,0)];                   
               %Clustering Coefficient
               clustco_wei.(thisBand) = [clustco_wei.(thisBand) ; iSegment, mean(clustering_coef_wu(refNet))];
               %Modularity
               [Ci, Q] = modularity_und(refNet, 1);
               modular_wei.(thisBand) = [modular_wei.(thisBand) ; iSegment, Q];
               %Small-World Propensity
               [SWP,delta_C,delta_L] = small_world_propensity(refNet);
               swp_wei.(thisBand) = [swp_wei.(thisBand) ; iSegment, SWP];
                           
               %Compute Stats on Unweighted Networks 
               if(~isfield(avgdeg_bin,thisBand))        avgdeg_bin.(thisBand) = [];          end              
               if(~isfield(trans_bin,thisBand))         trans_bin.(thisBand) = [];           end
               if(~isfield(characterpath_bin,thisBand)) characterpath_bin.(thisBand) = [];   end
               if(~isfield(assort_bin,thisBand))        assort_bin.(thisBand) = [];          end
               if(~isfield(clustco_bin,thisBand))       clustco_bin.(thisBand) = [];         end
               if(~isfield(modular_bin,thisBand))       modular_bin.(thisBand) = [];         end
               if(~isfield(swp_bin,thisBand))           swp_bin.(thisBand) = [];             end               
                                            
               refNet = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Networks.Ref_Binarized;
 
               %Average Degree
               avgdeg_bin.(thisBand) = [avgdeg_bin.(thisBand) ; iSegment, mean(degrees_und(refNet))];
               %Transitivity
               trans_bin.(thisBand) = [trans_bin.(thisBand) ; iSegment, transitivity_bu(refNet)];
               %Characteristic Path Length
               characterpath_bin.(thisBand) = [characterpath_bin.(thisBand) ; iSegment, charpath(distance_bin(refNet), 0, 0)];                   
               %Assortitivity Coefficient                   
               assort_bin.(thisBand) = [assort_bin.(thisBand) ; iSegment, assortativity_bin(refNet,0)];                   
               %Clustering Coefficient
               clustco_bin.(thisBand) = [clustco_bin.(thisBand) ; iSegment, mean(clustering_coef_bu(refNet))];
               %Modularity
               [Ci, Q] = modularity_und(refNet, 1);
               modular_bin.(thisBand) = [modular_bin.(thisBand) ; iSegment, Q];
               %Small-World Propensity
               [SWP,delta_C,delta_L] = small_world_propensity(refNet);
               swp_bin.(thisBand) = [swp_bin.(thisBand) ; iSegment, SWP];               
               
            end
        end       
       %Store Network Statistics is ECoG_conn object
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Thresholded.AvgDeg = avgdeg_wei;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Thresholded.Transitivity = trans_wei;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Thresholded.CharPathLen = characterpath_wei;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Thresholded.Assortitivity = assort_wei;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Thresholded.ClusteringCoef = clustco_wei;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Thresholded.Modularity = modular_wei;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Thresholded.SWP = swp_wei;

       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Binarized.AvgDeg = avgdeg_bin;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Binarized.Transitivity = trans_bin;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Binarized.CharPathLen = characterpath_bin;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Binarized.Assortitivity = assort_bin;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Binarized.ClusteringCoef = clustco_bin;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Binarized.Modularity = modular_bin;
       ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.Binarized.SWP = swp_bin;
       disp(['Completed: ', thisName, ' ', thisCond])
    end
end

save([outPath 'ECoG_wPLI_long_' num2str(segment_l) '.mat'],'-v7.3','ECoG_conn','batchParams');


%% Plot Graph Metrics

bands = {'alpha','beta','delta','gamma','theta'};
nbands = int32(length(bands));
ncol = 2;
nrow = double(idivide(nbands,ncol,'ceil'));
patients = fieldnames(ECoG_conn.WPLI);

%'AvgDeg', 'Transitivity', 'CharPathLen', 'Assortivity', 'ClusteringCoef', 'Modularity', or 'SWP'
thisMetric = 'Modularity';

for iPatient = 1:length(patients)
    figure(iPatient);
    thisPatient = patients{iPatient};
    for iBand = 1:length(bands)    
        band = bands{iBand};
        subplot(nrow,ncol,iBand);
        plotmetric(ECoG_conn, thisPatient, 'w', thisMetric, band);
        title(band);
    end
    suptitle([thisPatient, ' ', thisMetric]);
end

%% Calculate Confidence Intervals and one-way-ANOVA
GraphStats = [];
GraphStatsMat = {'Patient', 'NetworkType', 'GraphMetric', 'FreqBand', 'p-Value'};
patientIDs = fieldnames(ECoG_conn.WPLI);
for iPatient = 1:length(patientIDs)
    thisName = patientIDs{iPatient};
    
    conditions = fieldnames(ECoG_conn.WPLI.(thisName));
    condAnova = [];
    for iCond = 1:length(conditions) 
        thisCond = conditions{iCond};        
        
        nets = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats);
        for iNet = 1:length(nets)
            thisNet = nets{iNet};
            
            metrics = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.(thisNet));
            for iMetric = 1:length(metrics)
                thisMetric = metrics{iMetric};
                
                bands = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.(thisNet).(thisMetric));
                for iBand = 1:length(bands)
                   thisBand = bands{iBand};             
                   
                   values = ECoG_conn.WPLI.(thisName).(thisCond).NetworkStats.(thisNet).(thisMetric).(thisBand);
                   values = values(:,2);
                   
                   mn = mean(values);                               % Mean
                   SEM = std(values)/sqrt(length(values));          % Standard Error
                   ts = tinv([0.025  0.975],length(values)-1);      % T-Score
                   CI = mn + ts*SEM;                                % Confidence Intervals
        
                   GraphStats.WPLI.(thisName).(thisNet).(thisCond).(thisMetric).(thisBand).Mean = mn;
                   GraphStats.WPLI.(thisName).(thisNet).(thisCond).(thisMetric).(thisBand).CI = CI;
                   
                   %If field doesn't exist, initialize it first then populate it
                   try
                       condAnova.WPLI.(thisNet).(thisMetric).(thisBand).Data = [condAnova.WPLI.(thisNet).(thisMetric).(thisBand).Data, values(1:length(condAnova.WPLI.(thisNet).(thisMetric).(thisBand).Data))];
                   catch             
                       condAnova.WPLI.(thisNet).(thisMetric).(thisBand).Data = [];
                       condAnova.WPLI.(thisNet).(thisMetric).(thisBand).Data = [condAnova.WPLI.(thisNet).(thisMetric).(thisBand).Data, values];
                   end                  
                end                               
            end           
        end        
    end
    
    nets = fieldnames(condAnova.WPLI);
    for iNet = 1:length(nets)
        thisNet = nets{iNet};

        metrics = fieldnames(condAnova.WPLI.(thisNet));
        for iMetric = 1:length(metrics)
            thisMetric = metrics{iMetric};

            bands = fieldnames(condAnova.WPLI.(thisNet).(thisMetric));
            for iBand = 1:length(bands)
               thisBand = bands{iBand};
               anova_data = condAnova.WPLI.(thisNet).(thisMetric).(thisBand).Data;
               [p, tbl, stats] = anova1(anova_data, conditions, 'off');
               GraphStats.WPLI.(thisName).(thisNet).ANOVA.(thisMetric).(thisBand).p = p;
               GraphStats.WPLI.(thisName).(thisNet).ANOVA.(thisMetric).(thisBand).tbl = tbl;
               GraphStats.WPLI.(thisName).(thisNet).ANOVA.(thisMetric).(thisBand).stats = stats;
               disp([thisName ' ' thisNet ' ' thisMetric ' ' thisBand '   ' 'p-Value: ' num2str(p)])
               GraphStatsMat = [GraphStatsMat; thisName, thisNet, thisMetric, thisBand, num2cell(p)];
            end                               
        end           
    end 
end  
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(cell2mat(GraphStatsMat(2:end,5)));
GraphStatsMat = [GraphStatsMat ['p.adj.FDR' num2cell(adj_p)']'];
adj_p = bonf_holm(cell2mat(GraphStatsMat(2:end,5)));
GraphStatsMat = [GraphStatsMat ['p.adj.Bonf' num2cell(adj_p)']'];

%% 3D Brain Plotting Demo
atlas = '/Users/cpcasey3/Documents/MATLAB/ToolBoxes/fieldtrip-20171031/template/anatomy/surface_pial_both.mat';
adj = ECoG_conn.WPLI.patient394R.OAAS5.seg_2.alpha.Networks.Ref_Thresholded;
plotbrain(adj, atlas, batchParams, 'patient394R', 0.25)

