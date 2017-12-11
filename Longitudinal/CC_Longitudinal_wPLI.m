%% Setup Environment
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

segment_l = 20;
trial_l = 4;

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
        %clear data_ECoG
        
        %Redefine data into 20s segments, each of which will be further subdivided into trials to calculate a wPLI value
        %for each segment
        
        cfg         = [];
        cfg.length  = segment_l; %trial length in sec
        cfg.overlap = params.(thisCond).trialOverlap; %1/4 segment overlap
        data_ECoGds   = ft_redefinetrial(cfg, data_ECoGds);
        
        %For each segment, divide into 4s trials and calculate wPLI for each segment
        for iSegment = 1:length(data_ECoGds.trial)
            seg_dat = [];
            seg_dat.trial = data_ECoGds.trial(iSegment);
            seg_dat.time = data_ECoGds.time(iSegment);
            seg_dat.fsample = data_ECoGds.fsample;
            seg_dat.label = data_ECoGds.label;
            seg_dat.sampleinfo = data_ECoGds.sampleinfo(iSegment,:);
            
            cfg         = [];
            cfg.length  = trial_l; %trial length in sec
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
                FFT_ECoG       = ft_freqanalysis(cfg, seg_dat);

                % compute connectivity as wpli
                cfg         = [];
                cfg.method  ='wpli_debiased';
                cfg.jackknife = 'yes';
                ECoG_conn.WPLI.(thisName).(thisCond).(['seg_' num2str(iSegment)]).(thisBand) = ft_connectivityanalysis(cfg,FFT_ECoG);
            end            
        end       

    end %Loop over conditions
    batchParams.(patientIDs{iPatient}) = params;
end
save([outPath 'ECoG_wPLI_long' setName '.mat'],'ECoG_conn','batchParams');


%% Calculate Percolation Threshold
adjmatname = 'wpli_debiasedspctrm';
for iPatient = 1:length(patientIDs)
    thisName = patientIDs{iPatient};
    %foreach condition, e.g. OAAS5
    conditions = fieldnames(ECoG_conn.WPLI.(thisName));
    for iCond = 1:length(conditions) 
        thisCond = conditions{iCond};
        
        segs = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond));
        for iSegment = 1:length(segs)
            thisSeg = segs{iSegment};
            
            %foreach frequency band, e.g. delta, beta, gamma
            bands = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg));
            for iBand = 1:length(bands)
               thisBand = bands{iBand};

               %Compute the Percolation Threshold for each adjacency matrix
               adjmat = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).(adjmatname);
               if(any(any(adjmat < 0))) %If network has positive and negative values, make negative values positive
                   adjmat(adjmat < 0) = adjmat(adjmat < 0).*-1; 
               end
               ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).PercThr = PercolationThr(adjmat, 'cutoff')*100;
               disp([thisCond, ' ', thisBand,' ', num2str(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).PercThr)])
            end           
        end
    end
    disp(['Completed: ', thisName, ' WPLI'])
end        
save([outPath 'ECoG_wPLI_long' setName '.mat'],'ECoG_conn','batchParams');


%% Build Networks
ECoG_conn = BuildNetworks(ECoG_conn);
save([outPath 'ECoG_conn_' setName '.mat'],'ECoG_conn','batchParams');

%% Add Percolation Networks
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
            bands = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg));
            for iBand = 1:length(bands)
               thisBand = bands{iBand};
               
               %Compute the Percolation Network
               %disp(thisBand)
               %disp(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand))
               adjmat = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).(adjmatname);
               if(any(any(adjmat < 0))) %If network has positive and negative values, make negative values positive
                   adjmat(adjmat < 0) = adjmat(adjmat < 0).*-1; 
               end
               percthr = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).PercThr ./ 100;
               ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Thresholded = threshold_proportional(adjmat, percthr);
               ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Binarized = weight_conversion(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Thresholded, 'binarize');
            end
        end
    end
end        

%% Calculate Network Stats (Avg. Degree, Transitivity, Char. Path, Assortitivity Coef.)
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
            bands = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg));
            for iBand = 1:length(bands)
               thisBand = bands{iBand};
               thisThr = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).PercThr;

               %Compute Stats on Weighted Networks 
               if(~isfield(avgdeg_wei,thisBand))        avgdeg_wei.(thisBand) = [];          end              
               if(~isfield(trans_wei,thisBand))         trans_wei.(thisBand) = [];           end
               if(~isfield(characterpath_wei,thisBand)) characterpath_wei.(thisBand) = [];   end
               if(~isfield(assort_wei,thisBand))        assort_wei.(thisBand) = [];          end
               if(~isfield(clustco_wei,thisBand))       clustco_wei.(thisBand) = [];         end
               if(~isfield(modular_wei,thisBand))       modular_wei.(thisBand) = [];         end
               if(~isfield(swp_wei,thisBand))           swp_wei.(thisBand) = [];             end               
                                           
               refNet = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Thresholded;
 
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
                                            
               refNet = ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).Binarized;
 
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




%% Make matrix of percolation values for histograms
percmat = [];
for iPatient = 1:length(patientIDs)
    thisName = patientIDs{iPatient};
    conditions = fieldnames(ECoG_conn.WPLI.(thisName));
    for iCond = 1:length(conditions) 
        thisCond = conditions{iCond};        
        segs = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond));
        for iSegment = 1:length(segs)
            thisSeg = segs{iSegment};
            bands = fieldnames(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg));
            for iBand = 1:length(bands)
               thisBand = bands{iBand};
               %percmat = [percmat; thisName, thisCond, thisSeg, thisBand, num2str(ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).PercThr)];
               percmat = [percmat, ECoG_conn.WPLI.(thisName).(thisCond).(thisSeg).(thisBand).PercThr];
            end           
        end
    end
end 
