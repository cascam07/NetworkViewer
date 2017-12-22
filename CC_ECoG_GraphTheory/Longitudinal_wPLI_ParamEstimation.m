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

segment_l = 20;
trial_l = 2;
rel_thr = 0.3; %relative threshold for most alert state

%% wPLI Parameter Screen
patientIDs = fieldnames(batchParams)';
patientSet = 1:2; %Run on first 2 subjects for test


VarTrack = [];
seglengths = [10, 20, 40, 70, 100, 150, 200, 300];
for segment_l = seglengths
    ParamTrack = [];
    ECoG_conn = [];
disp(['Beginning Segment Length ' segment_l ' seconds'])
for iPatient = patientSet
    thisName = patientIDs{iPatient};
    params = batchParams.(patientIDs{iPatient});
    conditions = fieldnames(params)';
    conditions = conditions(1); %only look at OAAS5 for now
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
        
        %Redefine data into segments, each of which will be further subdivided into trials to calculate a wPLI value
        %for each segment
        
        cfg         = [];
        cfg.length  = segment_l; %trial length in ms
        cfg.overlap = params.(thisCond).trialOverlap; %1/4 segment overlap
        cfg.feedback   = 'no';
        data_ECoGds   = ft_redefinetrial(cfg, data_ECoGds);
        
        %For each segment, divide into trials and calculate wPLI for each segment
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
                               
                %Extract data for tracking
                ParamTrack.(thisName).(thisCond).(thisBand).(['segl_' num2str(segment_l)])(:,:,counter) = ECoG_conn.WPLI.(thisName).(thisCond).(['seg_' num2str(iSegment)]).(thisBand).wpli_debiasedspctrm;

            end            
        end   
        
        for iBand = 1:length(bands)
            thisBand = bands{iBand};
            A = ParamTrack.(thisName).(thisCond).(thisBand).(['segl_' num2str(segment_l)]);
            trial_var = mean(mean(A, 'omitnan'));
            VarTrack.(thisName).(thisCond).(thisBand).(['segl_' num2str(segment_l)]).Value = squeeze(trial_var)';
            VarTrack.(thisName).(thisCond).(thisBand).(['segl_' num2str(segment_l)]).Time = cellfun(@(x) x(1),squeeze(data_ECoGds.time));
            disp(['Seg Length: ' num2str(segment_l) '   Patient: ' thisName '   Condition: ' thisCond '   Band: ' thisBand])
        end
        
        
    end %Loop over conditions
    batchParams.(patientIDs{iPatient}) = params;
end

end

%% Plot WPLI means as function of time and segment length
bands = {'alpha','beta','delta','gamma','theta'};
nbands = int32(length(bands));
ncol = 2;
nrow = double(idivide(nbands,ncol,'ceil'));
patients = fieldnames(VarTrack);
for iPatient = 1:length(patients)
    figure(iPatient);
    thisPatient = patients{iPatient};
    for iBand = 1:length(bands)    
        band = bands{iBand};
        subplot(nrow,ncol,iBand);
        plotmeanvar(VarTrack, thisPatient, band,'time');
        title(band);
    end
    suptitle(thisPatient);

end

