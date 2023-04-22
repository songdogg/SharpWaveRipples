function anIn = spike_analysis(anIn)
%% Function to analyse single neuron activity (spikes)    


    %% Specify parameters
    if ~isfield(anIn,'rippleDir') || ~isfield(anIn,'rippleInputDir') ||~isfield(anIn,'taskDataDir') || ~isfield(anIn,'spikeOutDir') || ~isfield(anIn,'spikeDataDir') || ~isfield(anIn,'outputDir')
        error('File directory not specified.'); end
    if ~isfield(anIn,'gaussWin');   anIn.gaussWin = 0.01; end               % Window for Gaussian smoothing (10ms)
    
    % Find all mat files in specified output folder.
    matFilePattern = fullfile(anIn.rippleDir, '*.mat');
    matFiles = dir(matFilePattern);
    clear matFilePattern
            
    for file_no = 1:length(matFiles)
    
        % Load spike and ripple data
        load(char(fullfile(anIn.rippleDir, matFiles(file_no).name)));
        load(char(fullfile(anIn.rippleInputDir,strcat(out.file,'-Input.mat'))));
        load(char(fullfile(anIn.spikeDataDir,strcat(out.file,'-spikeData.mat'))));
        load(char(fullfile(anIn.taskDataDir,strcat(out.file,'-TaskData.mat'))));
        load(char(fullfile(anIn.outputDir,strcat(out.file,'-TaskAnalysis.mat'))));
        
        disp(strcat('Currently running: ', out.file))
    
        % Predefine cell for subsequent labelling of spikes
        spikeDat.timesWLabel = cell(length(spikeDat.times),1);
    
        % Calculate window to save data
        sampWin = round(-in.saveBuffer*in.downFs):round(in.saveBuffer*in.downFs);
    
        % Make vector of zeros at downsample rate for duration of recording
        spikes = zeros(round(in.downFs*out.recordingDur),length(spikeDat.times));
    
        % Calculate index of spike times
        spikeIdx = cellfun(@(x) round(x*in.downFs), spikeDat.times, 'UniformOutput', false);
    
                % Adapted code from task_analysis (used to label spikes)
                studyTotal               = sum(arts.Study.subAss,2);
                testTotal                = arts.Test.assAcc1(:) + arts.Test.assAcc2(:);
                studyHist                = histcounts(studyTotal,0:3);
                anOut.noStudy2           = studyHist(3); 
                anOut.noStudy1           = studyHist(2);
                anOut.noStudy0           = studyHist(1); clear studyHist
                testHist                 = histcounts(testTotal,0:3);
                anOut.noTestOld          = sum(~isnan(arts.Test.Times(:,4)));
                anOut.noTestNew          = sum(arts.Test.OldNew,'all');
                anOut.noTestOld2         = testHist(3);
                anOut.noTestOld1         = testHist(2);
                anOut.noTestOld0         = testHist(1); 
                anOut.noTestOld01        = anOut.noTestOld0 + anOut.noTestOld1; clear testHist
            
                allRTs                   = arts.Test.assRT1(:) + arts.Test.assRT2(:);
                anOut.noTestOldFast      = sum(allRTs<nanmedian(allRTs));
                anOut.noTestOldSlow      = sum(allRTs>=nanmedian(allRTs));
            
                oldNewResponseTime       = arts.Test.Times(:,3) + arts.Test.recRT(:)./1000;
                associativeResponseTime1 = arts.Test.Times(:,4) + arts.Test.assRT1(:)./1000;
                associativeResponseTime2 = arts.Test.Times(:,5) + arts.Test.assRT2(:)./1000;
        
        for cellNo = 1:length(spikeDat.times)
            % Indicate when the spikes are in spikes vector
            spikes(spikeIdx{cellNo},cellNo) = 1;
    
                % Labelling spikes adapted from ripple labelling (cont'd)
                % Indicate which part of the recording the spikes are in
                spikeDat.timesWLabel{cellNo} = struct('SpikeTime', num2cell(spikeDat.times{cellNo}));
    
                % Adapted code from task_analysis (used to label spikes)
                % Each column is a label
                studyBaseline = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Study.Times(:,1) & ...
                                             [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Study.Times(:,2),1));
                [spikeDat.timesWLabel{cellNo}.StudyBaseline] = studyBaseline{:}; clear studyBaseline
                studyTask = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Study.Times(:,2) & ...
                                         [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Study.Times(:,2)+anIn.studyDur,1));
                [spikeDat.timesWLabel{cellNo}.StudyTask] = studyTask{:}; clear studyTask
                studyTask2 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Study.Times(studyTotal(:)==2,2) & ...
                                          [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Study.Times(studyTotal(:)==2,2)+anIn.studyDur,1));
                [spikeDat.timesWLabel{cellNo}.StudyTask2] = studyTask2{:}; clear studyTask2
                studyTask1 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Study.Times(studyTotal(:)==1,2) & ...
                                          [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Study.Times(studyTotal(:)==1,2)+anIn.studyDur,1));
                [spikeDat.timesWLabel{cellNo}.StudyTask1] = studyTask1{:}; clear studyTask1
                studyTask0 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Study.Times(studyTotal(:)==0,2) & ...
                                          [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Study.Times(studyTotal(:)==0,2)+anIn.studyDur,1));
                [spikeDat.timesWLabel{cellNo}.StudyTask0] = studyTask0{:}; clear studyTask0
        
                cons = num2cell([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Study.Times(end,2)+anIn.consBuff & ...
                                [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(1,1)-anIn.consBuff);
                [spikeDat.timesWLabel{cellNo}.Consolidation] = cons{:}; clear cons
        
                testBaseline = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(:,1) & ...
                                            [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(:,2),1));
                [spikeDat.timesWLabel{cellNo}.TestBaseline] = testBaseline{:}; clear testBaseline
                testTask = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(:,2) & ...
                                        [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(:,3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTask] = testTask{:}; clear testTask
                testTask2 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(testTotal(:)==2,2) & ...
                                         [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(testTotal(:)==2,3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTask2] = testTask2{:}; clear testTask2
                testTask1 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(testTotal(:)==1,2) & ...
                                         [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(testTotal(:)==1,3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTask1] = testTask1{:}; clear testTask1
                testTask0 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(testTotal(:)==0,2) & ...
                                         [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(testTotal(:)==0,3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTask0] = testTask0{:}; clear testTask0
                testTask01 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(testTotal(:)~=2,2) & ...
                                          [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(testTotal(:)~=2,3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTask01] = testTask01{:}; clear testTask01
                
                testTaskOld = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(arts.Test.OldNew==0,2) & ...
                                           [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(arts.Test.OldNew==0,3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTaskOld] = testTaskOld{:}; clear testTaskOld
                testTaskNew = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(arts.Test.OldNew==1,2) & ...
                                           [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(arts.Test.OldNew==1,3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTaskNew] = testTaskNew{:}; clear testTaskNew
        
                
                testTaskFastRT = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(allRTs<nanmedian(allRTs),2) & ...
                                              [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(allRTs<nanmedian(allRTs),3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTaskFastRT] = testTaskFastRT{:}; clear testTaskFastRT
                testTaskSlowRT = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > arts.Test.Times(allRTs>=nanmedian(allRTs),2) & ...
                                              [spikeDat.timesWLabel{cellNo}.SpikeTime] < arts.Test.Times(allRTs>=nanmedian(allRTs),3),1)); 
                [spikeDat.timesWLabel{cellNo}.TestTaskSlowRT] = testTaskSlowRT{:}; clear testTaskSlowRT
        
                testBefRec = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > oldNewResponseTime(:)-anIn.respDur & ...
                                          [spikeDat.timesWLabel{cellNo}.SpikeTime] < oldNewResponseTime(:),1)); 
                [spikeDat.timesWLabel{cellNo}.TestBefRec] = testBefRec{:}; clear testBefRec
                testBefAss1 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > associativeResponseTime1(:)-anIn.respDur & ...
                                           [spikeDat.timesWLabel{cellNo}.SpikeTime] < associativeResponseTime1(:),1)); 
                [spikeDat.timesWLabel{cellNo}.TestBefAss1] = testBefAss1{:}; clear testBefAss1
                testBefAss2 = num2cell(any([spikeDat.timesWLabel{cellNo}.SpikeTime] > associativeResponseTime2(:)-anIn.respDur & ...
                                           [spikeDat.timesWLabel{cellNo}.SpikeTime] < associativeResponseTime2(:),1)); 
                [spikeDat.timesWLabel{cellNo}.TestBefAss2] = testBefAss2{:}; clear testBefAss2
    
        end, clear cellNo spikeIdx allRTs associativeResponseTime2 associativeResponseTime1 oldNewResponseTime studyTotal testTotal
    
        % Smooth spikes vector with a 10ms Gaussian
        smoothSpikes = smoothdata(spikes,1,'gaussian',anIn.gaussWin*in.downFs);
    
        % Calculate z-score
        spikeZsc = zscore(smoothSpikes,0,1); clear smoothSpikes
    
    
        %% Cycle through channels
        for c = 1:length(out.Channels)
            
            % Specify temporary variables to store wanted data
            smSpikeZsc  = zeros(length(out.Channels{c}),length(sampWin),length(spikeDat.times));
            nonSmSpike  = zeros(length(out.Channels{c}),length(sampWin),length(spikeDat.times));
            phaseOutput = cell(length(out.Channels{c}),length(spikeDat.times));
    
            for cellNo = 1:length(spikeDat.times)
                        
                % Extract a 1s window of smoothed Zsc and nonsmoothed spikes around peak of each ripple
                sampBase        = cellfun(@(x) round(x*in.downFs)+sampWin,{out.Channels{c}.PeakTime},'uni',false);
                smSpikeZscCell  = cellfun(@(x) spikeZsc(x,cellNo), sampBase,'UniformOutput',false); 
                nonSmSpikeCell  = cellfun(@(x) logical(spikes(x,cellNo)),sampBase,'UniformOutput',false);clear sampBase
                phaseOutputCell = cellfun(@(x,y) y(x), nonSmSpikeCell, {out.Channels{c}.Phase},'UniformOutput',false);
    
                smSpikeZsc(:,:,cellNo) = cell2mat(smSpikeZscCell)'; clear smSpikeZscCell
                nonSmSpike(:,:,cellNo) = cell2mat(nonSmSpikeCell)'; clear nonSmSpikeCell
                phaseOutput(:,cellNo)  = phaseOutputCell'; clear phaseOutputCell
    
            end, clear cellNo
    
            % Add outputs to out
            smSpikeZsc = num2cell(smSpikeZsc,[2 3]);
            nonSmSpike = num2cell(nonSmSpike,[2 3]);
            smSpikeZsc = cellfun(@squeeze,smSpikeZsc,'UniformOutput',false);
            nonSmSpike = cellfun(@squeeze,nonSmSpike,'UniformOutput',false);
    
            [out.Channels{c}.SmoothSpikeZsc] = smSpikeZsc{:}; clear smSpikeZsc
            [out.Channels{c}.Spikes]         = nonSmSpike{:}; clear nonSmSpike
    
            for row = 1:length(out.Channels{c})
                out.Channels{c}(row).PhaseOutput = phaseOutput(row,:);
            end, clear row phaseOutput
    
        end, clear c spikeZsc spikes sampWin
    
        %% Save output
        save(fullfile(anIn.spikeOutDir, strcat(spikeDat.fname, '-OutputWSpike')), 'out', '-v7.3');
        save(fullfile(anIn.spikeOutDir, strcat(spikeDat.fname, '-LabelledSpikes')), 'spikeDat', '-v7.3');
    
    end, clear file_no matFiles gaussWin

end