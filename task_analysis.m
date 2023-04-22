function anIn = task_analysis(anIn)
%% Function to analyse ripples identified by ripple_analysis.m

    %% Specify default durations (in seconds)
    if ~isfield(anIn,'rippleDir') || ~isfield(anIn,'taskDataDir') || ~isfield(anIn,'outputDir') || ~isfield(anIn,'allPtsDir')
        error('File directory not specified.'); end
    if ~isfield(anIn,'studyfixDur');        anIn.studyfixDur     = 2; end   % Duration (s) of Study fixation period
    if ~isfield(anIn,'studyDur');           anIn.studyDur        = 6; end   % Duration (s) of Study period
    if ~isfield(anIn,'studyBlankDur');      anIn.studyBlankDur   = 2; end   % Duration (s) of Study blank period
    if ~isfield(anIn,'testfixDur');         anIn.testfixDur      = 0.5; end % Duration (s) of Test fixation period
    if ~isfield(anIn,'studyRetrieveDur');   anIn.testRetrieveDur = 3; end   % Duration (s) of Test retrieve (Cue picture) period
    if ~isfield(anIn,'respDur');            anIn.respDur         = 1; end   % Duration (s) of time before Response periods to save
    if ~isfield(anIn,'consBuff');           anIn.consBuff        = 10; end  % Buffer for at start and end of consolidation period
    if ~isfield(anIn,'corThresh');          anIn.corThresh       = 2; end   % Threshold for number of correct/wrong trials
    
    % Find all mat files in specified output folder.
    matFilePattern = fullfile(anIn.rippleDir, '*.mat');
    matFiles       = dir(matFilePattern);
    clear matFilePattern
    
    % Define cell to save all patient data (allPts)
    noPts = length(matFiles);
    allPts.Patients = cell(noPts,1);
    propStructAll = struct('rippleFreq',zeros(noPts,1),'meanAmp',zeros(noPts,1),'medianAmp',zeros(noPts,1),...
                           'meanDur',zeros(noPts,1));
    durStructAll  = struct('StudyBaseline',propStructAll,'StudyTask',propStructAll,'StudyTask2',propStructAll,'StudyTask0',propStructAll,'Consolidation',propStructAll,...
                           'TestBaseline',propStructAll,'TestTask',propStructAll,'TestTask2',propStructAll,'TestTask0',propStructAll,'TestTask01',propStructAll,...
                           'TestTaskOld',propStructAll,'TestTaskNew',propStructAll,'TestTaskFastRT',propStructAll,'TestTaskSlowRT',propStructAll,...
                           'TestBefRec',propStructAll,'TestBefAss1',propStructAll,'TestBefAss2',propStructAll, ...
                           'StudyFirst',propStructAll,'StudySub',propStructAll);
    clear propStructAll
    allPts.AxChan = durStructAll; clear durStructAll
    
    for file_no = 1:length(matFiles)
        
        % Save file names
        anOut.fname = matFiles(file_no).name;
    
        % Load Output file and Task Data file
        load(char(fullfile(anIn.rippleDir, matFiles(file_no).name)));
        load(char(fullfile(anIn.taskDataDir,strcat(out.file,'-TaskData.mat'))));

        disp(strcat('Currently running: ', out.file))
        
        % Save overall ripple frequency and overall average ripple frequency
        anOut.overallRippleFreq     = out.RippleFreq;
        anOut.overallAvgRippleFreq  = out.avgRippleFreq;
        anOut.overallRippleDur      = out.RippleDur;
        anOut.overallAvgRippleDur   = out.avgRippleDur;
        anOut.overallRippleAmp      = out.RippleAmp;
        anOut.overallAvgRippleAmp   = out.avgRippleAmp;

        % Define output (anOut)
        noChan = length(out.Channels);
        propStruct = struct('rippleFreq',nan(noChan,1),'meanAmp',nan(noChan,1),'medianAmp',nan(noChan,1),...
                            'meanDur',nan(noChan,1));
        durStruct  = struct('StudyBaseline',propStruct,'StudyTask',propStruct,'StudyTask2',propStruct,'StudyTask0',propStruct,'Consolidation',propStruct,...
                            'TestBaseline',propStruct,'TestTask',propStruct,'TestTask2',propStruct,'TestTask0',propStruct,'TestTask01',propStruct,...
                            'TestTaskOld',propStruct,'TestTaskNew',propStruct,'TestTaskFastRT',propStruct,'TestTaskSlowRT',propStruct,...
                            'TestBefRec',propStruct,'TestBefAss1',propStruct,'TestBefAss2',propStruct, ...
                            'StudyFirst',propStruct,'StudySub',propStruct);
        clear propStruct
        anOut.analysisByChan = durStruct; clear durStruct
        propStructAx = struct('rippleFreq',[],'meanAmp',[],'medianAmp',[],'meanDur',[]);
        durStructAx  = struct('StudyBaseline',propStructAx,'StudyTask',propStructAx,'StudyTask2',propStructAx,'StudyTask0',propStructAx,'Consolidation',propStructAx,...
                              'TestBaseline',propStructAx,'TestTask',propStructAx,'TestTask2',propStructAx,'TestTask0',propStructAx,'TestTask01',propStructAx,...
                              'TestTaskOld',propStructAx,'TestTaskNew',propStructAx,'TestTaskFastRT',propStructAx,'TestTaskSlowRT',propStructAx,...
                              'TestBefRec',propStructAx,'TestBefAss1',propStructAx,'TestBefAss2',propStructAx, ...
                              'StudyFirst',propStructAx,'StudySub',propStructAx);
        anOut.analysisAxChan = durStructAx; 
    
        %% Task characteristics
        % Find subAss total (2 = both correct, 0 = both wrong)
        studyTotal = sum(arts.Study.subAss,2);
        testTotal  = arts.Test.assAcc1(:) + arts.Test.assAcc2(:);
        
        % Sanity check: nansum(studyTotal)==nansum(testTotal)
    
        % Specify number of trials
        anOut.noStudy = length(arts.Study.Times);

        % histcounts(studyTotal,0:3) will give the number of trials in each condition (0 / 1 / 2 correct)
        studyHist           = histcounts(studyTotal,0:3);
        anOut.noStudy2      = studyHist(3); 
        anOut.noStudy1      = studyHist(2);
        anOut.noStudy0      = studyHist(1); clear studyHist
        anOut.noStudyFirst  = sum(arts.Study.Pres(:)==1);
        anOut.noStudySub    = sum(arts.Study.Pres(:)~=1);
    
        anOut.noTest        = length(arts.Test.Times);
        testHist            = histcounts(testTotal,0:3);
        anOut.noTestOld     = sum(~isnan(arts.Test.Times(:,4)));
        anOut.noTestNew     = sum(arts.Test.OldNew,'all');
        anOut.noTestOld2    = testHist(3);
        anOut.noTestOld1    = testHist(2);
        anOut.noTestOld0    = testHist(1); 
        anOut.noTestOld01   = anOut.noTestOld0 + anOut.noTestOld1; clear testHist
    
        allRTs = arts.Test.assRT1(:) + arts.Test.assRT2(:); % Total reaction time
        anOut.noTestOldFast = sum(allRTs<nanmedian(allRTs));
        anOut.noTestOldSlow = sum(allRTs>=nanmedian(allRTs));
    
        % Calculate time when responses were made
        oldNewResponseTime       = arts.Test.Times(:,3) + arts.Test.recRT(:)./1000;
        associativeResponseTime1 = arts.Test.Times(:,4) + arts.Test.assRT1(:)./1000;
        associativeResponseTime2 = arts.Test.Times(:,5) + arts.Test.assRT2(:)./1000;
    
        % Create vector of individual durations corresponding to labels
        indDur = [anIn.studyfixDur anIn.studyDur anIn.studyDur anIn.studyDur nan ...
                  anIn.testfixDur anIn.testRetrieveDur anIn.testRetrieveDur anIn.testRetrieveDur anIn.testRetrieveDur...
                  anIn.testRetrieveDur anIn.testRetrieveDur anIn.testRetrieveDur anIn.testRetrieveDur ...
                  anIn.respDur anIn.respDur anIn.respDur ...
                  anIn.studyDur anIn.studyDur];
        % Calculate duration of consolidation period
        indDur(5) = (arts.Test.Times(1,1)-anIn.consBuff) - (arts.Study.Times(end,2)+anIn.consBuff);
    
        % Create vector of number of each trial corresponding to labels
        noTrials = [anOut.noStudy anOut.noStudy anOut.noStudy2 anOut.noStudy0 1 ...
                    anOut.noTest anOut.noTest anOut.noTestOld2 anOut.noTestOld0 anOut.noTestOld01...
                    anOut.noTestOld anOut.noTestNew anOut.noTestOldFast anOut.noTestOldSlow ...
                    anOut.noTest anOut.noTestOld anOut.noTestOld ...
                    anOut.noStudyFirst anOut.noStudySub];
        
        % Create vector of total duration corresponding to labels 
        durs = indDur .* noTrials;
    
        % Check patient has at least corThresh trials each where testTotal==0 and testTotal==2
        if anOut.noTestOld0 < anIn.corThresh || anOut.noTestOld2 < anIn.corThresh
            warning(strcat("Patient ", matFiles(file_no).name," has fewer than ",num2str(anIn.corThresh), " correct/wrong trials."))
        end
    
    
        %% Labeling ripples
        for c = 1 : noChan

            % Check that channel has at least the minimum rate of ripples
            % Any channel that do not fulfill this criterion will have nan recorded in the corresponding out.RippleFreq element
            if ~isnan(out.RippleFreq(c))
                
                % Labeling study trial number
                studyTrialNo = zeros(1,length(out.Channels{c}));
                for row = 1:length(arts.Study.Times)       
                    studyTrialMask = int8([out.Channels{c}.StartTime] > arts.Study.Times(row,2) & ...
                                          [out.Channels{c}.EndTime] < arts.Study.Times(row,2)+anIn.studyDur);
                    studyTrialNo(studyTrialMask==1) = row; clear studyTrialMask
                end, clear row
                studyTrialNo = num2cell(studyTrialNo);
                [out.Channels{c}.StudyTrialNo] = studyTrialNo{:}; clear studyTrialNo
        
                % Labeling test trials based on study trial numbers 
                testTrialNo = zeros(1,length(out.Channels{c}));
                for row = 1:length(arts.Test.Times)       
                    testTrialMask = int8([out.Channels{c}.StartTime] > arts.Test.Times(row,2) & ...
                                         [out.Channels{c}.EndTime] < arts.Test.Times(row,3));
                    testTrialNo(testTrialMask==1) = row; clear testTrialMask
                end, clear row
                testTrialNo = num2cell(testTrialNo);
                [out.Channels{c}.TestTrialNo] = testTrialNo{:}; clear testTrialNo
        
                % Label each ripple by portion of task it is in
                % Each column is a label
                studyBaseline   = num2cell(any([out.Channels{c}.StartTime] > arts.Study.Times(:,1) & ...
                                               [out.Channels{c}.EndTime] < arts.Study.Times(:,2),1));
                [out.Channels{c}.StudyBaseline] = studyBaseline{:}; clear studyBaseline
                studyTask       = num2cell(any([out.Channels{c}.StartTime] > arts.Study.Times(:,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Study.Times(:,2)+anIn.studyDur,1));
                [out.Channels{c}.StudyTask] = studyTask{:}; clear studyTask
                studyTask2      = num2cell(any([out.Channels{c}.StartTime] > arts.Study.Times(studyTotal(:)==2,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Study.Times(studyTotal(:)==2,2)+anIn.studyDur,1));
                [out.Channels{c}.StudyTask2] = studyTask2{:}; clear studyTask2
                studyTask1      = num2cell(any([out.Channels{c}.StartTime] > arts.Study.Times(studyTotal(:)==1,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Study.Times(studyTotal(:)==1,2)+anIn.studyDur,1));
                [out.Channels{c}.StudyTask1] = studyTask1{:}; clear studyTask1
                studyTask0      = num2cell(any([out.Channels{c}.StartTime] > arts.Study.Times(studyTotal(:)==0,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Study.Times(studyTotal(:)==0,2)+anIn.studyDur,1));
                [out.Channels{c}.StudyTask0] = studyTask0{:}; clear studyTask0

                studyFirst      = num2cell(any([out.Channels{c}.StartTime] > arts.Study.Times(arts.Study.Pres(:)==1,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Study.Times(arts.Study.Pres(:)==1,2)+anIn.studyDur,1));
                [out.Channels{c}.StudyFirst] = studyFirst{:}; clear studyFirst
                studySub        = num2cell(any([out.Channels{c}.StartTime] > arts.Study.Times(arts.Study.Pres(:)~=1,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Study.Times(arts.Study.Pres(:)~=1,2)+anIn.studyDur,1));
                [out.Channels{c}.StudySub] = studySub{:}; clear studySub
        
                cons            = num2cell([out.Channels{c}.StartTime] > arts.Study.Times(end,2)+anIn.consBuff & ...
                                           [out.Channels{c}.EndTime] < arts.Test.Times(1,1)-anIn.consBuff);
                [out.Channels{c}.Consolidation] = cons{:}; clear cons
        
                testBaseline    = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(:,1) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(:,2),1));
                [out.Channels{c}.TestBaseline] = testBaseline{:}; clear testBaseline
                testTask        = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(:,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(:,3),1)); 
                [out.Channels{c}.TestTask] = testTask{:}; clear testTask
                testTask2       = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(testTotal(:)==2,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(testTotal(:)==2,3),1)); 
                [out.Channels{c}.TestTask2] = testTask2{:}; clear testTask2
                testTask1       = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(testTotal(:)==1,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(testTotal(:)==1,3),1)); 
                [out.Channels{c}.TestTask1] = testTask1{:}; clear testTask1
                testTask0       = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(testTotal(:)==0,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(testTotal(:)==0,3),1)); 
                [out.Channels{c}.TestTask0] = testTask0{:}; clear testTask0
                testTask01      = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(testTotal(:)~=2,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(testTotal(:)~=2,3),1)); 
                [out.Channels{c}.TestTask01] = testTask01{:}; clear testTask01
                
                testTaskOld     = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(arts.Test.OldNew==0,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(arts.Test.OldNew==0,3),1)); 
                [out.Channels{c}.TestTaskOld] = testTaskOld{:}; clear testTaskOld
                testTaskNew     = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(arts.Test.OldNew==1,2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(arts.Test.OldNew==1,3),1)); 
                [out.Channels{c}.TestTaskNew] = testTaskNew{:}; clear testTaskNew
        
                
                testTaskFastRT  = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(allRTs<nanmedian(allRTs),2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(allRTs<nanmedian(allRTs),3),1)); 
                [out.Channels{c}.TestTaskFastRT] = testTaskFastRT{:}; clear testTaskFastRT
                testTaskSlowRT  = num2cell(any([out.Channels{c}.StartTime] > arts.Test.Times(allRTs>=nanmedian(allRTs),2) & ...
                                               [out.Channels{c}.EndTime] < arts.Test.Times(allRTs>=nanmedian(allRTs),3),1)); 
                [out.Channels{c}.TestTaskSlowRT] = testTaskSlowRT{:}; clear testTaskSlowRT
        
                testBefRec      = num2cell(any([out.Channels{c}.StartTime] > oldNewResponseTime(:)-anIn.respDur & ...
                                               [out.Channels{c}.EndTime] < oldNewResponseTime(:),1)); 
                [out.Channels{c}.TestBefRec] = testBefRec{:}; clear testBefRec
                testBefAss1     = num2cell(any([out.Channels{c}.StartTime] > associativeResponseTime1(:)-anIn.respDur & ...
                                               [out.Channels{c}.EndTime] < associativeResponseTime1(:),1)); 
                [out.Channels{c}.TestBefAss1] = testBefAss1{:}; clear testBefAss1
                testBefAss2     = num2cell(any([out.Channels{c}.StartTime] > associativeResponseTime2(:)-anIn.respDur & ...
                                               [out.Channels{c}.EndTime] < associativeResponseTime2(:),1)); 
                [out.Channels{c}.TestBefAss2] = testBefAss2{:}; clear testBefAss2
        
                %% Extract ripple attributes based on label
                labels = fieldnames(anOut.analysisByChan);
                for labelNo = 1:length(labels)
                    targetField = labels{labelNo};
                    anOut.analysisByChan.(targetField).rippleFreq(c) = sum([out.Channels{c}.(targetField)]) / durs(labelNo);
                    
                    % Extract amplitude in terms of the Z-scored amplitude
                    anOut.analysisByChan.(targetField).meanAmp(c) = nanmean([out.Channels{c}([out.Channels{c}.(targetField)]).PeakZScore]);
                    anOut.analysisByChan.(targetField).medianAmp(c) = nanmedian([out.Channels{c}([out.Channels{c}.(targetField)]).PeakZScore]);
                    anOut.analysisByChan.(targetField).meanDur(c) = nanmean([out.Channels{c}([out.Channels{c}.(targetField)]).Duration]);
                end, clear labelNo targetField
            end
        end, clear c noChan lables allRTs
    
        %% Analyse across channels
        durLabels = fieldnames(anOut.analysisByChan);
        for durLabelNo = 1:length(durLabels)
            durTargetField  = durLabels{durLabelNo};
            propLabels      = fieldnames(anOut.analysisByChan.(durTargetField));
    
            for propLabelNo = 1:length(propLabels)
                propTargetField = propLabels{propLabelNo};
                anOut.analysisAxChan.(durTargetField).(propTargetField) = nanmean(anOut.analysisByChan.(durTargetField).(propTargetField));
                allPts.AxChan.(durTargetField).(propTargetField)(file_no) = nanmean(anOut.analysisByChan.(durTargetField).(propTargetField));
            end, clear propLabelNo durTargetField propTargetField
    
        end, clear durLabelNo durLabels propLabels
    
        %% Save output
        save(fullfile(anIn.outputDir, strcat(arts.file, '-TaskAnalysis')), 'anOut', '-v7.3'); 
    
        allPts.Patients{file_no} = anOut; 
    
    end, clear file_no matFiles
    
    %% Analyse across patients
    allPts.Overall  = durStructAx; clear durStructAx propStructAx
    durLabels       = fieldnames(allPts.Overall);
        for durLabelNo = 1:length(durLabels)
            durTargetField = durLabels{durLabelNo};
            propLabels     = fieldnames(allPts.Overall.(durTargetField));
    
            for propLabelNo = 1:length(propLabels)

                propTargetField = propLabels{propLabelNo};
                allPts.Overall.(durTargetField).(propTargetField) = nanmean(allPts.AxChan.(durTargetField).(propTargetField));
                
            end, clear propLabelNo durTargetField propTargetField
    
        end, clear durLabelNo durLabels propLabels
    
    %% Save output for analysis across patients
    save(fullfile(anIn.allPtsDir, 'Analysis Across Patients'), 'allPts', '-v7.3'); 

end