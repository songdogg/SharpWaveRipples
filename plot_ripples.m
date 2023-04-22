%% Script to plot various figures 

% Toggles for which figures to plot
plotLFP       = false;
plotSpikes    = false;
plotHist      = false;
plotXCorr     = false;
plotSpikeCorr = false;
plotPhases    = false;

load 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Input\20190122-065346-004-Input.mat'
load 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Output\Analysis\TaskSpikeInput.mat'

% Input files to plot
inputDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Input';
outputDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Output';
spikeOutDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Spike Data\Spike Output';

% Find all mat files in specified output folder.
matFilePattern = fullfile(outputDir, '*.mat');
matFiles = dir(matFilePattern);
clear matFilePattern

% Define time axis
plottimeaxis = linspace(in.plotWin(1), in.plotWin(2), (in.plotWin(2)-in.plotWin(1))*in.downFs+1);
plotIdx = (in.saveBuffer*in.downFs+1) + in.plotWin.*in.downFs;

%% Plot Grand average ripple
if plotLFP
    allLFP  = []; 
    allFilt = [];
    allAmp  = [];
    allPow  = [];
    baseWin = [-1 -0.5];                                                    % Window for calculating baseline power (s)
    baseIdx = baseWin .* in.downFs + round(mean(plotIdx));                  % Convert raw time to time index

    for file_no = 1:length(matFiles)
        % Load Output file
        load(char(fullfile(outputDir, matFiles(file_no).name)));
    
        % Extracting relevant information from Channels
        for channel = 1:length(out.Channels)
            if length(out.Channels{channel}) > (in.minRipple * out.recordingDur / 60)
                rawSignalmat    = [out.Channels{channel}(:).RawSignal]';
                rawSignalmat    = rawSignalmat(:,plotIdx(1):plotIdx(2));
                allLFP          = vertcat(allLFP,rawSignalmat); clear rawSignalmat
    
                filtSignalmat   = [out.Channels{channel}(:).FiltSignal]';
                filtSignalmat   = filtSignalmat(:,plotIdx(1):plotIdx(2)); 
                allFilt         = vertcat(allFilt,filtSignalmat); clear filtSignalmat
    
                ampmat          = [out.Channels{channel}(:).ZScore]';
                ampmat          = ampmat(:,plotIdx(1):plotIdx(2)); 
                allAmp          = vertcat(allAmp,ampmat); clear ampmat
    
                powmat          = cat(3,out.Channels{channel}(:).Power);
                % Extract ripple power
                rippowmat       = log(powmat(:,plotIdx(1):plotIdx(2),:));
                % Extract baseline power
                basepowmat      = mean(log(powmat(:,baseIdx(1):baseIdx(2),:)),2);
                % Find difference in power
                diffpowmat      = rippowmat - basepowmat;
                allPow          = cat(3,allPow,diffpowmat); clear powmat

            end
        end, clear channel
    end, clear file_no

    % Plots
    nRipples = num2str(length(allLFP));
    figure
    subplot(2,2,1)
    shadedErrorBar(plottimeaxis, mean(allLFP,1), std(allLFP)./sqrt(size(allLFP,1))), hold on
    title({'Raw LFP'})
    axis square, ylabel('Voltage (µV)'), xlabel('Time (s)')
    plot(zeros(length(-10:15),1), -10:15, 'r--'), hold off
    subplot(2,2,2)
    shadedErrorBar(plottimeaxis, mean(allFilt,1), std(allFilt)./sqrt(size(allFilt,1))), hold on
    title({'Filtered LFP'})
    axis square, ylabel('Voltage (µV)'), xlabel('Time (s)')
    plot(zeros(length(-1:1),1), -1:1, 'r--'), hold off
    subplot(2,2,3)
    shadedErrorBar(plottimeaxis, mean(allAmp,1), std(allAmp)./sqrt(size(allAmp,1))), hold on
    title({'Amplitude'})
    axis square, ylabel('Z-score'), xlabel('Time (s)')
    plot(zeros(length(0:4),1), 0:4, 'r--'), hold off
    subplot(2,2,4)
    contourf(plottimeaxis, in.morletFreq, mean(allPow,3), 50,'linestyle','none'), hold on
    axis square, xlabel('Time (s)'), ylabel('Frequency (Hz)'), title({'Relative Power'})
    c = colorbar; c.Label.String = 'Power (dB)'; c.Label.FontSize = 12; clear c
    plot(zeros(length(0:200),1), 0:200, 'r--'), hold off
end

%% Plot instantaneous firing rates of spikes
if plotSpikes
    allSameChanSpikes = []; 
    allDiffChanSpikes = [];
    for file_no = 1:length(matFiles)
        % Load Output files
        load(char(fullfile(outputDir, matFiles(file_no).name)));
        load(char(fullfile(spikeOutDir,strcat(out.file,'-OutputWSpike.mat'))));
        load(char(fullfile(spikeOutDir,strcat(out.file,'-LabelledSpikes.mat'))));
        load(char(fullfile(inputDir,strcat(out.file,'-Input.mat'))));
    
        % Extracting relevant information from Channels
        for channel = 1:length(out.Channels)
            if ~isnan(out.RippleFreq(channel)) % Check channel is included in analysis
                sameMask            = in.chans2keep(channel)==spikeDat.chanClustID(:,1);

                chanSpikes          = cat(3,out.Channels{1}.SmoothSpikeZsc);
                sameChanSpikes      = chanSpikes(:,sameMask,:);
                sameChanSpikes      = reshape(sameChanSpikes,size(sameChanSpikes,1),[]); % Reduce back into 2-dim
                diffChanSpikes      = chanSpikes(:,~sameMask,:); clear chanSpikes
                diffChanSpikes      = reshape(diffChanSpikes,size(diffChanSpikes,1),[]);

                allSameChanSpikes   = [allSameChanSpikes sameChanSpikes];
                allDiffChanSpikes   = [allDiffChanSpikes diffChanSpikes];
            end
        end, clear channel
    end, clear file_no

    allSameChanSpikes = allSameChanSpikes(plotIdx(1):plotIdx(2),:);
    allDiffChanSpikes = allDiffChanSpikes(plotIdx(1):plotIdx(2),:);

    % Plots
    nSameChan = num2str(length(allSameChanSpikes));
    nDiffChan = num2str(length(allDiffChanSpikes));

    figure
    subplot(1,2,2)
    shadedErrorBar(plottimeaxis, mean(allSameChanSpikes,2), std(allSameChanSpikes,[],2)./sqrt(size(allSameChanSpikes,2)),'lineProps','m'), hold on
    title({'Single cell instaneous firing rate', '(same channel)',strcat('(n=',nSameChan,' combinations)')})
    axis square, ylabel('Z-score'), xlabel('Time (s)'), ylim([-0.02 0.1])
    plot(zeros(2,1), [-0.02 0.1], 'r--'), hold off
    subplot(1,2,1)
    shadedErrorBar(plottimeaxis, mean(allDiffChanSpikes,2), std(allDiffChanSpikes,[],2)./sqrt(size(allDiffChanSpikes,2)),'lineProps','b'), hold on
    title({'Single cell instaneous firing rate', '(different channel)',strcat('(n=',nDiffChan,' combinations)')})
    axis square, ylabel('Z-score'), xlabel('Time (s)'), ylim([-0.02 0.1])
    plot(zeros(2,1), [-0.02 0.1], 'r--'), hold off

end

%% Plot histograms of ripple characteristics
if plotHist
    allRF  = [];
    allMA  = [];
    allPA  = [];
    allDur = [];

    for file_no = 1:length(matFiles)
        % Load Output file
        load(char(fullfile(outputDir, matFiles(file_no).name)));

        % Extract ripple frequency in each channel
        allRF = [allRF out.RippleFreq(~isnan(out.RippleFreq))];
        
        % Extract relevant information from channels
        for channel = 1:length(out.Channels)
            if ~isnan(out.RippleFreq(channel))
                % Ripple mean amplitude
                allMA  = [allMA [out.Channels{channel}.MeanIntensity]];
                % Ripple peak amplitude
                allPA  = [allPA [out.Channels{channel}.MaxIntensity]];
                % Ripple duration
                allDur = [allDur [out.Channels{channel}.Duration]];
            end
        end, clear channel

    end, clear file_no
    
    n_chans = num2str(length(allRF));
    n_ripples = num2str(length(allMA));

    figure
    subplot(2,2,1)
    edges = linspace(0,0.7,29);
    histogram(allRF,edges,'Normalization','probability'), hold on
    plot(repmat(median(allRF),1,length(0:10)), (0:10)/100, 'r--','LineWidth',2), 
    plot(repmat(prctile(allRF,25),length(0:10),1), (0:10)/100, 'r--','LineWidth',1), 
    plot(repmat(prctile(allRF,75),length(0:10),1), (0:10)/100, 'r--','LineWidth',1), hold off
    axis square, yticklabels(yticks*100), title({"Ripple rate",strcat('(n=',n_chans,' channels)')})
    ylabel('Frequency (%)'), xlabel('Ripple rate (Hz)')

    subplot(2,2,2)
    histogram(allDur,'Normalization','probability'), hold on
    plot(repmat(median(allDur),length(0:15),1), (0:15)/100, 'r--','LineWidth',2), 
    plot(repmat(prctile(allDur,25),length(0:15),1), (0:15)/100, 'r--','LineWidth',1), 
    plot(repmat(prctile(allDur,75),length(0:15),1), (0:15)/100, 'r--','LineWidth',1), hold off
    axis square, yticklabels(yticks*100), title({"Duration of ripples",strcat('(n=',n_ripples,' ripples)')})
    ylabel("Frequency (%)"), xlabel("Duration (s)"), xlim([0 0.3])
    
    subplot(2,2,3)
    histogram(allMA,'Normalization','probability'), hold on
    plot(repmat(median(allMA),length(0:15),1), (0:15)/100, 'r--','LineWidth',2), 
    plot(repmat(prctile(allMA,25),length(0:15),1), (0:15)/100, 'r--','LineWidth',1), 
    plot(repmat(prctile(allMA,75),length(0:15),1), (0:15)/100, 'r--','LineWidth',1), hold off
    axis square, yticklabels(yticks*100), title({"Mean Amplitude of ripple",strcat('(n=',n_ripples,' ripples)')})
    ylabel("Frequency (%)"), xlabel("Mean Amplitude (Z-score)"), xlim([0 6])
    
    subplot(2,2,4)
    histogram(allPA,'Normalization','probability'), hold on
    plot(repmat(median(allPA),length(0:30),1), (0:30)/100, 'r--','LineWidth',2), 
    plot(repmat(prctile(allPA,25),length(0:30),1), (0:30)/100, 'r--','LineWidth',1), 
    plot(repmat(prctile(allPA,75),length(0:30),1), (0:30)/100, 'r--','LineWidth',1), hold off
    axis square, yticklabels(yticks*100), title({"Maximum Amplitude of ripple",strcat('(n=',n_ripples,' ripples)')})
    ylabel("Frequency (%)"), xlabel("Maximum Amplitude (Z-score)"), xlim([2 10])
end

%% Plot Cross-correlation of ripples
if plotXCorr
    allXCorr = [];
    allACorr = [];
    for file_no = 1:length(matFiles)
        % Load Output file
        load(char(fullfile(outputDir, matFiles(file_no).name)));
        
        % Extract XCorr
        for chan1 = 1:size(out.XCorr,1)-1
            for chan2 = 2:size(out.XCorr,2)
                if chan1~=chan2 && ~isnan(out.RippleFreq(chan1)) && ~isnan(out.RippleFreq(chan2))
                    allXCorr = [allXCorr [out.XCorr{chan1,chan2}]];
                end
            end, clear chan2
        end, clear chan1
        
        % Extract AutoCorr
        for chan = 1:size(out.XCorr,1)
            if ~isnan(out.RippleFreq(chan))
                allACorr = [allACorr [out.XCorr{chan,chan}]];
            end
        end, clear chan
    end, clear file_no
    
    % Set correlations at 0 lag = nan
    allXCorr(ceil(size(allXCorr,1)/2),:) = nan;
    allACorr(ceil(size(allACorr,1)/2),:) = nan;
    
    n_combiX = num2str(size(allXCorr,2));
    n_combiA = num2str(size(allACorr,2));
    
    figure
    subplot(1,2,1)
    shadedErrorBar(out.corrLags,mean(allXCorr,2),std(allXCorr,[],2)./sqrt(size(allXCorr,2)),'lineProps','b'), hold on
%     plot(out.corrLags,mean(allXCorr,2),'b'), hold on
    plot(zeros(2,1), [0 2.5e-6], 'r--'), hold off
    title({'Ripple cross-channel correlation',strcat('(n=',n_combiX,' combinations)')})
    xlabel('Time lag (s)'), ylabel({'X-correlation','(unbiased)'}), xlim([-250 250]), ylim([0 2.5e-6])
    axis square, xticklabels(xticks/1000)
    subplot(1,2,2)
    shadedErrorBar(out.corrLags,mean(allACorr,2),std(allACorr,[],2)./sqrt(size(allACorr,2)),'lineProps','m'), hold on
%     plot(out.corrLags,mean(allACorr,2),'m'), hold on
    plot(zeros(2,1), [0 2.5e-6], 'r--'), hold off
    title({'Ripple channel auto-correlation',strcat('(n=',n_combiA,' combinations)')})
    xlabel('Time lag (s)'), ylabel({'Auto-correlation','(unbiased)'}), xlim([-250 250]), ylim([0 2.5e-6])
    axis square, xticklabels(xticks/1000)
end

%% Plot cross-correlations between spikes and ripples
if plotSpikeCorr
    allSameChan = [];
    allDiffChan = [];

    for file_no = 1:length(matFiles)
        % Load spike output file
        load(char(fullfile(outputDir, matFiles(file_no).name)));
        load(char(fullfile(spikeOutDir,strcat(out.file,'-OutputWSpike.mat'))));
        load(char(fullfile(spikeOutDir,strcat(out.file,'-LabelledSpikes.mat'))));
        load(char(fullfile(inputDir,strcat(out.file,'-Input.mat'))));

        for cellNo = 1:length(spikeDat.times)
            spikes      = zeros(round(in.downFs*out.recordingDur),1);
            % Calculate index of spike times
            spikeIdx    = round(spikeDat.times{cellNo}*in.downFs);
            % Indicate when the spikes are in spikes vector
            spikes(spikeIdx) = 1;
            
            % Cycle through XCorr function with each channel
            for c = 1 : length(out.Channels)
                % Ripple vector
                ripplehist = zeros(round(out.recordingDur*in.downFs),1);
                % Append ripple start times to array for xcorr
                ripplehist(round([out.Channels{c}(:).StartTime].*in.downFs)) = 1;
                
                % Find X-corr
                r = xcorr(spikes, ripplehist, 'unbiased', round(in.plotWin(2)*in.downFs));
                if spikeDat.chanClustID(cellNo,1)==in.chans2keep(c)
                    % channel is the same as cell
                    allSameChan = [allSameChan r];
                else % other channels
                    allDiffChan = [allDiffChan r];
                end
            end, clear c
        end, clear cellNo
    end, clear file_no

    n_combiDiffSpikes = num2str(size(allDiffChan,2));
    n_combiSameSpikes = num2str(size(allSameChan,2));
    corrLags = -round(in.plotWin(2)*in.downFs):round(in.plotWin(2)*in.downFs);

    figure
    subplot(1,2,1)
    shadedErrorBar(corrLags,flip(mean(allDiffChan,2)),flip(std(allDiffChan,[],2))./sqrt(size(allDiffChan,2)),'lineProps','b'), hold on
    plot(zeros(2,1), [1e-6 9e-6], 'r--'), hold off
    title({'Spike cross-correlation with SWRs','(different channel)',strcat('(n=',n_combiDiffSpikes,' combinations)')})
    xlabel('Time lag (s)'), ylabel({'X-correlation','(unbiased)'}), xlim([-250 250]), ylim([1e-6 9e-6])
    axis square, xticklabels(xticks/1000)
    subplot(1,2,2)
    shadedErrorBar(corrLags,flip(mean(allSameChan,2)),flip(std(allSameChan,[],2))./sqrt(size(allSameChan,2)),'lineProps','m'), hold on
    plot(zeros(2,1), [1e-6 9e-6], 'r--'), hold off
    title({'Spike cross-correlation with SWRs','(same channel)',strcat('(n=',n_combiSameSpikes,' combinations)')})
    xlabel('Time lag (s)'), ylabel({'X-correlation','(unbiased)'}), xlim([-250 250]), ylim([1e-6 9e-6])
    axis square, xticklabels(xticks/1000)

end

%% Extract phases of ripples
if plotPhases
    allSameChanPhase = {};
    allDiffChanPhase = {};
    for file_no = 1:length(matFiles)
        % Load output file
        load(char(fullfile(outputDir, matFiles(file_no).name)));
        load(char(fullfile(spikeOutDir,strcat(out.file,'-OutputWSpike.mat'))));
        load(char(fullfile(spikeOutDir,strcat(out.file,'-LabelledSpikes.mat'))));
        load(char(fullfile(inputDir,strcat(out.file,'-Input.mat'))));
        
        for c = 1:length(out.Channels)
            if ~isnan(out.RippleFreq(c))
                % Extract phase
                chanPhase = vertcat(out.Channels{c}.PhaseOutput);
                cellSameChanPhase = [];
                cellDiffChanPhase = [];
                for cellNo = 1:size(chanPhase,2)
                    if spikeDat.chanClustID(cellNo,1)==in.chans2keep(c)
                        cellSameChanPhase = [cellSameChanPhase,chanPhase(:,cellNo)];
                    else 
                        cellDiffChanPhase = [cellDiffChanPhase,chanPhase(:,cellNo)];
                    end
                end, clear cellNo
                if ~isempty(cellSameChanPhase)
                    allSameChanPhase = [allSameChanPhase {cellSameChanPhase}];
                end
                allDiffChanPhase = [allDiffChanPhase {cellDiffChanPhase}];
            end
        end, clear c
    end, clear file_no
end