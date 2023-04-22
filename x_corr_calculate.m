%% Calculation and plotting of number of ripples/spikes within specified window of other ripples

outputDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Output';
spikeDataDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Spike Data';
inputDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Input';
corrWin   = 0.1;    % Finding other ripples within 0.1s (100ms) of ripples
ripCorr   = true;  % Toggle for finding ripple cross correlation
spikeCorr = true;  % Toggle for finding spike-ripple cross correlation

% Find all mat files in specified output folder.
matFilePattern = fullfile(outputDir, '*.mat');
matFiles       = dir(matFilePattern);
clear matFilePattern

if ripCorr
    allSameChanRipplesNo = 0;
    allSameChanRipples   = [];
    allDiffChanRipplesNo = 0;
    allDiffChanRipples   = [];
    
    for file_no = 1:length(matFiles)
        % Load Output file
        load(char(fullfile(outputDir, matFiles(file_no).name)));
        for c = 1:length(out.Channels)
            if ~isnan(out.RippleFreq(c))
                ChanRipples = zeros(length(out.Channels{c}),length(out.Channels));
        
                for r = 1:length(out.Channels{c})
                    for cc = 1:length(out.Channels)
                        if ~isnan(out.RippleFreq(cc))
                            % Find number of ripples in the same chan within corrWin seconds of current ripple
                            ChanRipples(r,cc) = sum([out.Channels{cc}.PeakTime]>out.Channels{c}(r).PeakTime-corrWin & [out.Channels{cc}.PeakTime] < out.Channels{c}(r).PeakTime+corrWin);
                        end
                    end, clear cc
                end, clear r
                
                % Save number of ripples in that channel that have at least one other ripple, either in same chan or other chan
                allSameChanRipplesNo    = allSameChanRipplesNo + nnz(ChanRipples(:,c)-1);
                allSameChanRipples      = vertcat(allSameChanRipples, ChanRipples(:,c)-1);
                ChanRipples(:,c)        = [];                                   % Remove column of the same channel
                allDiffChanRipplesNo    = allDiffChanRipplesNo + nnz(any(ChanRipples,2));
                allDiffChanRipples      = vertcat(allDiffChanRipples, sum(ChanRipples,2)./size(ChanRipples,2));
    
            end
        end, clear c
    end, clear file_no
    
    n_ripples = num2str(length(allDiffChanRipples));
    
    figure
    subplot(1,2,1)
    edges = linspace(0,1,11);
    histogram(allDiffChanRipples,edges,'Normalization','probability','FaceColor','b'), % hold on
%     plot(repmat(median(allDiffChanRipples),1,length(0:40)), (0:40)/100, 'r--','LineWidth',2), 
%     plot(repmat(prctile(allDiffChanRipples,25),length(0:40),1), (0:40)/100, 'r--','LineWidth',1), 
%     plot(repmat(prctile(allDiffChanRipples,75),length(0:40),1), (0:40)/100, 'r--','LineWidth',1), hold off
    axis square, yticklabels(yticks*100), title({"Number of SWRs","(different channel)",strcat('(n=',n_ripples,' ripples)')})
    ylabel('Frequency (%)'), xlabel('Number of SWRs (average)'), xlim([0 1])
    
    subplot(1,2,2)
    histogram(allSameChanRipples,'Normalization','probability','FaceColor','m'), % hold on
    % plot(repmat(median(allSameChanRipples),length(0:15),1), (0:15)/100, 'r--','LineWidth',2), 
    % plot(repmat(prctile(allSameChanRipples,25),length(0:15),1), (0:15)/100, 'r--','LineWidth',1), 
    % plot(repmat(prctile(allSameChanRipples,75),length(0:15),1), (0:15)/100, 'r--','LineWidth',1), hold off
    axis square, yticklabels(yticks*100), title({"Number of SWRs","(same channel)",strcat('(n=',n_ripples,' ripples)')})
    ylabel("Frequency (%)"), xlabel("Number of SWRs (average)"), xticklabels({'0','1','2','3',''})
end

if spikeCorr
    allSameChanSpikesNo = 0;
    allSameChanSpikes   = [];
    allDiffChanSpikesNo = 0;
    allDiffChanSpikes   = [];

    for file_no = 1:length(matFiles)
        % Load Output file
        load(char(fullfile(outputDir, matFiles(file_no).name)));
        load(char(fullfile(spikeDataDir,strcat(out.file,'-spikeData.mat'))));
    
        for c = 1:length(out.Channels)
            if ~isnan(out.RippleFreq(c))
                ChanSpikes = zeros(length(out.Channels{c}),length(spikeDat.times));
        
                for r = 1:length(out.Channels{c})
                    for cc = 1:length(spikeDat.times)
                        
                            % Find number of spikes within corrWin seconds of current ripple
                            ChanSpikes(r,cc) = sum(spikeDat.times{cc}>out.Channels{c}(r).PeakTime-corrWin & spikeDat.times{cc} < out.Channels{c}(r).PeakTime+corrWin);
                        
                    end, clear cc
                end, clear r
                
                % Save number of ripples in that channel that have at least one spike, either in same chan or other chan
                load(char(fullfile(inputDir,strcat(out.file,'-Input.mat'))));

                allSameChanSpikesNo     = allSameChanSpikesNo + nnz(any(ChanSpikes(:,1:size(ChanSpikes,2)==in.chans2keep(c)),2));
                allSameChanSpikes       = vertcat(allSameChanSpikes, sum(ChanSpikes(:,1:size(ChanSpikes,2)==in.chans2keep(c)),2)./size(ChanSpikes(:,1:size(ChanSpikes,2)==in.chans2keep(c)),2));
                ChanSpikes(:,1:size(ChanSpikes,2)==in.chans2keep(c)) = [];  % Remove column(s) of the same channel
                allDiffChanSpikesNo     = allDiffChanSpikesNo + nnz(any(ChanSpikes,2));
                allDiffChanSpikes       = vertcat(allDiffChanSpikes, sum(ChanSpikes,2)./size(ChanSpikes,2));

            end
        end, clear c
    end, clear file_no
    
    % Remove nans from allSameChanSpikes
    allSameChanSpikes = allSameChanSpikes(~isnan(allSameChanSpikes));

    n_ripples_diff = num2str(length(allDiffChanSpikes));
    n_ripples_same = num2str(length(allSameChanSpikes));
    
    figure
    subplot(1,2,1)
    edges = linspace(0,10,11);
    histogram(allDiffChanSpikes,edges,'Normalization','probability','FaceColor','b'), % hold on
%     plot(repmat(median(allDiffChanSpikes),length(0:60),1), (0:60)/100, 'r--','LineWidth',2), 
%     plot(repmat(prctile(allDiffChanSpikes,25),length(0:60),1), (0:60)/100, 'r--','LineWidth',1), 
%     plot(repmat(prctile(allDiffChanSpikes,75),length(0:60),1), (0:60)/100, 'r--','LineWidth',1), hold off
    axis square, yticklabels(yticks*100), title({"Number of spikes","(different channel)",strcat('(n=',n_ripples_diff,' ripples)')})
    ylabel('Frequency (%)'), xlabel('Number of spikes (average)'), ylim([0 0.70])
    
    subplot(1,2,2)
    histogram(allSameChanSpikes,edges,'Normalization','probability','FaceColor','m'), % hold on
    % plot(repmat(median(allSameChanRipples),length(0:15),1), (0:15)/100, 'r--','LineWidth',2), 
    % plot(repmat(prctile(allSameChanRipples,25),length(0:15),1), (0:15)/100, 'r--','LineWidth',1), 
    % plot(repmat(prctile(allSameChanRipples,75),length(0:15),1), (0:15)/100, 'r--','LineWidth',1), hold off
    axis square, yticklabels(yticks*100), title({"Number of spikes","(same channel)",strcat('(n=',n_ripples_same,' ripples)')})
    ylabel("Frequency (%)"), xlabel("Number of spikes (average)"), ylim([0 0.70])
end