function [in,out] = ripple_analysis(in)
%% Function to import raw .ns6 datafile and identify ripples.    

    %% Specify default parameters (as fields of the 'in' structure)

    if ~isfield(in,'dataFile')                                              % Directory for raw ns6file
        error('No input data file has been specified.'); end
    if ~isfield(in,'downFs');     in.downFs     = 1000; end                 % Downsample factor for resample function
    if ~isfield(in,'rippleBand'); in.rippleBand = [80 120]; end             % Target freq for SWR
    if ~isfield(in,'fOrder');     in.fOrder     = 5; end                    % Order for Butterworth filter
    if ~isfield(in,'chans2keep'); in.chans2keep = 1:2:31; end               % Change this to 1:2:31 for all channels
    if ~isfield(in,'mindur');     in.mindur     = 0.020; end                % Minimum duration for ripple
    if ~isfield(in,'zThresh');    in.zThresh    = 0; end                    % Z-threshold for ripple
    if ~isfield(in,'peakZ');      in.peakZ      = 3; end                    % Cut off for ripple peak
    if ~isfield(in,'morletFreq'); in.morletFreq = 1:3:200; end              % Frequency for morletTF input
    if ~isfield(in,'startEndEx'); in.startEndEx = 1; end                    % Duration (s) to exclude at start and end of recording
    if ~isfield(in,'exclude');    in.exclude    = []; end                   % Timepoint to exclude (artefacts)
    if ~isfield(in,'exWin');      in.exWin      = 1; end                    % +-Window to exclude ripples around artefacts
    if ~isfield(in,'plotWin');    in.plotWin    = [-0.25 0.25]; end         % Plot window for ripples +-1s
    if (~isfield(in,'saveBuffer') || in.saveBuffer<in.startEndEx)
                                  in.saveBuffer = in.startEndEx; end        % +-Window for saving ripple data (minumum = in.startEndEx)
    if ~isfield(in,'TFXCToggle'); in.TFXCToggle = false; end                % Toggle for enabling running of morletTF and XCorr
    if ~isfield(in,'minRipple');  in.minRipple  = 1; end                    % Minimum ripple rate (per min) for channels to be included in analysis

    
    % Pre-define cell for channel output
    out.Channels        = cell(1,length(in.chans2keep));
    out.ExcludedRipples = cell(1,length(in.chans2keep));
    out.RippleDur       = nan(length(in.chans2keep),1);
    out.RippleAmp       = nan(length(in.chans2keep),1);
    
    %% Identification and analysis of ripples
    % Load the raw data
    sample = openNSx(in.dataFile,'uV'); 
    
    % Extract sampling frequency
    samplerate = sample.MetaTags.SamplingFreq;

    % Extract sample duration
    tDur = sample.MetaTags.DataDurationSec;
    
    % Generate a band-pass filter
    [b,a] = butter(in.fOrder, in.rippleBand./(in.downFs/2));
    
    % Cycle through each channel
    for c = 1 : length(in.chans2keep)
        
        % Extract relevant channel based on channel ID
        sampledata = double(sample.Data([sample.ElectrodesInfo.ElectrodeID]==in.chans2keep(c),:));
        
        % Downsample the signal
        % Transpose data due to functions working in columns
        dsdata = resample(sampledata',in.downFs,samplerate); clear sampledata
        
        % Filter in ripple band
        filterdata = filtfilt(b,a,dsdata);
        
        % Hilbert transform to extract amplitude and phase
        hdata = abs(hilbert(filterdata)); 
        pdata = angle(hilbert(filterdata));

        % Generate a mask around all artefacts
        artMask = true(size(hdata));
        artWin  = -in.exWin*in.downFs : in.exWin*in.downFs;
        artInd  = round(in.exclude*in.downFs);
        if size(artInd,1) == 1
            artInd = artInd';
        end
        artInd  = repmat(artInd,1,length(artWin)) + repmat(artWin,length(artInd),1); clear artWin
        artInd(artInd<1) = 1;
        artInd(artInd>length(artMask)) = length(artMask);
        artInd  = unique(artInd);
        artMask(artInd) = false; clear artInd
        
        % Calculate z-score of filtered data with mean and std of non-masked data
        zsc = (hdata - mean(hdata(artMask))) ./ std(hdata(artMask)); clear artMask
        
        % Generate binary image of z-scores
        binImage = zsc >= in.zThresh;
        
        % Find ripples with duration > target duration (ie area > target duration)
        allRipples    = regionprops(binImage, zsc, "Area", "PixelIdxList","MaxIntensity","MeanIntensity","PixelValues");
        include       = [allRipples.Area] >= (in.mindur*in.downFs) & [allRipples.MaxIntensity] >= in.peakZ;
        targetRipples = allRipples(include); clear allRipples binImage include 
        
        % Find start and end times of ripples
        startIdx    = cellfun(@min,{targetRipples(:).PixelIdxList});
        startTime   = num2cell(startIdx./in.downFs);
        endIdx      = cellfun(@max,{targetRipples(:).PixelIdxList});
        endTime     = num2cell(endIdx./in.downFs);
        [targetRipples.StartTime] = startTime{:};
        [targetRipples.EndTime]   = endTime{:}; clear startTime endTime
        
        % Find peak times of ripples
        [~,peakIdx] = cellfun(@max,{targetRipples(:).PixelValues});
        peakIdx     = startIdx + peakIdx - 1;
        peakTime    = num2cell(peakIdx./in.downFs);
        [targetRipples.PeakTime] = peakTime{:}; clear peakTime peakIdx startIdx
        
        % Find ripple duration and midpoint
        duration = num2cell(([targetRipples.EndTime] - [targetRipples.StartTime])');
        midpoint = num2cell((([targetRipples.EndTime] + [targetRipples.StartTime])./2)');
        [targetRipples.Duration] = duration{:}; 
        [targetRipples.Midpoint] = midpoint{:}; clear duration midpoint
        
        % Exclude any ripples within specified seconds of the start or end of recording
        targetRipples([targetRipples.StartTime] < in.startEndEx | [targetRipples.EndTime] > (length(dsdata)/in.downFs-in.startEndEx),:) = [];
        
        % Output and exclude any ripples within exWin seconds of specified timepoint (artefacts)
        exStart     = in.exclude - in.exWin; 
        exEnd       = in.exclude + in.exWin;
        toExclude   = any([targetRipples.EndTime] > exStart(:) & [targetRipples.StartTime] < exEnd(:),1);
        out.ExcludedRipples{c} = targetRipples(toExclude,:);
        targetRipples(toExclude,:) = []; clear exStart exEnd toExclude

        % Output various attributes
        sampWin     = round(-in.saveBuffer*in.downFs):round(in.saveBuffer*in.downFs);
        sampBase    = cellfun(@(x) round(x*in.downFs)+sampWin,{targetRipples(:).PeakTime},'uni',false); clear sampWin
        timeBase    = cellfun(@(x) x./in.downFs,sampBase,'uni',false);
        rawSignal   = cellfun(@(x) dsdata(x),sampBase,'uni',false);
        filtSignal  = cellfun(@(x) filterdata(x),sampBase,'uni',false);
        amplitude   = cellfun(@(x) hdata(x),sampBase,'uni',false); 
        phase       = cellfun(@(x) pdata(x),sampBase,'uni',false); 
        zScore      = cellfun(@(x) zsc(x),sampBase,'uni',false); 
        peakzScore  = cellfun(@(x) x(in.saveBuffer*in.downFs+1),zScore,'uni',false); 

        [targetRipples.TimeBase]    = timeBase{:}; clear timeBase
        [targetRipples.RawSignal]   = rawSignal{:}; clear rawSignal
        [targetRipples.FiltSignal]  = filtSignal{:}; clear filtSignal
        [targetRipples.Amplitude]   = amplitude{:}; clear amplitude
        [targetRipples.Phase]       = phase{:}; clear phase
        [targetRipples.ZScore]      = zScore{:}; clear zScore
        [targetRipples.PeakZScore]  = peakzScore{:}; clear peakzScore
        
        % Save mean duration and amplitude
        out.RippleDur(c) = nanmean([targetRipples.Duration]);
        out.RippleAmp(c) = nanmean([targetRipples.PeakZScore]); 
        
        for i = 1:size(targetRipples,1)
            if in.TFXCToggle

                % Time frequency using morletTF - power output and average across ripples
                [pow, TFphase]           = morletTF(targetRipples(i).RawSignal', cell2mat(sampBase(i))./in.downFs, in.morletFreq);
                targetRipples(i).Power   = pow;
                targetRipples(i).TFPhase = TFphase;

            end
        end, clear pow phase BufferStartIdx BufferEndIdx i dsdata filterdata hdata zsc
        
        % Remove superfluous outputs
        targetRipples = rmfield(targetRipples,{'Area','PixelIdxList','PixelValues'});

        % Store data for that channel
        out.Channels{c} = targetRipples; clear targetRipples
        
    end, clear a b c samplerate
    
    % Store recording duration, compute ripple frequency for each channel
    out.recordingDur    = tDur;
    out.RippleFreq      = cellfun(@length,out.Channels) ./ tDur; clear sample
    out.RippleFreq(out.RippleFreq < (in.minRipple/60)) = nan;
    out.RippleDur(isnan(out.RippleFreq)) = nan;
    out.RippleAmp(isnan(out.RippleFreq)) = nan;
    out.avgRippleFreq   = nanmean(out.RippleFreq);
    out.avgRippleDur    = nanmean(out.RippleDur);
    out.avgRippleAmp    = nanmean(out.RippleAmp);

    %% Calculate Xcorr across channels
    
    if in.TFXCToggle

        % Append ripple start times to array for xcorr
        ripplehist = zeros(length(out.Channels),ceil(tDur*in.downFs));
        for c = 1 : length(out.Channels)
            ripplehist(c,ceil([out.Channels{c}(:).StartTime].*in.downFs)) = 1;
        end, clear c
        
        % Find xcorr
        [r, lags]    = xcorr(ripplehist', 'unbiased', ceil(0.5*in.downFs));
        out.corrLags = lags; clear lags
        
        % Note r(:,2) and r(:,n+1) are mirror images. Both are xcorr of C1/C2.
        % Include autocorrelation, ie each channel with itself
        % So, we keep columns 1:n, n+2:2n, 2n+3:3n,..., n^2
        n = length(out.Channels);
        corrCell = cell([n n]);
        
        for c = 0:n-1
            for i = (c*n+1+c):((c+1)*n)
                
                chan1 = c+1;
                chan2 = rem(i,n);
                if chan2 == 0; chan2 = n; end
                corrCell{chan1, chan2} = r(:,i);
        
            end, clear i chan1 chan2
        end, clear c r n
        
        % Save XCorr data
        out.XCorr = corrCell; clear corrCell
    end

end