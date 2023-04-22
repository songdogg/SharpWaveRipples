%% Wrapper script for task_analysis.m and spike_analysis.m

%% Specify input directories of files to analyse
anIn.rippleDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Output';
anIn.rippleInputDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Input';
anIn.taskDataDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Microwire Task Data';
anIn.outputDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Output\Analysis';
anIn.allPtsDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Output\Analysis\Across Pts';
anIn.spikeDataDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Spike Data';
anIn.spikeOutDir = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Spike Data\Spike Output';

%% Run task analysis
anIn = task_analysis(anIn);

%% Run spike analysis
anIn = spike_analysis(anIn);

%% Save anIn
save(fullfile(anIn.outputDir, 'TaskSpikeInput'), 'anIn', '-v7.3');