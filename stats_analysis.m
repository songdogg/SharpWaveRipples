%% Statistical anaylsis for comparisons between task periods

%% Load datafile and define key variable
load 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Output\Analysis\Across Pts\Analysis Across Patients.mat'

n_pts = length(allPts.Patients);

%% ANOVA

taskPhaseRippleFreqAnova    = anova_stats(allPts.AxChan,'rippleFreq');
taskPhaseMeanAmpAnova       = anova_stats(allPts.AxChan,'meanAmp');
taskPhaseMeanDurAnova       = anova_stats(allPts.AxChan,'meanDur');

Plot ANOVA bar chart
figure
subplot(3,3,1)  % Ripple frequency overall by phase
plot_anova(taskPhaseRippleFreqAnova,n_pts)
title('Ripple Rate'), ylabel('Rate (Hz)')

subplot(3,3,2)  % Mean amplitude overall by task phase
plot_anova(taskPhaseMeanAmpAnova,n_pts)
title('Mean Amplitude'), ylabel('Z-score')

subplot(3,3,3)  % Mean duration overall by task phase
plot_anova(taskPhaseMeanDurAnova,n_pts)
title('Mean Duration'), ylabel('Duration (s)')

%% Paired t-tests

[pairTRF,diffRF] = pairt_stats(allPts.AxChan,'rippleFreq');
[pairTMA,diffMA] = pairt_stats(allPts.AxChan,'meanAmp');
[pairTMD,diffMD] = pairt_stats(allPts.AxChan,'meanDur');

subplot(3,3,4)  % Ripple frequency: split by task characteristics
plot_pairt(diffRF,n_pts)
title('Ripple Rate'), ylabel('Rate Difference (Hz)')

subplot(3,3,5)  % Mean amplitude: split by task characteristics
plot_pairt(diffMA,n_pts)
title('Mean Amplitude'), ylabel('Z-score Difference')

subplot(3,3,6)  % Mean duration: split by task characteristics
plot_pairt(diffMD,n_pts)
title('Mean Duration'), ylabel('Duration Difference (s)')

%% t-test for time before task
ratioRF = ratiot_stats(allPts,'rippleFreq','overallAvgRippleFreq');
ratioMA = ratiot_stats(allPts,'meanAmp','overallAvgRippleAmp');
ratioMD = ratiot_stats(allPts,'meanDur','overallAvgRippleDur');

subplot(3,3,7)  % Ripple frequency: time before rec/ass1/2
plot_ratiot(ratioRF,n_pts)
title('Ripple Rate')

subplot(3,3,8)  % Mean amplitude: time before rec/ass1/2
plot_ratiot(ratioMA,n_pts)
title('Mean Amplitude')

subplot(3,3,9)  % Mean duration: time before rec/ass1/2
plot_ratiot(ratioMD,n_pts)
title('Mean Duration')

%% Functions used in this script
function anova_table = anova_stats(structure,field)    % field is a string
    anova_table = [structure.StudyBaseline.(field) ...
                   structure.StudyTask.(field) ...
                   structure.Consolidation.(field) ...
                   structure.TestBaseline.(field) ...
                   structure.TestTask.(field)];
end

function [pairt_table,diff] = pairt_stats(structure,field)
    pairt_table = [structure.StudyFirst.(field)     structure.StudySub.(field)...
                   structure.StudyTask0.(field)     structure.StudyTask2.(field)...
                   structure.TestTask0.(field)      structure.TestTask2.(field)...
                   structure.TestTaskOld.(field)    structure.TestTaskNew.(field)...
                   structure.TestTaskSlowRT.(field) structure.TestTaskFastRT.(field)];
    diff = nan(length(structure.StudyTask.(field)),5);
    for n = 1:5
        diff(:,n) = pairt_table(:,2*n) - pairt_table(:,2*n-1);
    end
end

function ratiot_table = ratiot_stats(structure,field,avgfield)
    ratiot_table = [structure.AxChan.TestBefRec.(field) ...
                    structure.AxChan.TestBefAss1.(field) ...
                    structure.AxChan.TestBefAss2.(field)];
    % Extract overall ripple rate from each patient
    avg = nan(length(structure.Patients),1);
    for pt = 1:length(structure.Patients)
        avg(pt) = structure.Patients{pt}.(avgfield);
    end, clear pt
    ratiot_table = ratiot_table ./ avg;
end

function plot_anova(anova_table,n_pts)
    bar(nanmean(anova_table)), hold on
    errorbar(nanmean(anova_table),nanstd(anova_table)./sqrt(n_pts),'k.')
    scatter(sort(repmat(1:5,1,size(anova_table,1))),anova_table(:),100,'k.'), hold off
    xticklabels({'Study Baseline','Study','Consolidation','Test Baseline','Test'})
end

function plot_pairt(pairt_diff,n_pts)
    bar(nanmean(pairt_diff)), hold on
    errorbar(nanmean(pairt_diff),nanstd(pairt_diff)./sqrt(n_pts),'k.')
    scatter(sort(repmat(1:5,1,size(pairt_diff,1))),pairt_diff(:),100,'k.'), hold off
    xticklabels({'Sub vs 1st (Study)','✓ vs ✗ (Study)','✓ vs ✗ (Test)','New vs Old (Test)','Fast vs Slow (Test)'})
end

function plot_ratiot(ratiot,n_pts)
    bar(nanmean(ratiot)), hold on
    errorbar(nanmean(ratiot),nanstd(ratiot)./sqrt(n_pts),'k.')
    scatter(sort(repmat(1:3,1,size(ratiot,1))),ratiot(:),100,'k.')
    plot(-1:5,ones(length(-1:5),1),'k-'), hold off
    xticklabels({'1s Bef Recog','1s Bef Assoc 1','1s Bef Assoc 2'})
    ylabel('Ratio'), xlim([0 4])
end