%% Wrapper script for ripple_analysis.m

%% Specify folders where data is located. Wrapper will loop through all .ns6 files in folder.
NS6Folder = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse';
InputFolder = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Input';
OutputFolder = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Output';
TaskFolder = 'C:\Users\chick\Desktop\Ripple analysis\Sample data to analyse\Microwire Task Data';

% Find all NS6 files in NS6folder.
NS6FilePattern = fullfile(NS6Folder, '*.ns6');
NS6Files = dir(NS6FilePattern);
clear NS6FilePattern NS6Folder

% Specify whether to run morletTF and XCorr
in.TFXCToggle = true;

%% Run ripple_analysis through all NS6 files in NS6folder.
for file_no = 1:length(NS6Files)
    
    % Directory for raw ns6file
    in.dataFile  = fullfile(NS6Files(file_no).folder, NS6Files(file_no).name);
    [~, name, ~] = fileparts(in.dataFile);

    % Load associated Task Data variable
    load(fullfile(TaskFolder, strcat(name,'-TaskData.mat')));

    % Specify channels to keep
    in.chans2keep = arts.channels;

    % Specify timepoints to exclude
    in.exclude = [arts.iedTimes{:}]; clear arts

    disp(strcat('Currently running: ',name))
    
    % Run analysis
    [in,out] = ripple_analysis(in);
    out.file = name;
    
    % Save input and output
    save(fullfile(OutputFolder, strcat(name, '-Output')), 'out', '-v7.3');
    save(fullfile(InputFolder, strcat(name, '-Input')), 'in', '-v7.3');

end, clear in file_no name NS6Files TaskFolder OutputFolder out InputFolder
