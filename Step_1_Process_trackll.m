% trackll .csv to trackerPar .mat conversion
% written by Sun Jay Yoo
% original anisotropy code by Anders Sejr Hansen
% [https://gitlab.com/anders.sejr.hansen/anisotropy]

data_path = ['.', filesep, 'Input', filesep];
output_path = ['.', filesep, 'trackerPar_reformatted', filesep];
%filename = "TFIIA_COMP3.csv";
frame_rate = 100; %Hz 
resolution = 16/150;

%%%%%%%%%%%%%%%%%%%%%%%% STEP 1: DEFINE DATASETS %%%%%%%%%%%%%%%%%%%%%%%%%%

% run on all CSV-files in the input_path directory
CSV_files=dir([data_path,'*.csv']);
FileNames = ''; %for saving the actual file name

for iter = 1:length(CSV_files)
    % get the relevant file name
    FileNames{iter} = CSV_files(iter).name(1:end-4);
end

%%%%%%%%%%%%%%%%%%% STEP 2: REFORMAT TO Anders FORMAT %%%%%%%%%%%%%%%%%%%%
%   loop over all CSV_files
%   convert to Anders format
%   save file in Anders format
disp(['Reformatting ', num2str(length(FileNames)), ' datasets to Anders format']); tic;

for FileIter = 1:length(FileNames)
    disp(['working on file ', num2str(FileIter), ' of ', num2str(length(FileNames)), ' total files']);
    % load the data
    clear dataset_1
    clear data_1
    dataset_1 = strcat(data_path, CSV_files(FileIter).name(1:end-4), '.csv');
    data_1 = readtable(dataset_1);
    data_1 = removevars(data_1, {'Var1'});
    data_1 = removevars(data_1, {'z'});
    trackedPar_QC = struct('xy',{},'Frame',{},'TimeStamp',{});
    struct_row = 1;
    trajectory_start = 1;
    for i = 1:size(data_1,1)
        if i == size(data_1,1)
        trackedPar_QC(struct_row).xy = table2array(data_1(trajectory_start:i,3:4)).*resolution;
        trackedPar_QC(struct_row).Frame = table2array(data_1(trajectory_start:i,2));
        trackedPar_QC(struct_row).TimeStamp = trackedPar_QC(struct_row).Frame/frame_rate;
        elseif struct_row ~= data_1.Trajectory(i)
        trackedPar_QC(struct_row).xy = table2array(data_1(trajectory_start:i-1,3:4)).*resolution;
        trackedPar_QC(struct_row).Frame = table2array(data_1(trajectory_start:i-1,2));
        trackedPar_QC(struct_row).TimeStamp = trackedPar_QC(struct_row).Frame/frame_rate;
        trajectory_start = i;
        struct_row = struct_row + 1;
        end
    end
    save([output_path,FileNames{FileIter},'_',num2str(frame_rate),'Hz.mat'], 'trackedPar_QC')
end
toc;