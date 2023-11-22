%   Batch_vbSPT_classify.m
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3
%   Dependent functions:
%       - vbSPT-1.1.3 and all associated functions.m
%       - InferFrameRateFromName.m
%       - EditRunInputFile_for_batch.m
clear; clc; close all;

%   DESCRIPTION
%   This file takes as input a list of MAT-files with SPT data, converts
%   them to the vbSPT readable format and then classifies them according to
%   the vbSPT protocol
%   vbSPT is a sophisticated HMM-implementation for SPT data. Please see
%   the citation below for full details on vbSPT: 
%   Person F,Lindén M, Unoson C, Elf J, Extracting intracellular reaction rates from single molecule tracking data, Nature Methods 10, 265?269 (2013). doi:10.1038/nmeth.2367

%%%%%%%%%%%%%%%%%%%%%%%% PREAMBLE %%%%%%%%%%%%%%%%%%%%%%%%
% add path for dependent functions:
addpath('Functions');
% Add paths that are needed to run the VB3 analysis 
dir0= ['.', filesep, 'Functions', filesep, 'HMM_classification', filesep, 'vbSPT-1.1.3'];
addpath(genpath([dir0 filesep '.' filesep 'VB3']))
addpath(genpath([dir0 filesep '.' filesep 'HMMcore']))
addpath(genpath([dir0 filesep '.' filesep 'Tools']))
addpath(genpath([dir0 filesep '.' filesep 'external']))
addpath([dir0 filesep '.'])
disp('Added local vbSPT paths');
disp('------------------------------------');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   DEFINE KEY PARAMETERS
MinTrajLength = 2;
maxHidden = 3;

input_path = ['.', filesep, 'trackerPar_reformatted', filesep];
path_reformatted = ['.', filesep, 'QC_data_reformatted', filesep];
path_classified = ['.', filesep, 'vbSPT_classified', filesep];

%%%%%%%%%%%%%%%%%%%%%%%% STEP 1: DEFINE DATASETS %%%%%%%%%%%%%%%%%%%%%%%%%%

% run on all MAT-files in the input_path directory
MAT_files = dir([input_path,'*.mat']);
FileNames = ''; %for saving the actual file name

for iter = 1:length(MAT_files)
    % get the relevant file name
    Filenames{iter} = MAT_files(iter).name(1:end-4);
end



%%%%%%%%%%%%%%%%%%% STEP 2: REFORMAT TO vbSPT FORMAT %%%%%%%%%%%%%%%%%%%%

%   loop over all MAT_files
%   convert to vbSPT format
%   save file in vbSPT format
disp(['Reformatting ', num2str(length(Filenames)), ' datasets to the vbSPT format']); tic;
for FileIter = 1:length(Filenames)
    disp(['working on file ', num2str(FileIter), ' of ', num2str(length(Filenames)), ' total files']);
    % load the data
    clear trackedPar_QC CellTracks
    load([input_path, Filenames{FileIter}, '.mat']);
    
    % CONVERT TO vbSPT FORMAT
    % make a cell array of trajectpries
    CellTracks = cell(1);
    iter = 1;
    
    % REMOVE GAPS FROM TRAJECTORIES: vbSPT cannot handle gaps
    for i = 1:length(trackedPar_QC)
        if size(trackedPar_QC(i).xy,1) >= MinTrajLength
            %Now slice up trajectory to get rid of gaps
            currFrame = trackedPar_QC(i).Frame;
            currTrack = trackedPar_QC(i).xy;
            %check for gaps - expect a difference of 1 between frames if there
            %are no gaps, so search for difference of 2:
            Idx = find(diff(currFrame)==2);
            %If there were no gaps, just save:
            if isempty(Idx)           
                CellTracks{iter} = currTrack;
                iter = iter + 1;
            else
                i
                %Modify Idx:
                Idx
                length(currFrame)
                Idx_slice = [0 Idx length(currFrame)];
                %So there are gaps:
                %Slice up trajectories if they satisfy MinTrajLength 
                for j=2:length(Idx_slice)
                    if size(currTrack(Idx_slice(j-1)+1:Idx_slice(j),:),1) >= MinTrajLength
                        CellTracks{iter} = currTrack(Idx_slice(j-1)+1:Idx_slice(j),:);
                        iter = iter + 1;
                    end
                end
            end
            
        end
    end
    clear currFrame currTrack iter i Idx Idx_slice
    
    % OK now you have a cell array in the vbSPT format, so now you can save
    % the file, but first get the frame rate 
    LagTime  = InferFrameRateFromName( Filenames{FileIter} );
    
    % SAVE TRAJECTORIES IN THE vbSPT REFORMATTED FORM
    % save the reformatted dataset
    save([path_reformatted,Filenames{FileIter}, '_reformatted.mat'], 'CellTracks', 'LagTime');
end
toc;


%%%%%%%%%%%%%%%%% STEP 3: PERFORM vbSPT CLASSIFICATION %%%%%%%%%%%%%%%%%%%%

%   loop over all MAT_files

%   edit the "vbSPT_RunInputFileBatch" according to files names, paths and
%   frame rates

%   Perform vbSPT classification

%   save the output


disp(['Performing vbSPT classification of ', num2str(length(Filenames)), ' datasets']); 
for FileIter = 1:length(Filenames)
    disp(['working on HMM-classification of file number ', num2str(FileIter), ' of ', num2str(length(Filenames)), ' total files']);
    % load the data
    clear LagTime CellTracks CellTrackViterbiClass
    load([path_reformatted,Filenames{FileIter}, '_reformatted.mat']);
    
    % UPDATE 2017-10-09
    % I am getting a lot of weird conflicts, where iteration n, get's parts
    % of its information from the RunInputFile from the previous run. So to
    % deal with this, make a new RunInputFile for each run. 
    % Copy the master-file to directory 'FolderWithRunInputFiles' and then
    % give each file a new number. 
    
    % first see if the file already exists:
    if exist(['./FolderWithRunInputFiles/vbSPT_RunInputFileBatch', num2str(FileIter), '.m'])
        delete(['./FolderWithRunInputFiles/vbSPT_RunInputFileBatch', num2str(FileIter), '.m']);
    end
    % now edit the RunInputFileBatch    
    % EDIT THE RunInputFile 
    inputfile = [path_reformatted,Filenames{FileIter}, '_reformatted.mat'];
    outputfile = [path_classified,Filenames{FileIter}, '_classified.mat'];
    % now edit the RunInputFile using a function
    Finish = EditRunInputFile_for_batch( 'vbSPT_RunInputFileBatch.m', inputfile, outputfile, LagTime );
    
    % now copy it to the directory:
    copyfile('vbSPT_RunInputFileBatch.m', ['./FolderWithRunInputFiles/vbSPT_RunInputFileBatch', num2str(FileIter), '.m']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% Run vbSPT on the data %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R=VB3_HMManalysis(['./FolderWithRunInputFiles/vbSPT_RunInputFileBatch', num2str(FileIter), '.m']);
    %Now distill out the key things:
    CellTrackViterbiClass = R.Wbest.est2.viterbi;

    %Display the output:
    disp('====================================================');
    %VB3_getResult('vbSPT_RunInputFileBatch');
    
    % save the key vbSPT results to a metadata structure array
    vbSPT_metadata = struct();
    vbSPT_metadata.NumberOfStates = R.Wbest.N;
    vbSPT_metadata.DiffusionConstants = R.Wbest.est.DdtMean/LagTime;
    vbSPT_metadata.Occupancy = R.Wbest.est.Ptot;
    vbSPT_metadata.TransitionMatrix = R.Wbest.est.Amean;
    vbSPT_metadata.LagTime = LagTime;
    
    % re-load CellTracks:
    load([path_reformatted,Filenames{FileIter}, '_reformatted.mat'], 'CellTracks');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% SAVE THE DATA %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save([path_classified,regexprep(Filenames{FileIter}, {'_100Hz'}, {''}), ...
        '_classified.mat'], 'CellTracks', 'CellTrackViterbiClass', 'vbSPT_metadata', 'LagTime');
end

% Clear 'QC_data' and 'QC_data_reformatted' folders
% delete([input_path,'*.mat']);
delete([path_reformatted,'*.mat']);
delete([path_classified,'*_100Hz_classified.log']);
delete([path_classified,'*_100Hz_classified.mat']);
delete([input_path,'*.mat']);


