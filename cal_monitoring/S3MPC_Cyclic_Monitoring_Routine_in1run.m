function S3MPC_Cyclic_Monitoring_Routine_in1run(mission, cycle)
% All monitoring routines for Cyclic Reports in one run
% Includes:
% - Downloading from ftp
% - Monitoring CAL L1B, Thermal and USO for the whole mission
% - Monitoring CAL L1B and Thermal for one cycle

%% Downloading files from ftp site
USOdir=S3MPCftp_CALdownload(mission, cycle); % USO folder is output for USO monitoring routine

%%
% Initializing folders locations
switch mission
    case 'S3A'
        folderL0  = 'C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_0_CAL\';
        folderL1  = 'C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_1_CAL\CYCLE_ALL\';
    case 'S3B'
        folderL0  = 'C:\Users\Pablo\Desktop\S3MPC\S3B\data\CAL\SR_0_CAL\';
        folderL1  = 'C:\Users\Pablo\Desktop\S3MPC\S3B\data\CAL\SR_1_CAL\CYCLE_ALL\';
end

%%
%%%%%%%%---- mission monitoring ----%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Monitor CAL1 L1b data %%%%%%%%%%%%%%%%%%%%
disp(datetime)
tic
disp(['_______ Monitoring CAL1 L1b for WHOLE MISSION _______ folder ' folderL1]);
switch mission
    case 'S3A'
        Monitoring_CAL_S3A(folderL1)  % run S3A L1b monitoring routine for the whole mission
    case 'S3B'
        Monitoring_CAL_S3B(folderL1)  % run S3B L1b monitoring routine for the whole mission
end
toc

%%%%%%%%%%%%%%%%%%%% Monitor Thermal data %%%%%%%%%%%%%%%%%%%%
disp(datetime)
tic
disp(['_______ Monitoring CAL1 Thermal for WHOLE MISSION _______ folder ' folderL0]);
switch mission
    case 'S3A'
        S3A_Thermal(0, cycle) % run S3A Thermal monitoring routine for the whole mission
    case 'S3B'
        S3B_Thermal(0, cycle) % run S3B Thermal monitoring routine for the whole mission
end
toc

% Copy mission results in new corresponding folder
switch mission
    case 'S3A'
        refname=dir([folderL1 'S3MPC_S3A_CAL1_SAR_NORM_data_*']); % get the reference for dates from a file.
        mkdir([folderL1(1:end-10) 'whole_mission\' refname.name(30:60) '_BOMtoC0' num2str(cycle) '\']); % create the new folder
        movefile([folderL1 'S3MPC_S3A_CAL*'],[folderL1(1:end-10) 'whole_mission\' refname.name(30:60) '_BOMtoC0' num2str(cycle) '\']); % copy the files to the new folder
    case 'S3B'
        refname=dir([folderL1 'S3MPC_S3B_CAL1_SAR_NORM_data_*']); % get the reference for dates from a file.
        mkdir([folderL1(1:end-10) 'whole_mission\' refname.name(30:60) '_endBOM_C0' num2str(cycle) '\']); % create the new folder
        movefile([folderL1 'S3MPC_S3B_CAL*'],[folderL1(1:end-10) 'whole_mission\' refname.name(30:60) '_endBOM_C0' num2str(cycle) '\']); % copy the files to the new folder
end

%%%%%%%%%%%%%%%%%%%% Monitor USO %%%%%%%%%%%%%%%%%%%%
filein=dir([USOdir '\*.DBL']);
S3AB_USO_Monitor([USOdir '\' filein.name],mission) % run USO monitoring routine for the whole mission


%%
%%%%%%%%---- cycle monitoring ----%%%%%%%%

disp('Selecting L1 files from the cycle');

% list cycle L1 files ------------------------------------------
listL1=dir([folderL1 'cal_ncread_S3*.mat']);
for i=length(listL1):-1:1
    if str2double(listL1(i,1).name(81:83)) ~= cycle  % to delete
        listL1(i,:) = [];
    end
end
% list cycle L1 folders -----------------------------------------
listL1f=dir([folderL1 mission '_SR_1_CAL*']);
for i=length(listL1f):-1:1
    if str2double(listL1f(i,1).name(70:72)) ~= cycle  % to delete
        listL1f(i,:) = [];
    end
end

%%%%%%%%%%%%%%%%%%%% Monitor CAL1 L1b data %%%%%%%%%%%%%%%%%%%%
disp(datetime)
tic
disp(['_______ Monitoring CAL1 L1b for CYCLE ' num2str(cycle) '. Mission ' mission '.']);
mkdir([folderL1(1:end-9) num2str(cycle,'%.3d')]); % create cycle folder
for i=1:length(listL1)
    movefile([folderL1 listL1(i,1).name],[folderL1(1:end-9) num2str(cycle,'%.3d') '\']); % move cycle L1 files
end
for i=1:length(listL1f)
    movefile([folderL1 listL1f(i,1).name],[folderL1(1:end-9) num2str(cycle,'%.3d') '\']); % move cycle L1 folders
end
switch mission
    case 'S3A'
        Monitoring_CAL_S3A([folderL1(1:end-9) num2str(cycle,'%.3d') '/']) % run S3A L1b monitoring routine for the cycle
    case 'S3B'
        Monitoring_CAL_S3B([folderL1(1:end-9) num2str(cycle,'%.3d') '/']) % run S3B L1b monitoring routine for the cycle
end
toc

%%%%%%%%%%%%%%%%%%%% Monitor Thermal data %%%%%%%%%%%%%%%%%%%%
disp(datetime)
tic
disp(['_______ Monitoring CAL1 Thermal for CYCLE ' num2str(cycle) '. Mission ' mission '.']);
switch mission
    case 'S3A'
        S3A_Thermal(cycle, cycle); % run S3A Thermal monitoring routine for the cycle
    case 'S3B'
        S3B_Thermal(cycle, cycle); % run S3B Thermal monitoring routine for the cycle
end

% move back files from cycle folder to general ALLCYCLE folder
movefile([folderL1(1:end-9) num2str(cycle,'%.3d') '\cal_ncread*.mat'], folderL1); % move back cycle L1 files
switch mission
    case 'S3A'
        movefile([folderL1(1:end-9) num2str(cycle,'%.3d') '\S3A_SR_1_CAL_*'], folderL1); % move back cycle L1 folders
    case 'S3B'
        movefile([folderL1(1:end-9) num2str(cycle,'%.3d') '\S3B_SR_1_CAL_*'], folderL1); % move back cycle L1 folders
end
toc

end