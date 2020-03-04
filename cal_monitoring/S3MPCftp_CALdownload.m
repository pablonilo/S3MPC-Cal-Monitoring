function USOdirout=S3MPCftp_CALdownload(mission, cycle)
% download S3 CAL data from S3MPC ftp site
% calling example (CUBA): S3MPCftp_CALdownload('S3A', 44)

%% Initializing folders locations
cyc=num2str(cycle);
ftp_folder_L1=['/' mission '/SRAL/SR_1/SR_1_CAL/Cycle0' cyc];
ftp_folder_L0=['/' mission '/SRAL/SR_0/SR_0_CAL/Cycle0' cyc];

if strcmp(mission,'S3A')
    localfolderL1='C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_1_CAL\CYCLE_ALL\';
    localfolderL0='C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_0_CAL\';
    localfolderUSO='C:\Users\Pablo\Desktop\S3MPC\data\CAL\USO\';
    ftp_folder_USO='/ADFs/all_current_ADFs/20160223';
elseif strcmp(mission,'S3B')
    localfolderL1='C:\Users\Pablo\Desktop\S3MPC\S3B\data\CAL\SR_1_CAL\CYCLE_ALL\';
    localfolderL0='C:\Users\Pablo\Desktop\S3MPC\S3B\data\CAL\SR_0_CAL\';
    localfolderUSO='C:\Users\Pablo\Desktop\S3MPC\S3B\data\CAL\USO\';
    ftp_folder_USO='/ADFs/all_current_ADFs/20180501';
end

%% setting the ftp transfer
disp('Opening connection to ACRI S3MPC ftp site');
ts = ftp('ftp.acri-cwa.fr','ftp_s3mpc-stm','hyu258::'); % Connect to the FTP server host and create the FTP object
binary(ts); % set to binary ftp transfer type

%% reading folders from ftp
disp('Creating folders list from ftp site');
dir_L1 = dir(ts,ftp_folder_L1); % list the existent L1 folders
dir_L0 = dir(ts,ftp_folder_L0); % list the existent L0 folders

% Avoid downloading L0 Autocal prods (too big)
disp('Avoiding L0 AutoCal folders (too big files)');
list_Autocal=[22 157 278];
for i=length(dir_L0):-1:1
    if ~isempty(find(str2double(dir_L0(i).name(74:76))==list_Autocal,1))
        dir_L0(i)=[]; % delete AutoCal product folder from list
    end
end

%% reading L1 and L0 CAL files from ftp
disp(['Downloading ' num2str(length(dir_L1)) ' products']);
for i=1:length(dir_L1)
 disp(['Downloading L1 CAL product nº ' num2str(i) '. Completed: ' num2str(100*i/length(dir_L1)) '% . foldername: ' dir_L1(i).name]);
    mget(ts,[ftp_folder_L1 '/' dir_L1(i).name],localfolderL1); % download L1
end
for i=1:length(dir_L0)
 disp(['Downloading L0 CAL product nº ' num2str(i) '. Completed: ' num2str(100*i/length(dir_L1)) '% . foldername: ' dir_L0(i).name]);
    mget(ts,[ftp_folder_L0 '/' dir_L0(i).name],localfolderL0); % download L0
end

% move files to root & delete ftp-like folder
disp('Relocating L1 and L0 CAL files in laptop folder');
movefile([localfolderL1 mission '\SRAL\SR_1\SR_1_CAL\Cycle0' cyc '\*'],localfolderL1);
rmdir([localfolderL1 mission],'s');
movefile([localfolderL0 mission '\SRAL\SR_0\SR_0_CAL\Cycle0' cyc '\*'],localfolderL0);
rmdir([localfolderL0 mission],'s');

%% USO file downloading
dir_USO=dir(ts,[ftp_folder_USO '/*USO*']);
for i=1:length(dir_USO)
    if strcmp(dir_L1(end).name(end-66:end-59), dir_USO(i).name(end-66:end-59))
        disp(['Downloading USO file: ' dir_USO(i).name]);
        mget(ts,dir_USO(i).name,localfolderUSO); % download USO coincident with last cycle date from L1 CAL files
        break;
    end
end
dir_USOzipfile=dir(ts,[dir_USO(i).name '/' mission '*']);
disp('Unzipping the USO binary file into the proper folder and removing the ftp-like folder');
untar([localfolderUSO dir_USOzipfile.name(2:end)], [localfolderUSO dir_USOzipfile.name(end-97:end-35)]);
rmdir([localfolderUSO 'ADFs'],'s'); % remove ftp-like folder 

USOdirout=[localfolderUSO dir_USOzipfile.name(end-97:end-35)]; % for Automatic Monitoring 

disp('Closing connection to ACRI S3MPC ftp site');
close(ts); % Close the connection to the S3MPC ftp server
end