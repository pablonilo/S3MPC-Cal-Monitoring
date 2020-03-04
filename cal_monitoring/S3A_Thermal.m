function S3A_ISP_TH = S3A_Thermal(in_cycle_ini,in_cycle_end)
%Thermal_gather1 Update of the S3A Thermal Repository
% inputs: first and last cycles.  in_cycle_ini: 0 if whole mission.
% output: The Thermal Repository is updated


%% select & read not readed ISP files
folderL0='C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_0_CAL\'; % S3A ISP folder
L0isp=dir([folderL0 'S3A_SR_0_CAL____*']); % get all L0 (folders) & matlab (files) ISP names
L0matpre=dir([folderL0 'S3A_SR_0_CAL____*.mat']); % get only matlab ISP files names

if size(L0matpre,1)~=0
    % list of mat L0 files
    for i_file=1:length(L0matpre)
%        matchL0.namemat(i_file,1:length(L0matpre(i_file,1).name))=L0matpre(i_file,1).name;
        matchL0.datemat(i_file,1)=datenum([L0matpre(i_file,1).name(17:24) L0matpre(i_file,1).name(26:31)],'yyyymmddHHMMSS');  % sensing start
    end
    % list of isp L0 folders
    count_ispfile=0;
    for i_file=1:length(L0isp)
        if L0isp(i_file,1).isdir
            count_ispfile=count_ispfile+1;
            matchL0.nameispdir(count_ispfile,1:length(L0isp(i_file,1).name))=L0isp(i_file,1).name;
            matchL0.dateispdir(count_ispfile,1)=datenum([L0isp(i_file,1).name(17:24) L0isp(i_file,1).name(26:31)],'yyyymmddHHMMSS');  % sensing start
        end
    end
    
    % initialize ISPs list to read
    ISPread(1,:)='no ISPs left to read';
    % fill the ISPs list to read
    count=0;
    for i_file=1:length(matchL0.nameispdir)
        matchL0.flag(i_file)=0;
        if ~isempty(find(matchL0.dateispdir(i_file)==matchL0.datemat,1)) % match L0 matlab file and ISP file names?
            matchL0.flag(i_file)=1; % match found --> do not read
        else
            count=count+1;
            ISPread(count,1:length(matchL0.nameispdir(i_file,:)))=matchL0.nameispdir(i_file,:);
        end
    end
    
else % no matlab L0 files, read all ISPs
    count=size(L0isp,1);
    for i=1:length(L0isp)
        ISPread(i,:) = L0isp(i).name;
    end
end

% Read the ISPs that left to read calling the ISP reading routine
if ~strcmp(ISPread(1,:),'no ISPs left to read')
    disp(['reading ' num2str(count) ' S3 SRAL ISP products from parent folder: ' folderL0]);
    for i_folder=1:length(ISPread(:,1))
        disp(['... reading file number ' num2str(i_folder) ' of ' num2str(length(ISPread(:,1)))]);
        S3_ISP_data_read([folderL0 ISPread(i_folder,:) '/ISPData.dat']); % read the ISP and save a matlab ISP file
        clear ans % delete the returned read structure
    end
else
    disp(['No S3 SRAL ISP products left to read from folder: ' folderL0]);
end

%% fill the repository

d=dir('C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_0_CAL\S3A_SR_0_CAL____20*.mat'); % list of matlab ISP files
load('C:\Users\Pablo\Desktop\S3MPC\Routine\code\S3A_ISP_TH.mat'); % open mission thermal data matlab repository

count_ini=size(S3A_ISP_TH.file.cycle,2); % initial fixed size of the thermal data repository

count=count_ini; % initialize record number of thermal data repository
for i=1:length(d)
    sens_start=datenum([d(i).name(17:24) d(i).name(26:31)],'yyyymmddHHMMSS');
    if isempty(find(S3A_ISP_TH.file.timestart==sens_start,1)) % if the date does not coincide with one in the historic
        load([d(i).folder '/' d(i).name]); % open ith matlab ISP file
        if isfield(L0,'cal1lrmiq') & isfield(L0,'cal1sar') & isfield(L0,'cal2')  % if there are the 3 CAL modes
            count=count+1;
            disp(['Adding new line to thermal repository from file: ' d(i).name(1:94)]);
            % cal1 lrmiq mode
            S3A_ISP_TH.cal1lrm_iq.mean_TH_KU_SSPA_THR1(count)       = mean(L0.cal1lrmiq.therm.TH_KU_SSPA_THR1(L0.cal1lrmiq.therm.TH_KU_SSPA_THR1~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_C_SSPA_THR2(count)        = mean(L0.cal1lrmiq.therm.TH_C_SSPA_THR2(L0.cal1lrmiq.therm.TH_C_SSPA_THR2~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_KU_RX_THR3(count)         = mean(L0.cal1lrmiq.therm.TH_KU_RX_THR3(L0.cal1lrmiq.therm.TH_KU_RX_THR3~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_C_RX_THR4(count)          = mean(L0.cal1lrmiq.therm.TH_C_RX_THR4(L0.cal1lrmiq.therm.TH_C_RX_THR4~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_FI_CNG_THR5(count)        = mean(L0.cal1lrmiq.therm.TH_FI_CNG_THR5(L0.cal1lrmiq.therm.TH_FI_CNG_THR5~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_X16_THR6 (count)          = mean(L0.cal1lrmiq.therm.TH_X16_THR6(L0.cal1lrmiq.therm.TH_X16_THR6~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_DCDC_RFU_THR7(count)      = mean(L0.cal1lrmiq.therm.TH_DCDC_RFU_THR7(L0.cal1lrmiq.therm.TH_DCDC_RFU_THR7~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_ref_can1(count)              = mean(L0.cal1lrmiq.therm.ref_can1(L0.cal1lrmiq.therm.ref_can1~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_PLL_THR8(count)           = mean(L0.cal1lrmiq.therm.TH_PLL_THR8(L0.cal1lrmiq.therm.TH_PLL_THR8~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_CHIRP_THR9(count)         = mean(L0.cal1lrmiq.therm.TH_CHIRP_THR9(L0.cal1lrmiq.therm.TH_CHIRP_THR9~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_FI_THR10(count)           = mean(L0.cal1lrmiq.therm.TH_FI_THR10(L0.cal1lrmiq.therm.TH_FI_THR10~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_ASIC_DFFT_THR11(count)    = mean(L0.cal1lrmiq.therm.TH_ASIC_DFFT_THR11(L0.cal1lrmiq.therm.TH_ASIC_DFFT_THR11~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_TH_DCDC_THR12(count)         = mean(L0.cal1lrmiq.therm.TH_DCDC_THR12(L0.cal1lrmiq.therm.TH_DCDC_THR12~=0));
            S3A_ISP_TH.cal1lrm_iq.mean_ref_can2(count)              = mean(L0.cal1lrmiq.therm.ref_can2(L0.cal1lrmiq.therm.ref_can2~=0));
            % cal1 sar mode
            S3A_ISP_TH.cal1sar.mean_TH_KU_SSPA_THR1(count)      = mean(L0.cal1sar.therm.TH_KU_SSPA_THR1(L0.cal1sar.therm.TH_KU_SSPA_THR1~=0));
            S3A_ISP_TH.cal1sar.mean_TH_C_SSPA_THR2(count)       = mean(L0.cal1sar.therm.TH_C_SSPA_THR2(L0.cal1sar.therm.TH_C_SSPA_THR2~=0));
            S3A_ISP_TH.cal1sar.mean_TH_KU_RX_THR3(count)        = mean(L0.cal1sar.therm.TH_KU_RX_THR3(L0.cal1sar.therm.TH_KU_RX_THR3~=0));
            S3A_ISP_TH.cal1sar.mean_TH_C_RX_THR4(count)         = mean(L0.cal1sar.therm.TH_C_RX_THR4(L0.cal1sar.therm.TH_C_RX_THR4~=0));
            S3A_ISP_TH.cal1sar.mean_TH_FI_CNG_THR5(count)       = mean(L0.cal1sar.therm.TH_FI_CNG_THR5(L0.cal1sar.therm.TH_FI_CNG_THR5~=0));
            S3A_ISP_TH.cal1sar.mean_TH_X16_THR6(count)          = mean(L0.cal1sar.therm.TH_X16_THR6(L0.cal1sar.therm.TH_X16_THR6~=0));
            S3A_ISP_TH.cal1sar.mean_TH_DCDC_RFU_THR7(count)     = mean(L0.cal1sar.therm.TH_DCDC_RFU_THR7(L0.cal1sar.therm.TH_DCDC_RFU_THR7~=0));
            S3A_ISP_TH.cal1sar.mean_ref_can1(count)             = mean(L0.cal1sar.therm.ref_can1(L0.cal1sar.therm.ref_can1~=0));
            S3A_ISP_TH.cal1sar.mean_TH_PLL_THR8(count)          = mean(L0.cal1sar.therm.TH_PLL_THR8(L0.cal1sar.therm.TH_PLL_THR8~=0));
            S3A_ISP_TH.cal1sar.mean_TH_CHIRP_THR9(count)        = mean(L0.cal1sar.therm.TH_CHIRP_THR9(L0.cal1sar.therm.TH_CHIRP_THR9~=0));
            S3A_ISP_TH.cal1sar.mean_TH_FI_THR10(count)          = mean(L0.cal1sar.therm.TH_FI_THR10(L0.cal1sar.therm.TH_FI_THR10~=0));
            S3A_ISP_TH.cal1sar.mean_TH_ASIC_DFFT_THR11(count)   = mean(L0.cal1sar.therm.TH_ASIC_DFFT_THR11(L0.cal1sar.therm.TH_ASIC_DFFT_THR11~=0));
            S3A_ISP_TH.cal1sar.mean_TH_DCDC_THR12(count)        = mean(L0.cal1sar.therm.TH_DCDC_THR12(L0.cal1sar.therm.TH_DCDC_THR12~=0));
            S3A_ISP_TH.cal1sar.mean_ref_can2(count)             = mean(L0.cal1sar.therm.ref_can2(L0.cal1sar.therm.ref_can2~=0));
            % cal2 mode
            S3A_ISP_TH.cal2.mean_TH_KU_SSPA_THR1(count)     = mean(L0.cal2.therm.TH_KU_SSPA_THR1(L0.cal2.therm.TH_KU_SSPA_THR1~=0));
            S3A_ISP_TH.cal2.mean_TH_C_SSPA_THR2(count)      = mean(L0.cal2.therm.TH_C_SSPA_THR2(L0.cal2.therm.TH_C_SSPA_THR2~=0));
            S3A_ISP_TH.cal2.mean_TH_KU_RX_THR3(count)       = mean(L0.cal2.therm.TH_KU_RX_THR3(L0.cal2.therm.TH_KU_RX_THR3~=0));
            S3A_ISP_TH.cal2.mean_TH_C_RX_THR4(count)        = mean(L0.cal2.therm.TH_C_RX_THR4(L0.cal2.therm.TH_C_RX_THR4~=0));
            S3A_ISP_TH.cal2.mean_TH_FI_CNG_THR5(count)      = mean(L0.cal2.therm.TH_FI_CNG_THR5(L0.cal2.therm.TH_FI_CNG_THR5~=0));
            S3A_ISP_TH.cal2.mean_TH_X16_THR6(count)         = mean(L0.cal2.therm.TH_X16_THR6(L0.cal2.therm.TH_X16_THR6~=0));
            S3A_ISP_TH.cal2.mean_TH_DCDC_RFU_THR7(count)    = mean(L0.cal2.therm.TH_DCDC_RFU_THR7(L0.cal2.therm.TH_DCDC_RFU_THR7~=0));
            S3A_ISP_TH.cal2.mean_ref_can1(count)            = mean(L0.cal2.therm.ref_can1(L0.cal2.therm.ref_can1~=0));
            S3A_ISP_TH.cal2.mean_TH_PLL_THR8(count)         = mean(L0.cal2.therm.TH_PLL_THR8(L0.cal2.therm.TH_PLL_THR8~=0));
            S3A_ISP_TH.cal2.mean_TH_CHIRP_THR9(count)       = mean(L0.cal2.therm.TH_CHIRP_THR9(L0.cal2.therm.TH_CHIRP_THR9~=0));
            S3A_ISP_TH.cal2.mean_TH_FI_THR10(count)         = mean(L0.cal2.therm.TH_FI_THR10(L0.cal2.therm.TH_FI_THR10~=0));
            S3A_ISP_TH.cal2.mean_TH_ASIC_DFFT_THR11(count)  = mean(L0.cal2.therm.TH_ASIC_DFFT_THR11(L0.cal2.therm.TH_ASIC_DFFT_THR11~=0));
            S3A_ISP_TH.cal2.mean_TH_DCDC_THR12(count)       = mean(L0.cal2.therm.TH_DCDC_THR12(L0.cal2.therm.TH_DCDC_THR12~=0));
            S3A_ISP_TH.cal2.mean_ref_can2(count)            = mean(L0.cal2.therm.ref_can2(L0.cal2.therm.ref_can2~=0));
            % gather names
            S3A_ISP_TH.file.name(count,:)=d(i).name(1:94);
            S3A_ISP_TH.file.cycle(count)=str2double(d(i).name(70:72));
            S3A_ISP_TH.file.relorb(count)=str2double(d(i).name(74:76));
            S3A_ISP_TH.file.timestart(count)=datenum([d(i).name(17:24) d(i).name(26:31)],'yyyymmddHHMMSS');
            S3A_ISP_TH.file.timeend(count)=datenum([d(i).name(33:40) d(i).name(42:47)],'yyyymmddHHMMSS');
            S3A_ISP_TH.file.timeproc(count)=datenum([d(i).name(49:56) d(i).name(58:63)],'yyyymmddHHMMSS');
        end
    end
end

% save the updated repository if updated
if count_ini<count
    disp('Saving the Thermal repository as C:\Users\Pablo\Desktop\S3MPC\Routine\code\S3A_ISP_TH.mat');
    save('C:\Users\Pablo\Desktop\S3MPC\Routine\code\S3A_ISP_TH.mat','S3A_ISP_TH');
end

%% PLOTTING
disp(['Plotting S3A Thermal series from cycle ' num2str(in_cycle_ini) ' to cycle ' num2str(in_cycle_end)])
% initialise figures format
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontName','Arial');
set(0,'defaultTextFontSize',16);
mida = get(0,'ScreenSize'); mida(3:4)=[1920,1080]; set(0,'defaultFigurePosition',mida);
set(0,'DefaultFigurePaperPositionMode','auto')

folderplot='C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_1_CAL\CYCLE_ALL\';

if in_cycle_ini==in_cycle_end
    folderplot=['C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_1_CAL\C0' num2str(in_cycle_end) '\'];
else
    folderplot='C:\Users\Pablo\Desktop\S3MPC\data\CAL\SR_1_CAL\CYCLE_ALL\';
end

timeplot=S3A_ISP_TH.file.timeend; % General time axis for plots

% filter the data for the cycles of interest
filt=find(S3A_ISP_TH.file.cycle>=in_cycle_ini & S3A_ISP_TH.file.cycle<=in_cycle_end);

% CAL1 LRM IQ
namef=['S3MPC_S3A_CAL1LRMIQ_Therm_'  datestr(timeplot(filt(1)),1) '_' datestr(timeplot(filt(end)),1)];
figure; hold all;
title(['S3A SRAL CAL1 LRM IQ THERM. Selected thermistors time series. From ' datestr(timeplot(filt(1)),1) ' to ' datestr(timeplot(filt(end)),1) '.']);
plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_KU_SSPA_THR1(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_KU_SSPA_THR1(filt(1)))/10,'o','DisplayName', ['TH-KU-SSPA-THR1, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_KU_SSPA_THR1(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_C_SSPA_THR2(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_C_SSPA_THR2(filt(1)))/10,'o','DisplayName', ['TH-C-SSPA-THR2, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_C_SSPA_THR2(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_ASIC_DFFT_THR11(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_ASIC_DFFT_THR11(filt(1)))/10,'o','DisplayName', ['ASIC-DFFT-THR11, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_ASIC_DFFT_THR11(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_KU_RX_THR3(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_KU_RX_THR3(filt(1)))/10,'o','DisplayName', ['TH-KU-RX-THR3, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_KU_RX_THR3(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_C_RX_THR4(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_C_RX_THR4(filt(1)))/10,'o','DisplayName', ['TH-C-RX-THR4, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_C_RX_THR4(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_FI_CNG_THR5(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_FI_CNG_THR5(filt(1)))/10,'o','DisplayName', ['TH-FI-CNG-THR5, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_FI_CNG_THR5(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_X16_THR6(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_X16_THR6(filt(1)))/10,'o','DisplayName', ['TH-X16-THR6, from ' num2str((filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_DCDC_RFU_THR7(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_DCDC_RFU_THR7(filt(1)))/10,'o','DisplayName', ['TH-DCDC-RFU-THR7, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_DCDC_RFU_THR7(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_ref_can1(filt)-S3A_ISP_TH.cal1lrm_iq.mean_ref_can1(filt(1)))/10,'o','DisplayName', ['mean-ref-can1, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_ref_can1(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_PLL_THR8(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_PLL_THR8(filt(1)))/10,'o','DisplayName', ['TH-PLL-THR8, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_PLL_THR8(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_CHIRP_THR9(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_CHIRP_THR9(filt(1)))/10,'o','DisplayName', ['TH-CHIRP-THR9, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_CHIRP_THR9(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_FI_THR10(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_FI_THR10(filt(1)))/10,'o','DisplayName', ['TH-FI-THR10, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_FI_THR10(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_ASIC_DFFT_THR11(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_ASIC_DFFT_THR11(filt(1)))/10,'o','DisplayName', ['ASIC-DFFT-THR11, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_ASIC_DFFT_THR11(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_TH_DCDC_THR12(filt)-S3A_ISP_TH.cal1lrm_iq.mean_TH_DCDC_THR12(filt(1)))/10,'o','DisplayName', ['TH-DCDC-THR12, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_TH_DCDC_THR12(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal1lrm_iq.mean_ref_can2(filt)-S3A_ISP_TH.cal1lrm_iq.mean_ref_can2(filt(1)))/10,'o','DisplayName', ['mean-ref-can2, from ' num2str(S3A_ISP_TH.cal1lrm_iq.mean_ref_can2(filt(1))/10) 'ºC'], 'linewidth', 2)
legend('show','location', 'northeast');
xlabel('Date'); ylabel('Thermistors Values [ºC]');
twide=timeplot(filt(end))-timeplot(filt(1)); 
if  twide < 90
    xticks([twide/6+timeplot(filt(1)) twide/6*2+timeplot(filt(1)) twide/6*3+timeplot(filt(1)) twide/6*4+timeplot(filt(1)) twide/6*5+timeplot(filt(1))]); % 5 ticks in time axis
    datetick('x',20,'keepticks') % mmmyydd for short time series
    xlim([timeplot(filt(1))-1 timeplot(filt(end))+1]);
else
    datetick('x',12,'keeplimits') % mmmyy for long time series
    xlim([timeplot(filt(1))-5 timeplot(filt(end))+5]);
end
saveas(gcf,[folderplot namef],'jpg');
saveas(gcf,[folderplot namef],'fig');
close(gcf);

% CAL1 SAR NORM
namef=['S3MPC_S3A_CAL1SAR_Therm_'  datestr(timeplot(filt(1)),1) '_' datestr(timeplot(filt(end)),1)];
figure; hold all;
title(['S3A SRAL CAL1 SAR THERM. Selected thermistors time series. From ' datestr(timeplot(filt(1)),1) ' to ' datestr(timeplot(filt(end)),1) '.']);
plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_KU_SSPA_THR1(filt)-S3A_ISP_TH.cal1sar.mean_TH_KU_SSPA_THR1(filt(1)))/10,'o','DisplayName', ['TH-KU-SSPA-THR1, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_KU_SSPA_THR1(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_C_SSPA_THR2(filt)-S3A_ISP_TH.cal1sar.mean_TH_C_SSPA_THR2(filt(1)))/10,'o','DisplayName', ['TH-C-SSPA-THR2, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_C_SSPA_THR2(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_ASIC_DFFT_THR11(filt)-S3A_ISP_TH.cal1sar.mean_TH_ASIC_DFFT_THR11(filt(1)))/10,'o','DisplayName', ['ASIC-DFFT-THR11, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_ASIC_DFFT_THR11(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_KU_RX_THR3(filt)-S3A_ISP_TH.cal1sar.mean_TH_KU_RX_THR3(filt(1)))/10,'o','DisplayName', ['TH-KU-RX-THR3, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_KU_RX_THR3(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_C_RX_THR4(filt)-S3A_ISP_TH.cal1sar.mean_TH_C_RX_THR4(filt(1)))/10,'o','DisplayName', ['TH-C-RX-THR4, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_C_RX_THR4(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_FI_CNG_THR5(filt)-S3A_ISP_TH.cal1sar.mean_TH_FI_CNG_THR5(filt(1)))/10,'o','DisplayName', ['TH-FI-CNG-THR5, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_FI_CNG_THR5(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_X16_THR6(filt)-S3A_ISP_TH.cal1sar.mean_TH_X16_THR6(filt(1)))/10,'o','DisplayName', ['TH-X16-THR6, from ' num2str((filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_DCDC_RFU_THR7(filt)-S3A_ISP_TH.cal1sar.mean_TH_DCDC_RFU_THR7(filt(1)))/10,'o','DisplayName', ['TH-DCDC-RFU-THR7, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_DCDC_RFU_THR7(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_ref_can1(filt)-S3A_ISP_TH.cal1sar.mean_ref_can1(filt(1)))/10,'o','DisplayName', ['mean-ref-can1, from ' num2str(S3A_ISP_TH.cal1sar.mean_ref_can1(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_PLL_THR8(filt)-S3A_ISP_TH.cal1sar.mean_TH_PLL_THR8(filt(1)))/10,'o','DisplayName', ['TH-PLL-THR8, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_PLL_THR8(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_CHIRP_THR9(filt)-S3A_ISP_TH.cal1sar.mean_TH_CHIRP_THR9(filt(1)))/10,'o','DisplayName', ['TH-CHIRP-THR9, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_CHIRP_THR9(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_FI_THR10(filt)-S3A_ISP_TH.cal1sar.mean_TH_FI_THR10(filt(1)))/10,'o','DisplayName', ['TH-FI-THR10, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_FI_THR10(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_ASIC_DFFT_THR11(filt)-S3A_ISP_TH.cal1sar.mean_TH_ASIC_DFFT_THR11(filt(1)))/10,'o','DisplayName', ['ASIC-DFFT-THR11, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_ASIC_DFFT_THR11(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_TH_DCDC_THR12(filt)-S3A_ISP_TH.cal1sar.mean_TH_DCDC_THR12(filt(1)))/10,'o','DisplayName', ['TH-DCDC-THR12, from ' num2str(S3A_ISP_TH.cal1sar.mean_TH_DCDC_THR12(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal1sar.mean_ref_can2(filt)-S3A_ISP_TH.cal1sar.mean_ref_can2(filt(1)))/10,'o','DisplayName', ['mean-ref-can2, from ' num2str(S3A_ISP_TH.cal1sar.mean_ref_can2(filt(1))/10) 'ºC'], 'linewidth', 2)
legend('show','location', 'northeast');
xlabel('Date'); ylabel('Thermistors Values [ºC]');
twide=timeplot(filt(end))-timeplot(filt(1)); 
if  twide < 90
    xticks([twide/6+timeplot(filt(1)) twide/6*2+timeplot(filt(1)) twide/6*3+timeplot(filt(1)) twide/6*4+timeplot(filt(1)) twide/6*5+timeplot(filt(1))]); % 5 ticks in time axis
    datetick('x',20,'keepticks') % mmmyydd for short time series
    xlim([timeplot(filt(1))-1 timeplot(filt(end))+1]);
else
    datetick('x',12,'keeplimits') % mmmyy for long time series
    xlim([timeplot(filt(1))-5 timeplot(filt(end))+5]);
end
saveas(gcf,[folderplot namef],'jpg');
saveas(gcf,[folderplot namef],'fig');
close(gcf);

% CAL2
namef=['S3MPC_S3A_CAL2_Therm_'  datestr(timeplot(filt(1)),1) '_' datestr(timeplot(filt(end)),1)];
figure; hold all;
title(['S3A SRAL CAL2 THERM. Selected thermistors time series. From ' datestr(timeplot(filt(1)),1) ' to ' datestr(timeplot(filt(end)),1) '.']);
plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_KU_SSPA_THR1(filt)-S3A_ISP_TH.cal2.mean_TH_KU_SSPA_THR1(filt(1)))/10,'o','DisplayName', ['TH-KU-SSPA-THR1, from ' num2str(S3A_ISP_TH.cal2.mean_TH_KU_SSPA_THR1(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_C_SSPA_THR2(filt)-S3A_ISP_TH.cal2.mean_TH_C_SSPA_THR2(filt(1)))/10,'o','DisplayName', ['TH-C-SSPA-THR2, from ' num2str(S3A_ISP_TH.cal2.mean_TH_C_SSPA_THR2(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_ASIC_DFFT_THR11(filt)-S3A_ISP_TH.cal2.mean_TH_ASIC_DFFT_THR11(filt(1)))/10,'o','DisplayName', ['ASIC-DFFT-THR11, from ' num2str(S3A_ISP_TH.cal2.mean_TH_ASIC_DFFT_THR11(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_KU_RX_THR3(filt)-S3A_ISP_TH.cal2.mean_TH_KU_RX_THR3(filt(1)))/10,'o','DisplayName', ['TH-KU-RX-THR3, from ' num2str(S3A_ISP_TH.cal2.mean_TH_KU_RX_THR3(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_C_RX_THR4(filt)-S3A_ISP_TH.cal2.mean_TH_C_RX_THR4(filt(1)))/10,'o','DisplayName', ['TH-C-RX-THR4, from ' num2str(S3A_ISP_TH.cal2.mean_TH_C_RX_THR4(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_FI_CNG_THR5(filt)-S3A_ISP_TH.cal2.mean_TH_FI_CNG_THR5(filt(1)))/10,'o','DisplayName', ['TH-FI-CNG-THR5, from ' num2str(S3A_ISP_TH.cal2.mean_TH_FI_CNG_THR5(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_X16_THR6(filt)-S3A_ISP_TH.cal2.mean_TH_X16_THR6(filt(1)))/10,'o','DisplayName', ['TH-X16-THR6, from ' num2str((filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_DCDC_RFU_THR7(filt)-S3A_ISP_TH.cal2.mean_TH_DCDC_RFU_THR7(filt(1)))/10,'o','DisplayName', ['TH-DCDC-RFU-THR7, from ' num2str(S3A_ISP_TH.cal2.mean_TH_DCDC_RFU_THR7(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_ref_can1(filt)-S3A_ISP_TH.cal2.mean_ref_can1(filt(1)))/10,'o','DisplayName', ['mean-ref-can1, from ' num2str(S3A_ISP_TH.cal2.mean_ref_can1(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_PLL_THR8(filt)-S3A_ISP_TH.cal2.mean_TH_PLL_THR8(filt(1)))/10,'o','DisplayName', ['TH-PLL-THR8, from ' num2str(S3A_ISP_TH.cal2.mean_TH_PLL_THR8(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_CHIRP_THR9(filt)-S3A_ISP_TH.cal2.mean_TH_CHIRP_THR9(filt(1)))/10,'o','DisplayName', ['TH-CHIRP-THR9, from ' num2str(S3A_ISP_TH.cal2.mean_TH_CHIRP_THR9(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_FI_THR10(filt)-S3A_ISP_TH.cal2.mean_TH_FI_THR10(filt(1)))/10,'o','DisplayName', ['TH-FI-THR10, from ' num2str(S3A_ISP_TH.cal2.mean_TH_FI_THR10(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_ASIC_DFFT_THR11(filt)-S3A_ISP_TH.cal2.mean_TH_ASIC_DFFT_THR11(filt(1)))/10,'o','DisplayName', ['ASIC-DFFT-THR11, from ' num2str(S3A_ISP_TH.cal2.mean_TH_ASIC_DFFT_THR11(filt(1))/10) 'ºC'], 'linewidth', 2)
%     plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_TH_DCDC_THR12(filt)-S3A_ISP_TH.cal2.mean_TH_DCDC_THR12(filt(1)))/10,'o','DisplayName', ['TH-DCDC-THR12, from ' num2str(S3A_ISP_TH.cal2.mean_TH_DCDC_THR12(filt(1))/10) 'ºC'], 'linewidth', 2)
plot(timeplot(filt),(S3A_ISP_TH.cal2.mean_ref_can2(filt)-S3A_ISP_TH.cal2.mean_ref_can2(filt(1)))/10,'o','DisplayName', ['mean-ref-can2, from ' num2str(S3A_ISP_TH.cal2.mean_ref_can2(filt(1))/10) 'ºC'], 'linewidth', 2)
legend('show','location', 'northeast');
xlabel('Date'); ylabel('Thermistors Values [ºC]');
twide=timeplot(filt(end))-timeplot(filt(1)); 
if  twide < 90
    xticks([twide/6+timeplot(filt(1)) twide/6*2+timeplot(filt(1)) twide/6*3+timeplot(filt(1)) twide/6*4+timeplot(filt(1)) twide/6*5+timeplot(filt(1))]); % 5 ticks in time axis
    datetick('x',20,'keepticks') % mmmyydd for short time series
    xlim([timeplot(filt(1))-1 timeplot(filt(end))+1]);
else
    datetick('x',12,'keeplimits') % mmmyy for long time series
    xlim([timeplot(filt(1))-5 timeplot(filt(end))+5]);
end
saveas(gcf,[folderplot namef],'jpg');
saveas(gcf,[folderplot namef],'fig');
close(gcf);

end