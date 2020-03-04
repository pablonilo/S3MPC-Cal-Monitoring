function out = trpf_visinsp(L2folder, lat1, lat2)

% script for the Tracker Performance study.
% visual inspection of wfms.
% examples:
% Baikal Lake
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\data\L2\tracking_performance\S3A_SR_2_LAN____20170303T032407_20170303T041222_20170328T213925_2895_015_061______LN3_O_NT_002.SEN3/', 53, 56)
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\data\L2\tracking_performance\S3B_SR_2_LAN____20190628T032553_20190628T041406_20190723T202618_2893_027_061______LN3_O_NT_003.SEN3/', 52.5, 56)
% Brasilian coast
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\data\L2\tracking_performance\S3A_SR_2_LAN____20170422T005411_20170422T014320_20170423T172551_2949_017_002______LN3_O_ST_002.SEN3/', -23.4, -22.9)
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\data\L2\tracking_performance\S3B_SR_2_LAN____20190624T005553_20190624T014501_20190719T220402_2948_027_002______LN3_O_NT_003.SEN3/', -23.7, -23.2)) 

% naming
name=L2folder(end-99:end-24);
namefig=[name(1:3) '. Cycle ' name(70:72) ', Rel.Orb. ' name(74:76)];

% read L2 enhanced product
L2=readanyNETCDF_V1([L2folder '/enhanced_measurement.nc']);

% change folders for saving
currentworkingfolder=pwd;
cd(L2folder); cd ..;

% constants
c=299792458;
BW=320e6;

% records selection
index = (L2.data.lat_20_c > lat1*1e6 & L2.data.lat_20_c < lat2*1e6);
    
% L2 data reading
out.wfm             = double(L2.data.waveform_20_plrm_ku(:,index));
out.lat             = double(L2.data.lat_20_c(index))*1e-6;
out.lon             = double(L2.data.lon_20_c(index))*1e-6;
out.alt             = double(L2.data.alt_20_c(index))*1e-4;
out.tracker_range   = double(L2.data.tracker_range_20_c(index))*1e-4;
out.surf_class_20_c   = double(L2.data.surf_class_20_c(index))*1e-4;

out.trwfmcor_m      = out.alt-out.tracker_range; 
out.trwfmcor_samp   = out.trwfmcor_m/(1/BW*c/2);

out.newX_ax = 1:128+ceil(abs(min(out.trwfmcor_samp)-max(out.trwfmcor_samp))); % new range bin axis corrected with tracker range
out.wfmd(out.newX_ax,1:length(out.lat))=NaN;
for  i=1:length(out.lat)
    out.wfmd((1:128)-round(out.trwfmcor_samp(i)-max(out.trwfmcor_samp)),i)=out.wfm(1:128,i);
end

kmlwriteline(name,out.lat,out.lon,'width',3)

save(['TRCKperformance_visualinsp_' name],'out');

% lat_apx=round(lat*1e3)/1e3; % latitudes with mdeg resolution
% lat_apx_min=min(lat_apx); lat_apx_max=max(lat_apx); % first and last latitudes
% lat_apx_array=(lat_apx_min:1e-3:lat_apx_max);
% wfmd_fix(1:128+ceil(abs(min(trwfmcor_samp)-max(trwfmcor_samp))),1:length(lat_apx_array))=NaN;
% for i=1:length(lat_apx_array)
%     orig_position=find(round(lat_apx_array(i)*1e3)/1e3 == lat_apx,1);
%     if ~isempty(orig_position)
%         wfmd_fix((1:128)-round(trwfmcor_samp(orig_position)-max(trwfmcor_samp)),i)=wfm(1:128,orig_position);
%     end
% end

%% plots
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',18);
set(0,'defaultTextFontName','Arial');
set(0,'defaultTextFontSize',18);
mida = get(0,'ScreenSize'); mida(3:4)=[1920,1080]; set(0,'defaultFigurePosition',mida);
set(0,'DefaultFigurePaperPositionMode','auto')

clims = [0 5e5];
figure;imagesc(out.lat,[max(out.trwfmcor_m)+128*(1/BW*c/2) min(out.trwfmcor_m)],out.wfmd,clims); % corrected for the tracker jumps
title(['Radargram corrected by Tracker Range. ' namefig]); xlabel('Latitude [deg.]'); ylabel('Elevation [m]'); 
saveas(gcf,['WFMS_corrbyTR__' name],'jpg'); saveas(gcf,['WFMS_corrbyTR__' name],'fig');
close;

% 
% figure; plot(lat,trwfmcor_samp,'.','linewidth', 3) % tracker jumps
% title('Altitude - Tracker Range in range bin units'); xlabel('Latitude [degrees]'); ylabel('Altitude - Tracker Range [range bins]'); 
% saveas(gcf,'ALTminusTRACKRANGE','jpg'); saveas(gcf,'ALTminusTRACKRANGE','fig');
% close;
% 
% clims = [0 1e5];
% figure;imagesc(lat,1:128,wfm,clims); % original wfm with tracker jumps
% title('Waveforms without Tracker Range correction'); xlabel('Latitude [degrees]'); ylabel('range bin'); 
% saveas(gcf,'WFMS_notcorrbyTR','jpg'); saveas(gcf,'WFMS_notcorrbyTR','fig');
% close;

cd(currentworkingfolder);
end