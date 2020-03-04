function trpf_visinsp_RSR(folder)
tic
% script for the Tracker Performance study.
% visual inspection of wfms.
% examples:
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\Routine\Docs\TRCK_performance\for RSR09\');
%
% old......
% Baikal Lake
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\data\L2\tracking_performance\S3A_SR_2_LAN____20170303T032407_20170303T041222_20170328T213925_2895_015_061______LN3_O_NT_002.SEN3/', 53, 56)
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\data\L2\tracking_performance\S3B_SR_2_LAN____20190628T032553_20190628T041406_20190723T202618_2893_027_061______LN3_O_NT_003.SEN3/', 52.5, 56)
% Brasilian coast
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\data\L2\tracking_performance\S3A_SR_2_LAN____20170422T005411_20170422T014320_20170423T172551_2949_017_002______LN3_O_ST_002.SEN3/', -23.4, -22.9)
% out = trpf_visinsp('C:\Users\Pablo\Desktop\S3MPC\data\L2\tracking_performance\S3B_SR_2_LAN____20190624T005553_20190624T014501_20190719T220402_2948_027_002______LN3_O_NT_003.SEN3/', -23.7, -23.2))

listdir_02A=dir([folder 'S3A*_002_*']);
listdir_02B=dir([folder 'S3B*_002_*']);
listdir_61A=dir([folder 'S3A*_061_*']);
listdir_61B=dir([folder 'S3B*_061_*']);

for i=1:length(listdir_02A)
    out = trpf_visinsp([folder listdir_02A(i).name,'/'], -23.4, -22.9);
    clear out;
end

for i=1:length(listdir_02B)
    out = trpf_visinsp([folder listdir_02B(i).name,'/'], -23.7, -23.2);
    clear out;
end

for i=1:length(listdir_61A)
    out = trpf_visinsp([folder listdir_61A(i).name,'/'], 53, 56);
    clear out;
end

for i=1:length(listdir_61B)
    out = trpf_visinsp([folder listdir_61B(i).name,'/'], 52.5, 56);
    clear out;
end

toc
end