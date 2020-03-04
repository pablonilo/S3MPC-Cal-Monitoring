function out=S3AB_USO_Monitor(filein,mission)
% USOread
% Read a DBL USO aux file and produce the frequency and range impact figures for the cyclic reports
% Call example: S3AB_USO_Monitor('C:\Users\Pablo\Desktop\S3MPC\data\CAL\USO\S3A_SR_1_USO_AX_20160223T195017_20181109T014037_20181109T084742___________________CNE_O_AL_001.SEN3\S3A_SR_1_USO_AX_20160223T195017_20181109T014037_20181109T075537___________________CNE_O_______.DBL','S3A');

disp(['Monitoring ' mission ' USO drift from file ' filein]);
%% reading the USO ASCII file

fid=fopen(filein); %open file
filespecs=dir(filein);
[fparts.pathstr,fparts.name,fparts.ext] = fileparts(filein);
% headersize=2905;
headersize=2969;

fseek(fid,headersize,'bof'); %jump whole header
for i=1:(filespecs.bytes-headersize)/66
    out.day(i)=fscanf(fid,'%5i',1);
    out.sec(i)=fscanf(fid,'%12f',1);
    fseek(fid,4,'cof');
    out.USOmeasured(i)=fscanf(fid,'%11f',1);
    fseek(fid,4,'cof');
    out.USOmodelled(i)=fscanf(fid,'%11f',1);
    fseek(fid,4,'cof');
    out.USOresidual(i)=fscanf(fid,'%11f',1);
    fseek(fid,1,'cof');
    out.val_flag(i)=fscanf(fid,'%1i',1);
end
fclose(fid);

%% dates comparison between filename and data 
% uso date from data
out.date=datevec(out.day+1+out.sec/86400);
out.enddate.year=out.date(end,1)+1950;
out.enddate.month=out.date(end,2);
out.enddate.day=out.date(end,3);
out.enddate.hour=out.date(end,4);
out.enddate.minute=out.date(end,5);
out.enddate.second=out.date(end,6); % not accounted for leap seconds from 1950 (36 in 2016, 37 in 2017)

% uso end date from filename
out.startdate_filename=filein(end-65-16:end-58-16);
out.enddate_filename=filein(end-65:end-58);
out.enddate.year_filename=filein(end-65:end-62);
out.enddate.month_filename=filein(end-61:end-60);
out.enddate.day_filename=filein(end-59:end-58);
out.enddate.hour_filename=filein(end-56:end-55);
out.enddate.minute_filename=filein(end-54:end-53);
out.enddate.second_filename=filein(end-52:end-51);

lastdatedata=datenum([filein(end-65:end-58) filein(end-56:end-51)],'yyyymmddHHMMSS');
lastdatefilename=datenum(out.enddate.year,out.enddate.month,out.enddate.day,out.enddate.hour,out.enddate.minute,out.enddate.second);
out.difflastdate_dataVSfilename=datevec(abs(lastdatefilename-lastdatedata));

%% USO range computation
mean_orbit_height=815000; nominal_USO=10e6;
out.USOrange=mean_orbit_height*out.USOmodelled/nominal_USO;

%% plotting settings
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',18);
set(0,'defaultTextFontName','Arial');
set(0,'defaultTextFontSize',1);
mida = get(0,'ScreenSize'); mida(3:4)=[1920,1080]; set(0,'defaultFigurePosition',mida);
set(0,'DefaultFigurePaperPositionMode','auto')

% Xaxisrangevalid=out.day(out.val_flag==1)+1+out.sec(out.val_flag==1)/86400+1950*365.25-14.5;
Xaxisrangevalid=out.day(out.val_flag==1)+1+out.sec(out.val_flag==1)/86400+1950*365.25;

% plot of USO frequency
namef=['S3MPC_SRAL_' mission '_USO_freq_lastdate_'  out.enddate_filename ];
[p,~] = polyfit(Xaxisrangevalid,out.USOmodelled(out.val_flag==1),1);
figure; plot(Xaxisrangevalid,out.USOmodelled(out.val_flag==1), '.-', 'linewidth',3,'markersize',6, 'DisplayName',['USO frequency modelled. Slope = ' num2str(p(1)*365.25) ' Hz/yr.']);
title([mission ' SRAL USO delta frequency modelled. Trend from '  out.startdate_filename ' to ' out.enddate_filename '.']);
xlabel('Date'); ylabel('USO frequency [Hz]');
legend('show');
datetick('x','mmmyy');
xlim([Xaxisrangevalid(1)-30 Xaxisrangevalid(end)+30]);
saveas(gcf,[fparts.pathstr '\' namef],'jpg');
saveas(gcf,[fparts.pathstr '\' namef],'fig');
close(gcf);

% plot of USO range
USOrangevalid=out.USOrange(out.val_flag==1);
namef=['S3MPC_SRAL_' mission '_USO_range_lastdate_'  out.enddate_filename ];
[p,~] = polyfit(Xaxisrangevalid,USOrangevalid,1);
figure; plot(Xaxisrangevalid,USOrangevalid*1000, '.-', 'linewidth',3,'markersize',6, 'DisplayName',['USO range. Slope = ' num2str(p(1)*365.25*1000) ' mm/yr.']);
title([mission ' SRAL USO Range. Trend from '  out.startdate_filename ' to ' out.enddate_filename '.']);
xlabel('Date'); ylabel('USO Range [mm]');
legend('show');
datetick('x','mmmyy');
xlim([Xaxisrangevalid(1)-30 Xaxisrangevalid(end)+30]);
ylim([USOrangevalid(1)*1000-0.5 USOrangevalid(end)*1000+2]);
saveas(gcf,[fparts.pathstr '\' namef],'jpg');
saveas(gcf,[fparts.pathstr '\' namef],'fig');
close(gcf);

%% save the output matlab file
save([fparts.pathstr '\' mission '_USOmatlabfile_upto' datestr(Xaxisrangevalid(end),'yyyymmmdd')],'out');

end
