function Monitoring_CAL_S3B( folderin )
%Monitoring_CAL: Routine phase S3B monitoring report
% calling: Monitoring_CAL_S3B('C:\Users\Pablo\Desktop\S3MPC\S3B\data\CAL\SR_1_CAL\CYCLE_ALL\');

load chd; load cst;

%% first read nc files & convert to matlab files

%%% make check of nc files already converted
calfolders=dir([folderin 'S3B_SR_1_CAL____*']);
calmatfiles=dir([folderin 'cal_ncread_S3B_SR_1_CAL____*.mat']);

if isempty(calmatfiles) && isempty(calfolders) % no matlab L1 files, nor nc folders --> close and return
    disp('No CAL1 L1B nc folders, nor matlab files');
    return
elseif ~isempty(calmatfiles) && ~isempty(calfolders) % matlab L1 files and nc folders --> find matches
    % list of mat L1 files
    for i_file=1:length(calmatfiles)
        matchL1.namemat(i_file,1:length(calmatfiles(i_file,1).name))=calmatfiles(i_file,1).name;
    end
    % list of isp L1 folders
    count_L1file=0;
    if ~isempty(calfolders)
        for i_file=1:length(calfolders)
            if calfolders(i_file,1).isdir
                count_L1file=count_L1file+1;
                matchL1.nameispdir(count_L1file,1:length(calfolders(i_file,1).name))=calfolders(i_file,1).name;
            end
        end
    else
        disp('No CAL1 L1B folders found')
    end
    
    % initialize L1b files list to read
    L1bread(1,:)='no L1B files left to read';
    % fill the list of L1b files to read
    count=0;
    for i_file=1:size(matchL1.nameispdir,1)
        for i_file1=1:size(calmatfiles,1)
            if strcmp(matchL1.nameispdir(i_file,1:76),matchL1.namemat(i_file1,12:75+12))
                matchL1.flag(i_file)=1;
                break;
            else
                matchL1.flag(i_file)=0;
            end
        end
        if matchL1.flag(i_file)==0
            count=count+1;
            L1bread(count,1:length(matchL1.nameispdir(i_file,:)))=matchL1.nameispdir(i_file,:);
        end
    end
    
    %%%%%%%%%% converting the non matched nc files
    disp(['Starting conversion of non matched S3 L1b CAL NetCDF files to MATLAB files from folder ' folderin]);
    for i_folder=1:size(L1bread,1)
        if isdir([folderin L1bread(i_folder,:)]) && strcmp('S3B_SR_1_CAL____20',L1bread(i_folder,1:18)) % only folders that begin as CAL L1b S3
            cal=readanyNETCDF_V1([folderin L1bread(i_folder,:) '/calibration.nc']);
            save([folderin 'cal_ncread_' L1bread(i_folder,1:94)], 'cal');
            disp(['converting from nc to mat file: ' L1bread(i_folder,1:94)]);
        end
    end
    
elseif  isempty(calmatfiles) && ~isempty(calfolders) % no matlab L1 files, but with nc folders --> read all nc files
    count=size(calfolders,1);
    for i=1:length(calfolders)
        L1bread(i,:) = calfolders(i).name;
    end
    
    %%%%%%%%%% converting all nc files
    disp(['Starting conversion of all S3 L1b CAL NetCDF files to MATLAB files from folder ' folderin]);
    for i_folder=1:size(L1bread,1)
        if isdir([folderin L1bread(i_folder,:)]) && strcmp('S3B_SR_1_CAL____20',L1bread(i_folder,1:18)) % only folders that begin as CAL L1b S3
            cal=readanyNETCDF_V1([folderin L1bread(i_folder,:) '/calibration.nc']);
            save([folderin 'cal_ncread_' L1bread(i_folder,1:94)], 'cal');
            disp(['converting from nc to mat file: ' L1bread(i_folder,1:94)]);
        end
    end
    
elseif  ~isempty(calmatfiles) && isempty(calfolders) % with matlab L1 files, but no nc folders --> read all nc files
    disp(['*No conversion needed* of S3 L1b CAL NetCDF files to MATLAB files. Only Matlab files have been found in folder ' folderin]);
    
end

dirmat=dir([folderin 'cal_ncread_S3B_SR_1_CAL____*.mat']); % resulting collection of matlab files
timemeas1=dirmat(1).name; timemeas2=dirmat(end).name;
timemeas_name = [timemeas1(28:42) '_' timemeas2(44:58)]; % for the figures filenames: first start sensing - last stop sensing
cycle_number=num2str(dirmat(end).name(81:83));

%% initialise figures format
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontName','Arial');
set(0,'defaultTextFontSize',16);
mida = get(0,'ScreenSize'); mida(3:4)=[1920,1080]; set(0,'defaultFigurePosition',mida);
set(0,'DefaultFigurePaperPositionMode','auto')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now read the variables to be monitored from the different calibration modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CAL1 LRM IQ read data

recnum=0;
for i_file=1:length(dirmat)
    
    load([folderin dirmat(i_file).name]);
    
    if isempty(cal.data.time_l1b_cal1_lrm_iq) % no data from CAL1 LRM IQ
        
        continue
        
    else
        
        recnum=recnum+1;
        
        timein = double(cal.data.time_l1b_cal1_lrm_iq)/86400+1;  % date in days from 2000
        [tv.Y, tv.M, tv.D, tv.H, tv.MN, tv.S]=datevec(timein); tv.Y=tv.Y+2000; % time vector-- correct pivot year
        lrmiq.time(recnum)=datenum(tv.Y,tv.M,tv.D,tv.H,tv.MN,tv.S); % new time from year 0
        lrmiq.lat(recnum) = double(cal.data.lat_l1b_cal1_lrm_iq) * 1e-6;    % latitude
        lrmiq.lon(recnum) = double(cal.data.lon_l1b_cal1_lrm_iq) * 1e-6;    % longitude
        
        %ku
        lrmiq.ptrdelku(recnum) = double(cal.data.diff_tx_rx_ku_l1b_cal1_lrm_iq) * 1e-4;                 % ku time delay diff_tx_rx in m
        lrmiq.ptrpowku(recnum) = 10*log10(double(cal.data.ptr_pow_ku_l1b_cal1_lrm_iq) * 1e-2);          % ku total power in dB
        lrmiq.ptrmaxpowku(recnum) = 10*log10(double(cal.data.ptr_max_ku_l1b_cal1_lrm_iq) * 1e-4);       % ku max power in dB
        lrmiq.ptrmaxdelku(recnum) = double(cal.data.dist_ptr_max_ku_l1b_cal1_lrm_iq) * 1e-6;            % ku time delay max in m
        lrmiq.ptrmaxfrqku(recnum) = double(cal.data.freq_ptr_max_ku_l1b_cal1_lrm_iq);                   % ku time delay max in Hz
        lrmiq.ptrwidthku(recnum) = double(cal.data.ptr_main_width_ku_l1b_cal1_lrm_iq);                  % ku main lobe width at -3 dB in Hz
        lrmiq.ptrsecpowku(recnum,1:60) = 10*log10(double(cal.data.ptr_sec_lobe_max_ku_l1b_cal1_lrm_iq') * 1e-4); % ku sec lobe power in dB
        lrmiq.ptrsecfrqku(recnum,1:60) = double(cal.data.ptr_sec_lobe_pos_ku_l1b_cal1_lrm_iq');         % ku sec lobe position in Hz
        lrmiq.ptrmedfrqku(recnum) = double(cal.data.freq_ptr_median_ku_l1b_cal1_lrm_iq);                % ku median position in Hz
        
        lrmiq.ptrdelku_f(recnum) = double(cal.data.flag_diff_tx_rx_ku_l1b_cal1_lrm_iq);         % ku time delay flag
        lrmiq.ptrpowku_f(recnum) = double(cal.data.flag_ptr_pow_ku_l1b_cal1_lrm_iq);         	% ku power flag
        lrmiq.ptrwidthku_f(recnum) = double(cal.data.flag_ptr_main_width_ku_l1b_cal1_lrm_iq);   % ku ptr width flag
        lrmiq.ptrsecku_f(recnum,1:60) = double(cal.data.flag_ptr_sec_lobe_ku_l1b_cal1_lrm_iq'); % ku sec lobe flag
        lrmiq.ptrmedku_f(recnum) = double(cal.data.flag_freq_ptr_median_ku_l1b_cal1_lrm_iq);    % ku median flag
        
        %c
        lrmiq.ptrdelc(recnum) = double(cal.data.diff_tx_rx_c_l1b_cal1_lrm_iq) * 1e-4;                   % c time delay diff_tx_rx in m
        lrmiq.ptrpowc(recnum) = 10*log10(double(cal.data.ptr_pow_c_l1b_cal1_lrm_iq) * 1e-2);            % c total power in dB
        lrmiq.ptrmaxpowc(recnum) = 10*log10(double(cal.data.ptr_max_c_l1b_cal1_lrm_iq) * 1e-4);         % c max power in dB
        lrmiq.ptrmaxdelc(recnum) = double(cal.data.dist_ptr_max_c_l1b_cal1_lrm_iq) * 1e-6;              % c time delay max in m
        lrmiq.ptrmaxfrqc(recnum) = double(cal.data.freq_ptr_max_c_l1b_cal1_lrm_iq);                     % c time delay max in Hz
        lrmiq.ptrwidthc(recnum) = double(cal.data.ptr_main_width_c_l1b_cal1_lrm_iq);                    % c main lobe width at -3 dB in Hz
        lrmiq.ptrsecpowc(recnum,1:60) = 10*log10(double(cal.data.ptr_sec_lobe_max_c_l1b_cal1_lrm_iq') * 1e-4); % c sec lobe power in dB
        lrmiq.ptrsecfrqc(recnum,1:60) = double(cal.data.ptr_sec_lobe_pos_c_l1b_cal1_lrm_iq');           % c sec lobe position in Hz
        lrmiq.ptrmedfrqc(recnum) = double(cal.data.freq_ptr_median_c_l1b_cal1_lrm_iq);                  % c median position in Hz
        
        lrmiq.ptrdelc_f(recnum) = double(cal.data.flag_diff_tx_rx_c_l1b_cal1_lrm_iq);           % c time delay flag
        lrmiq.ptrpowc_f(recnum) = double(cal.data.flag_ptr_pow_c_l1b_cal1_lrm_iq);              % c power flag
        lrmiq.ptrwidthc_f(recnum) = double(cal.data.flag_ptr_main_width_c_l1b_cal1_lrm_iq);     % c ptr width flag
        lrmiq.ptrsecc_f(recnum,1:60) = double(cal.data.flag_ptr_sec_lobe_c_l1b_cal1_lrm_iq');   % c sec lobe flag
        lrmiq.ptrmedc_f(recnum) = double(cal.data.flag_freq_ptr_median_c_l1b_cal1_lrm_iq);      % c median flag
        
    end
    
end

if recnum>0
    save([folderin 'S3MPC_S3B_CAL1_LRM_IQ_data_' timemeas_name],'lrmiq');   % save the data into a mat file
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Plotting CAL1 LRM IQ data');
    
    %%%%%%%% CAL1 LRM IQ calibration areas plot ------------------------------------------------
    namef=['S3MPC_S3B_CAL1_LRM_IQ_Calibration_Areas_'  timemeas_name];
    figure
    title(['S3B SRAL CAL1 LRM Calibration Areas from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1)  '.'],'fontsize', 14);
    axesm eckert4;
    framem; gridm;
    axis off
    geoshow('landareas.shp', 'FaceColor', 'black');
    geoshow(lrmiq.lat,lrmiq.lon,'displaytype','point','Marker','o','markerfacecolor','r','markersize',8)
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL1 Time Delay ------------------------------------------------
    
    % Ku
    namef=['S3MPC_S3B_CAL1_LRM_Ku_TimeDelay_'  timemeas_name];
    Xaxisrange=lrmiq.time(lrmiq.ptrdelku_f==0)/365.25;
    [p,~] = polyfit(Xaxisrange,lrmiq.ptrdelku(lrmiq.ptrdelku_f==0),1);
    stdevdat=std(lrmiq.ptrdelku(lrmiq.ptrdelku_f==0));
    meandat=mean(lrmiq.ptrdelku(lrmiq.ptrdelku_f==0));
    figure; plot(lrmiq.time(lrmiq.ptrdelku_f==0),lrmiq.ptrdelku(lrmiq.ptrdelku_f==0), 'o-','linewidth',3,'markersize',6, 'DisplayName',['Mean = ' num2str(meandat) ' m. Slope = ' num2str(p(1)*1000) ' mm/yr. Stdev= ' num2str(stdevdat*1000) ' mm.']);
    title(['S3B SRAL CAL1 LRM Ku Time Delay. Trend from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    xlabel('Date'); ylabel('Time Delay [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C
    namef=['S3MPC_S3B_CAL1_LRM_C_TimeDelay_'  timemeas_name];
    Xaxisrange=lrmiq.time(lrmiq.ptrdelc_f==0)/365.25;
    [p,~] = polyfit(Xaxisrange,lrmiq.ptrdelc(lrmiq.ptrdelc_f==0),1);
    stdevdat=std(lrmiq.ptrdelc(lrmiq.ptrdelc_f==0));
    meandat=mean(lrmiq.ptrdelc(lrmiq.ptrdelc_f==0));
    figure; plot(lrmiq.time(lrmiq.ptrdelc_f==0),lrmiq.ptrdelc(lrmiq.ptrdelc_f==0), 'o-','linewidth',3,'markersize',6, 'DisplayName',['Mean = ' num2str(meandat) ' m. Slope = ' num2str(p(1)*1000) ' mm/yr. Stdev= ' num2str(stdevdat*1000) ' mm.']);
    title(['S3B SRAL CAL1 LRM C Time Delay. Trend from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    xlabel('Date'); ylabel('Time Delay [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL1 ALL Delays (max, med, tx-rx) ------------------------------
    
    % Ku
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_Ku/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_Ku; % Hz to m
    namef=['S3MPC_S3B_CAL1_LRM_Ku_ALLDelays_'  timemeas_name];
    Xaxisrange  = lrmiq.time(lrmiq.ptrdelku_f==0)/365.25;
    Xaxisrange1 = lrmiq.time(lrmiq.ptrmedku_f==0)/365.25;
    [p_rxtx,~] = polyfit(Xaxisrange,   lrmiq.ptrdelku(lrmiq.ptrdelku_f==0),1);
    [p_max,~]  = polyfit(Xaxisrange,lrmiq.ptrmaxdelku(lrmiq.ptrdelku_f==0),1);
    [p_med,~]  = polyfit(Xaxisrange1,lrmiq.ptrmedfrqku(lrmiq.ptrmedku_f==0)*freq2dist,1);
    stdev_rxtx = std(lrmiq.ptrdelku(lrmiq.ptrdelku_f==0));
    stdev_max  = std(lrmiq.ptrmaxdelku(lrmiq.ptrdelku_f==0));
    stdev_med  = std(lrmiq.ptrmedfrqku(lrmiq.ptrmedku_f==0)*freq2dist);
    mean_rxtx  = mean(lrmiq.ptrdelku(lrmiq.ptrdelku_f==0));
    mean_max   = mean(lrmiq.ptrmaxdelku(lrmiq.ptrdelku_f==0));
    mean_med   = mean(lrmiq.ptrmedfrqku(lrmiq.ptrmedku_f==0)*freq2dist);
    figure; hold all;
    plot(lrmiq.time(lrmiq.ptrdelku_f==0), lrmiq.ptrdelku(lrmiq.ptrdelku_f==0),              'go-','linewidth',3,'markersize',6,'DisplayName',['Delay DifTrv Tx-Rx.   Mean= ' num2str(mean_rxtx) ' m. Slope=' num2str(p_rxtx(1)*1000) ' mm/yr. Stdev=' num2str(stdev_rxtx*1000) ' mm']);
    plot(lrmiq.time(lrmiq.ptrmedku_f==0),lrmiq.ptrmedfrqku(lrmiq.ptrmedku_f==0)*freq2dist, 'bo-','linewidth',3,'markersize',6,'DisplayName',['Delay PTR Median.    Mean= ' num2str(mean_med) ' m. Slope='   num2str(p_med(1)*1000) ' mm/yr. Stdev=' num2str(stdev_med*1000) ' mm']);
    plot(lrmiq.time(lrmiq.ptrdelku_f==0), lrmiq.ptrmaxdelku(lrmiq.ptrdelku_f==0),           'ro-','linewidth',3,'markersize',6,'DisplayName',['Delay PTR maximum. Mean= ' num2str(mean_max) ' m. Slope='  num2str(p_max(1)*1000) ' mm/yr. Stdev=' num2str(stdev_max*1000) ' mm']);
    title(['S3B SRAL CAL1 LRM Ku Delays. Trends from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    xlabel('Date'); ylabel('Time Delay [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_C/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_C; % Hz to m
    namef=['S3MPC_S3B_CAL1_LRM_C_ALLDelays_'  timemeas_name];
    Xaxisrange  = lrmiq.time(lrmiq.ptrdelc_f==0)/365.25;
    Xaxisrange1 = lrmiq.time(lrmiq.ptrmedc_f==0)/365.25;
    [p_rxtx,~] = polyfit(Xaxisrange,   lrmiq.ptrdelc(lrmiq.ptrdelc_f==0),1);
    [p_max,~]  = polyfit(Xaxisrange,lrmiq.ptrmaxdelc(lrmiq.ptrdelc_f==0),1);
    [p_med,~]  = polyfit(Xaxisrange1,lrmiq.ptrmedfrqc(lrmiq.ptrmedc_f==0)*freq2dist,1);
    stdev_rxtx = std(lrmiq.ptrdelc(lrmiq.ptrdelc_f==0));
    stdev_max  = std(lrmiq.ptrmaxdelc(lrmiq.ptrdelc_f==0));
    stdev_med  = std(lrmiq.ptrmedfrqc(lrmiq.ptrmedc_f==0)*freq2dist);
    mean_rxtx  = mean(lrmiq.ptrdelc(lrmiq.ptrdelc_f==0));
    mean_max   = mean(lrmiq.ptrmaxdelc(lrmiq.ptrdelc_f==0));
    mean_med   = mean(lrmiq.ptrmedfrqc(lrmiq.ptrmedc_f==0)*freq2dist);
    figure; hold all;
    plot(lrmiq.time(lrmiq.ptrdelc_f==0), lrmiq.ptrdelc(lrmiq.ptrdelc_f==0),          'go-','linewidth',3,'markersize',6,'DisplayName',['Delay DifTrv Tx-Rx.   Mean= ' num2str(mean_rxtx) ' m. Slope=' num2str(p_rxtx(1)*1000) ' mm/yr. Stdev=' num2str(stdev_rxtx*1000) ' mm']);
    plot(lrmiq.time(lrmiq.ptrmedc_f==0),lrmiq.ptrmedfrqc(lrmiq.ptrmedc_f==0)*freq2dist, 'bo-','linewidth',3,'markersize',6,'DisplayName',['Delay PTR Median.    Mean= ' num2str(mean_med) ' m. Slope='   num2str(p_med(1)*1000) ' mm/yr. Stdev=' num2str(stdev_med*1000) ' mm']);
    plot(lrmiq.time(lrmiq.ptrdelc_f==0), lrmiq.ptrmaxdelc(lrmiq.ptrdelc_f==0), 'ro-','linewidth',3,'markersize',6,'DisplayName',['Delay PTR maximum. Mean= ' num2str(mean_max) ' m. Slope='  num2str(p_max(1)*1000) ' mm/yr. Stdev=' num2str(stdev_max*1000) ' mm']);
    title(['S3B SRAL CAL1 LRM C Delays. Trends from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    xlabel('Date'); ylabel('Time Delay [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL1 Total & Max Power ------------------------------------------------
    
    % Ku
    namef=['S3MPC_S3B_CAL1_LRM_Ku_Power_' timemeas_name];
    Xaxisrange = lrmiq.time(lrmiq.ptrpowku_f==0)/365.25;
    Y1 = lrmiq.ptrpowku(lrmiq.ptrpowku_f==0);
    Y2 = lrmiq.ptrmaxpowku(lrmiq.ptrpowku_f==0);
    [p0,~] = polyfit(Xaxisrange,Y1,1);
    [p1,~] = polyfit(Xaxisrange,Y2,1);
    stdev0=std(Y1);
    stdev1=std(Y2);
    meanv0=mean(Y1);
    meanv1=mean(Y2);
    figure;
    [AX,H1,H2] = plotyy(lrmiq.time(lrmiq.ptrpowku_f==0),Y1,lrmiq.time(lrmiq.ptrpowku_f==0),Y2);
    set(get(AX(1),'Ylabel'),'String','Integrated Power [dB]');
    set(get(AX(2),'Ylabel'),'String','Maximum Power [dB]');
    set(H1,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Integrated Power. Mean= ' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
    set(H2,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Maximum Power.  Mean= ' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    title(['S3B SRAL CAL1 LRM Ku Integrated & Max Power. Trends from ' datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    xlabel('Date');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick(AX(1),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        datetick(AX(2),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim(AX(1),[lrmiq.time(1)-1 lrmiq.time(end)+1]);
        xlim(AX(2),[lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick(AX(1),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        datetick(AX(2),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim(AX(1),[lrmiq.time(1)-5 lrmiq.time(end)+5]);
        xlim(AX(2),[lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    P2Pplus=max([peak2peak(Y1) peak2peak(Y2)]);
    set(AX(1), 'YLim', [min(Y1) max(Y1)+0.5*P2Pplus]);
    set(AX(2), 'YLim', [min(Y2)-0.5*P2Pplus max(Y2)]);
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    clear('Y1','Y2','Xaxisrange');
    close(gcf);
    
    % C
    namef=['S3MPC_S3B_CAL1_LRM_C_Power_' timemeas_name];
    Xaxisrange = lrmiq.time(lrmiq.ptrpowc_f==0)/365.25;
    Y1 = lrmiq.ptrpowc(lrmiq.ptrpowc_f==0);
    Y2 = lrmiq.ptrmaxpowc(lrmiq.ptrpowc_f==0);
    [p0,~] = polyfit(Xaxisrange,Y1,1);
    [p1,~] = polyfit(Xaxisrange,Y2,1);
    stdev0=std(Y1);
    stdev1=std(Y2);
    meanv0=mean(Y1);
    meanv1=mean(Y2);
    figure;
    [AX,H1,H2] = plotyy(lrmiq.time(lrmiq.ptrpowc_f==0),Y1,lrmiq.time(lrmiq.ptrpowc_f==0),Y2);
    set(get(AX(1),'Ylabel'),'String','Integrated Power [dB]');
    set(get(AX(2),'Ylabel'),'String','Maximum Power [dB]');
    set(H1,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Integrated Power. Mean= ' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
    set(H2,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Maximum Power.  Mean= ' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    title(['S3B SRAL CAL1 LRM C Integrated & Max Power. Trends from ' datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    xlabel('Date');
    legend('show');
    set(AX(2),'XTickLabel',[]); % keep only one X axis
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick(AX(1),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        datetick(AX(2),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim(AX(1),[lrmiq.time(1)-1 lrmiq.time(end)+1]);
        xlim(AX(2),[lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick(AX(1),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        datetick(AX(2),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim(AX(1),[lrmiq.time(1)-5 lrmiq.time(end)+5]);
        xlim(AX(2),[lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    P2Pplus=max([peak2peak(Y1) peak2peak(Y2)]);
    set(AX(1), 'YLim', [min(Y1) max(Y1)+0.5*P2Pplus]);
    set(AX(2), 'YLim', [min(Y2)-0.5*P2Pplus max(Y2)]);
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    clear('Y1','Y2','Xaxisrange');
    close(gcf);
    
    %%%%%%%% CAL1 Width ------------------------------------------------
    
    % Ku
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_Ku/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_Ku; % Hz to m
    namef=['S3MPC_S3B_CAL1_LRM_Ku_Width_'  timemeas_name];
    Xaxisrange = lrmiq.time(lrmiq.ptrwidthku_f==0)/365.25;
    [p,~] = polyfit(Xaxisrange,lrmiq.ptrwidthku(lrmiq.ptrwidthku_f==0)*freq2dist,1);
    stdev0=std(lrmiq.ptrwidthku(lrmiq.ptrwidthku_f==0)*freq2dist);
    meanv0=mean(lrmiq.ptrwidthku(lrmiq.ptrwidthku_f==0)*freq2dist);
    figure; plot(lrmiq.time(lrmiq.ptrwidthku_f==0),lrmiq.ptrwidthku(lrmiq.ptrwidthku_f==0)*freq2dist, 'o-','linewidth',3,'markersize',6,'DisplayName',['Mean= ' num2str(meanv0) ' m. Slope = ' num2str(p(1)*1000) ' mm/yr. Stdev = ' num2str(stdev0*1000) ' mm.']);
    title(['S3B SRAL CAL1 LRM Ku PTR Width. Trend from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    xlabel('Date'); ylabel('PTR Width [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_C/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_C; % Hz to m
    namef=['S3MPC_S3B_CAL1_LRM_C_Width_'  timemeas_name];
    Xaxisrange = lrmiq.time(lrmiq.ptrwidthc_f==0)/365.25;
    [p,~] = polyfit(Xaxisrange,lrmiq.ptrwidthc(lrmiq.ptrwidthc_f==0)*freq2dist,1);
    stdev0=std(lrmiq.ptrwidthc(lrmiq.ptrwidthc_f==0)*freq2dist);
    meanv0=mean(lrmiq.ptrwidthc(lrmiq.ptrwidthc_f==0)*freq2dist);
    figure; plot(lrmiq.time(lrmiq.ptrwidthc_f==0),lrmiq.ptrwidthc(lrmiq.ptrwidthc_f==0)*freq2dist, 'o-','linewidth',3,'markersize',6,'DisplayName',['Mean= ' num2str(meanv0) ' m. Slope = ' num2str(p(1)*1000) ' mm/yr. Stdev = ' num2str(stdev0*1000) ' mm.']);
    title(['S3B SRAL CAL1 LRM C PTR Width. Trend from '  datestr(Xaxisrange(1)*365.25,1) ' to ' datestr(Xaxisrange(end)*365.25,1) '.']);
    xlabel('Date'); ylabel('PTR Width [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL1 sec lobes ------------------------------------------------
    
    % Ku distribution
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_Ku/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_Ku; % Hz to m
    namef=['S3MPC_S3B_CAL1_LRM_Ku_Sec_Lobes_Distrib'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 LRM Ku Secondary Lobes Distribution. From '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    for i_rec=1:recnum
        plot(lrmiq.ptrsecfrqku(i_rec,lrmiq.ptrsecku_f(i_rec,:)==0)*freq2dist,lrmiq.ptrsecpowku(i_rec,lrmiq.ptrsecku_f(i_rec,:)==0) - lrmiq.ptrmaxpowku(i_rec), 'o-','linewidth',3,'markersize',6);
    end
    xlabel('Range [m]'); ylabel('Power wrt Main Lobe [dB]');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C distribution
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_C/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_C; % Hz to m
    namef=['S3MPC_S3B_CAL1_LRM_C_Sec_Lobes_Distrib'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 LRM C Secondary Lobes Distribution. From '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    for i_rec=1:recnum
        plot(lrmiq.ptrsecfrqc(i_rec,lrmiq.ptrsecc_f(i_rec,:)==0)*freq2dist,lrmiq.ptrsecpowc(i_rec,lrmiq.ptrsecc_f(i_rec,:)==0) - lrmiq.ptrmaxpowc(i_rec), 'o-','linewidth',3,'markersize',6);
    end
    xlabel('Range [m]'); ylabel('Power wrt Main Lobe [dB]');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % Ku Power Trend
    namef=['S3MPC_S3B_CAL1_LRM_Ku_Sec_Lobes_Power_Trend'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 LRM Ku Secondary Lobes Power Trend from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    seclobepow_annualslope_ku(1:20)=0;
    seclobepow_std_ku(1:20)=0;
    for i_lobe=1:20
        if i_lobe<11
            num_lobe=i_lobe-11;
        else
            num_lobe=i_lobe-10;
        end
        Xaxisrange = lrmiq.time(lrmiq.ptrsecku_f(:,i_lobe)==0)'/365.25; 
        [p,~] = polyfit(Xaxisrange,lrmiq.ptrsecpowku(lrmiq.ptrsecku_f(:,i_lobe)==0,i_lobe) - lrmiq.ptrmaxpowku(recnum),1);
        seclobepow_annualslope_ku(i_lobe)=p(1);
        seclobepow_std_ku(i_lobe)=std(lrmiq.ptrsecpowku(lrmiq.ptrsecku_f(:,i_lobe)==0,i_lobe));
        plot(lrmiq.time(lrmiq.ptrsecku_f(:,i_lobe)==0),lrmiq.ptrsecpowku(lrmiq.ptrsecku_f(:,i_lobe)==0,i_lobe) - lrmiq.ptrmaxpowku(recnum), 'o-','linewidth',3,'markersize',6,'DisplayName',['SL ' num2str(num_lobe)]);
    end
    xlabel('Date'); ylabel('Power wrt Main Lobe [dB]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C Power Trend
    namef=['S3MPC_S3B_CAL1_LRM_C_Sec_Lobes_Power_Trend'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 LRM C Secondary Lobes Power Trend from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    seclobepow_annualslope_c(1:20)=0;
    seclobepow_std_c(1:20)=0;
    for i_lobe=1:20
        if i_lobe<11
            num_lobe=i_lobe-11;
        else
            num_lobe=i_lobe-10;
        end
        Xaxisrange = lrmiq.time(lrmiq.ptrsecc_f(:,i_lobe)==0)'/365.25; 
        [p,~] = polyfit(Xaxisrange,lrmiq.ptrsecpowc(lrmiq.ptrsecc_f(:,i_lobe)==0,i_lobe) - lrmiq.ptrmaxpowc(recnum),1);
        seclobepow_annualslope_c(i_lobe)=p(1);
        seclobepow_std_c(i_lobe)=std(lrmiq.ptrsecpowc(lrmiq.ptrsecc_f(:,i_lobe)==0,i_lobe));
%         disp(['Power --> SL C' num2str(num_lobe) ': yearslope=' num2str(p(1))]);
        plot(lrmiq.time(lrmiq.ptrsecc_f(:,i_lobe)==0),lrmiq.ptrsecpowc(lrmiq.ptrsecc_f(:,i_lobe)==0,i_lobe) - lrmiq.ptrmaxpowc(recnum), 'o-','linewidth',3,'markersize',6,'DisplayName',['SL ' num2str(num_lobe)]);
    end
    xlabel('Date'); ylabel('Power wrt Main Lobe [dB]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([lrmiq.time(1)-1 lrmiq.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([lrmiq.time(1)-5 lrmiq.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % Ku & C sec lobes POWER std & slope
    namef=['S3MPC_S3B_CAL1_LRM_Sec_Lobes_Power_std_slope'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 LRM Secondary Lobes Power Stdev & Inter-Annual Slope from '  datestr(lrmiq.time(1),1) ' to ' datestr(lrmiq.time(end),1) '.']);
    plot([(-10:-1) (1:10)], seclobepow_annualslope_ku, 'o-','linewidth',3,'markersize',6,'DisplayName','seclobe annual slope Ku');
    plot([(-10:-1) (1:10)], seclobepow_std_ku*100, 'o-','linewidth',3,'markersize',6,'DisplayName','seclobe std Ku');
    plot([(-10:-1) (1:10)], seclobepow_annualslope_c, 'o-','linewidth',3,'markersize',6,'DisplayName','seclobe annual slope C');
    plot([(-10:-1) (1:10)], seclobepow_std_c*100, 'o-','linewidth',3,'markersize',6,'DisplayName','seclobe std C');
    xlabel('Secondary Lobe Index'); ylabel('Stdev [dBx10^-^2] & Inter-Annual Slope [dB/year]');
    legend('show');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
end

%% CAL1 SAR NORM

recnum=0;
for i_file=1:length(dirmat)
    
    load([folderin dirmat(i_file).name]);
    
    if isempty(cal.data.time_l1b_cal1_sar_norm) % no data from CAL1 SAR NORM
        
        continue
        
    else
        
        recnum=recnum+1;
        
        timein = double(cal.data.time_l1b_cal1_sar_norm)/86400+1;  % date in days from 2000
        [tv.Y, tv.M, tv.D, tv.H, tv.MN, tv.S]=datevec(timein); tv.Y=tv.Y+2000; % time vector-- correct pivot year
        sar.time(recnum)=datenum(tv.Y,tv.M,tv.D,tv.H,tv.MN,tv.S); % new time from year 0
        sar.lat(recnum) = double(cal.data.lat_l1b_cal1_sar_norm) * 1e-6;    % latitude
        sar.lon(recnum) = double(cal.data.lon_l1b_cal1_sar_norm) * 1e-6;    % longitude
        
        
        %ku
        sar.ptrdelku(recnum) = double(cal.data.diff_tx_rx_ku_l1b_cal1_sar_norm) * 1e-4;                 % ku time delay diff_tx_rx in m
        sar.ptrpowku(recnum) = 10*log10(double(cal.data.ptr_pow_ku_l1b_cal1_sar_norm) * 1e-2);          % ku total power in dB
        sar.ptrmaxpowku(recnum) = 10*log10(double(cal.data.ptr_max_ku_l1b_cal1_sar_norm) * 1e-4);       % ku max power in dB
        sar.ptrmaxdelku(recnum) = double(cal.data.dist_ptr_max_ku_l1b_cal1_sar_norm) * 1e-6;            % ku time delay max in m
        sar.ptrmaxfrqku(recnum) = double(cal.data.freq_ptr_max_ku_l1b_cal1_sar_norm);                   % ku time delay max in Hz
        sar.ptrwidthku(recnum) = double(cal.data.ptr_main_width_ku_l1b_cal1_sar_norm);                  % ku main lobe width at -3 dB in Hz
        sar.ptrsecpowku(recnum,1:60) = 10*log10(double(cal.data.ptr_sec_lobe_max_ku_l1b_cal1_sar_norm') * 1e-4); % ku sec lobe power in dB
        sar.ptrsecfrqku(recnum,1:60) = double(cal.data.ptr_sec_lobe_pos_ku_l1b_cal1_sar_norm');         % ku sec lobe position in Hz
        sar.ptrmedfrqku(recnum) = double(cal.data.freq_ptr_median_ku_l1b_cal1_sar_norm);                % ku median position in Hz
        sar.burstpower(recnum,1:64) = 10*log10(double(cal.data.burst_power_cor_l1b_cal1_sar_norm') * 1e-4);% burst power array in dB
        sar.burstphase(recnum,1:64) = double(cal.data.burst_phase_cor_l1b_cal1_sar_norm') * 1e-4;       % burst phase array in radians
        
        sar.ptrdelku_f(recnum) = double(cal.data.flag_diff_tx_rx_ku_l1b_cal1_sar_norm);         % ku time delay flag
        sar.ptrpowku_f(recnum) = double(cal.data.flag_ptr_pow_ku_l1b_cal1_sar_norm);         	% ku power flag
        sar.ptrwidthku_f(recnum) = double(cal.data.flag_ptr_main_width_ku_l1b_cal1_sar_norm);   % ku ptr width flag
        sar.ptrsecku_f(recnum,1:60) = double(cal.data.flag_ptr_sec_lobe_ku_l1b_cal1_sar_norm'); % ku sec lobe flag
        sar.ptrmedku_f(recnum) = double(cal.data.flag_freq_ptr_median_ku_l1b_cal1_sar_norm);    % ku median flag
        
        %c
        sar.ptrdelc(recnum) = double(cal.data.diff_tx_rx_c_l1b_cal1_sar_norm) * 1e-4;                   % c time delay diff_tx_rx in m
        sar.ptrpowc(recnum) = 10*log10(double(cal.data.ptr_pow_c_l1b_cal1_sar_norm) * 1e-2);            % c total power in dB
        sar.ptrmaxpowc(recnum) = 10*log10(double(cal.data.ptr_max_c_l1b_cal1_sar_norm) * 1e-4);         % c max power in dB
        sar.ptrmaxdelc(recnum) = double(cal.data.dist_ptr_max_c_l1b_cal1_sar_norm) * 1e-6;              % c time delay max in m
        sar.ptrmaxfrqc(recnum) = double(cal.data.freq_ptr_max_c_l1b_cal1_sar_norm);                     % c time delay max in Hz
        sar.ptrwidthc(recnum) = double(cal.data.ptr_main_width_c_l1b_cal1_sar_norm);                    % c main lobe width at -3 dB in Hz
        sar.ptrsecpowc(recnum,1:60) = 10*log10(double(cal.data.ptr_sec_lobe_max_c_l1b_cal1_sar_norm') * 1e-4); % c sec lobe power in dB
        sar.ptrsecfrqc(recnum,1:60) = double(cal.data.ptr_sec_lobe_pos_c_l1b_cal1_sar_norm');           % c sec lobe position in Hz
        sar.ptrmedfrqc(recnum) = double(cal.data.freq_ptr_median_c_l1b_cal1_sar_norm);                  % c median position in Hz
        
        sar.ptrdelc_f(recnum) = double(cal.data.flag_diff_tx_rx_c_l1b_cal1_sar_norm);           % c time delay flag
        sar.ptrpowc_f(recnum) = double(cal.data.flag_ptr_pow_c_l1b_cal1_sar_norm);              % c power flag
        sar.ptrwidthc_f(recnum) = double(cal.data.flag_ptr_main_width_c_l1b_cal1_sar_norm);     % c ptr width flag
        sar.ptrsecc_f(recnum,1:60) = double(cal.data.flag_ptr_sec_lobe_c_l1b_cal1_sar_norm');   % c sec lobe flag
        sar.ptrmedc_f(recnum) = double(cal.data.flag_freq_ptr_median_c_l1b_cal1_sar_norm);      % c median flag
        
    end
    
end

if recnum>0
    save([folderin 'S3MPC_S3B_CAL1_SAR_NORM_data_' timemeas_name],'sar');   % save the data into a mat file
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Plotting CAL1 SAR data');
    
    %%%%%%%% CAL1 SAR NORM calibration areas plot ------------------------------------------------
    sartime=sar.time(sar.time~=0);
    namef=['S3MPC_S3B_CAL1_SAR_Calibration_Areas_'  timemeas_name];
    figure
    title(['S3B SRAL CAL1 SAR Calibration Areas from '  datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.'],'fontsize', 14);
    axesm eckert4;
    framem; gridm;
    axis off
    geoshow('landareas.shp', 'FaceColor', 'black');
    geoshow(sar.lat(sar.time~=0),sar.lon(sar.time~=0),'displaytype','point','Marker','o','markerfacecolor','r','markersize',8)
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL1 Time Delay ------------------------------------------------
    
    % Ku
    namef=['S3MPC_S3B_CAL1_SAR_Ku_TimeDelay_'  timemeas_name];
    Xaxisrange = sar.time(sar.ptrdelku_f==0)/365.25;
    [p,~] = polyfit(Xaxisrange,sar.ptrdelku(sar.ptrdelku_f==0),1);
    stdevdat=std(sar.ptrdelku(sar.ptrdelku_f==0));
    meandat=mean(sar.ptrdelku(sar.ptrdelku_f==0));
    figure; plot(sar.time(sar.ptrdelku_f==0),sar.ptrdelku(sar.ptrdelku_f==0), 'o-','linewidth',3,'markersize',6, 'DisplayName',['Mean = ' num2str(meandat) ' m. Slope = ' num2str(p(1)*1000) ' mm/yr. Stdev= ' num2str(stdevdat*1000) ' mm.']);
    title(['S3B SRAL CAL1 SAR Ku Time Delay. Trend from '  datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.']);
    xlabel('Date'); ylabel('Time Delay [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25 < 90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
   end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C
    namef=['S3MPC_S3B_CAL1_SAR_C_TimeDelay_'  timemeas_name];
    Xaxisrange = sar.time(sar.ptrdelc_f==0)/365.25;
    [p,~] = polyfit(Xaxisrange,sar.ptrdelc(sar.ptrdelc_f==0),1);
    stdevdat=std(sar.ptrdelc(sar.ptrdelc_f==0));
    meandat=mean(sar.ptrdelc(sar.ptrdelc_f==0));
    figure; plot(sar.time(sar.ptrdelc_f==0),sar.ptrdelc(sar.ptrdelc_f==0), 'o-','linewidth',3,'markersize',6, 'DisplayName',['Mean = ' num2str(meandat) ' m. Slope = ' num2str(p(1)*100) ' mm/yr. Stdev= ' num2str(stdevdat*1000) ' mm.']);
    title(['S3B SRAL CAL1 SAR C Time Delay. Trend from '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Date'); ylabel('Time Delay [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL1 ALL Delays (max, med, tx-rx) ------------------------------
    
    % Ku
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_Ku/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_Ku; % Hz to m
    namef=['S3MPC_S3B_CAL1_SAR_Ku_ALLDelays_'  timemeas_name];
    Xaxisrange  = sar.time(sar.ptrdelku_f==0)/365.25;
    Xaxisrange1 = sar.time(sar.ptrmedku_f==0)/365.25;
    [p_rxtx,~] = polyfit(Xaxisrange, sar.ptrdelku(sar.ptrdelku_f==0),1);
    [p_max,~]  = polyfit(Xaxisrange, sar.ptrmaxdelku(sar.ptrdelku_f==0),1);
    [p_med,~]  = polyfit(Xaxisrange1,sar.ptrmedfrqku(sar.ptrmedku_f==0)*freq2dist,1);
    stdev_rxtx = std(sar.ptrdelku(sar.ptrdelku_f==0));
    stdev_max  = std(sar.ptrmaxdelku(sar.ptrdelku_f==0));
    stdev_med  = std(sar.ptrmedfrqku(sar.ptrmedku_f==0)*freq2dist);
    mean_rxtx =  mean(sar.ptrdelku(sar.ptrdelku_f==0));
    mean_max  =  mean(sar.ptrmaxdelku(sar.ptrdelku_f==0));
    mean_med  =  mean(sar.ptrmedfrqku(sar.ptrmedku_f==0)*freq2dist);
    figure; hold all;
    plot(sar.time(sar.ptrdelku_f==0), sar.ptrdelku(sar.ptrdelku_f==0),             'go-','linewidth',3,'markersize',6,'DisplayName',['Delay DifTrv Tx-Rx.  Mean=' num2str(mean_rxtx) ' m. Slope=' num2str(p_rxtx(1)*1000) ' mm/yr. Stdev=' num2str(stdev_rxtx*1000) ' mm']);
    plot(sar.time(sar.ptrmedku_f==0),sar.ptrmedfrqku(sar.ptrmedku_f==0)*freq2dist, 'bo-','linewidth',3,'markersize',6,'DisplayName',['Delay PTR Median.   Mean=' num2str(mean_med) ' m. Slope='   num2str(p_med(1)*1000) ' mm/yr. Stdev=' num2str(stdev_med*1000) ' mm']);
    plot(sar.time(sar.ptrdelku_f==0), sar.ptrmaxdelku(sar.ptrdelku_f==0),          'ro-','linewidth',3,'markersize',6,'DisplayName',['Delay PTR maximum. Mean=' num2str(mean_max) ' m. Slope='  num2str(p_max(1)*1000) ' mm/yr. Stdev=' num2str(stdev_max*1000) ' mm']);
    title(['S3B SRAL CAL1 SAR Ku Delays. Trends from '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Date'); ylabel('Time Delay [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_C/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_C; % Hz to m
    namef=['S3MPC_S3B_CAL1_SAR_C_ALLDelays_'  timemeas_name];
    Xaxisrange  = sar.time(sar.ptrdelc_f==0)/365.25;
    Xaxisrange1 = sar.time(sar.ptrmedc_f==0)/365.25;
    [p_rxtx,~] = polyfit(Xaxisrange ,sar.ptrdelc(sar.ptrdelc_f==0),1);
    [p_max,~]  = polyfit(Xaxisrange ,sar.ptrmaxdelc(sar.ptrdelc_f==0),1);
    [p_med,~]  = polyfit(Xaxisrange1,sar.ptrmedfrqc(sar.ptrmedc_f==0)*freq2dist,1);
    stdev_rxtx = std(sar.ptrdelc(sar.ptrdelc_f==0));
    stdev_max  = std(sar.ptrmaxdelc(sar.ptrdelc_f==0));
    stdev_med  = std(sar.ptrmedfrqc(sar.ptrmedc_f==0)*freq2dist);
    mean_rxtx =  mean(sar.ptrdelc(sar.ptrdelc_f==0));
    mean_max  =  mean(sar.ptrmaxdelc(sar.ptrdelc_f==0));
    mean_med  =  mean(sar.ptrmedfrqc(sar.ptrmedc_f==0)*freq2dist);
    figure; hold all;
    plot(sar.time(sar.ptrdelc_f==0), sar.ptrdelc(sar.ptrdelc_f==0),          'go-','linewidth',3,'markersize',6,'DisplayName',['Delay DifTrv Tx-Rx.   Mean=' num2str(mean_rxtx) ' m. Slope=' num2str(p_rxtx(1)*1000) ' mm/yr. Stdev=' num2str(stdev_rxtx*1000) ' mm']);
    plot(sar.time(sar.ptrmedc_f==0),sar.ptrmedfrqc(sar.ptrmedc_f==0)*freq2dist, 'bo-','linewidth',3,'markersize',6,'DisplayName',['Delay PTR Median.    Mean=' num2str(mean_med) ' m. Slope='   num2str(p_med(1)*1000) ' mm/yr. Stdev=' num2str(stdev_med*1000) ' mm']);
    plot(sar.time(sar.ptrdelc_f==0), sar.ptrmaxdelc(sar.ptrdelc_f==0), 'ro-','linewidth',3,'markersize',6,'DisplayName',['Delay PTR maximum. Mean=' num2str(mean_max) ' m. Slope='  num2str(p_max(1)*1000) ' mm/yr. Stdev=' num2str(stdev_max*1000) ' mm']);
    title(['S3B SRAL CAL1 SAR C Delays. Trends from '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Date'); ylabel('Time Delay [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL1 Total & Max Power ------------------------------------------------
    
    % Ku
    namef=['S3MPC_S3B_CAL1_SAR_Ku_Power_' timemeas_name];
    Xaxisrange = sar.time(sar.ptrpowku_f==0)/365.25;
    Y1 = sar.ptrpowku(sar.ptrpowku_f==0);
    Y2 = sar.ptrmaxpowku(sar.ptrpowku_f==0);
    [p0,~] = polyfit(Xaxisrange,Y1,1);
    [p1,~] = polyfit(Xaxisrange,Y2,1);
    stdev0=std(Y1);
    stdev1=std(Y2);
    meanv0=mean(Y1);
    meanv1=mean(Y2);
    figure;
    [AX,H1,H2] = plotyy(sar.time(sar.ptrpowku_f==0),Y1,sar.time(sar.ptrpowku_f==0),Y2);
    set(get(AX(1),'Ylabel'),'String','Integrated Power [dB]');
    set(get(AX(2),'Ylabel'),'String','Maximum Power [dB]');
    set(H1,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Integrated Power. Mean=' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
    set(H2,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Maximum Power. Mean=' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    title(['S3B SRAL CAL1 SAR Ku Integrated & Max Power. Trends from ' datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Date');
    legend('show');
    set(AX(2),'XTickLabel',[]); % keep only one X axis
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick(AX(1),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        datetick(AX(2),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim(AX(1),[sar.time(1)-1 sar.time(end)+1]);
        xlim(AX(2),[sar.time(1)-1 sar.time(end)+1]);
    else
        datetick(AX(1),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        datetick(AX(2),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim(AX(1),[sar.time(1)-5 sar.time(end)+5]);
        xlim(AX(2),[sar.time(1)-5 sar.time(end)+5]);
    end
    P2Pplus=max([peak2peak(Y1) peak2peak(Y2)]);
    set(AX(1), 'YLim', [min(Y1) max(Y1)+0.5*P2Pplus]);
    set(AX(2), 'YLim', [min(Y2)-0.5*P2Pplus max(Y2)]);
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    clear('Y1','Y2','Xaxisrange');
    close(gcf);    
    
%     % power model for SAR Ku long term trend estimation
%     namef=['S3MPC_S3B_CAL1_SAR_Ku_Power_modelfitting' timemeas_name];
%     Xaxisrange = sar.time(sar.ptrpowku_f==0)/365.25;
%     Y1 = sar.ptrpowku(sar.ptrpowku_f==0);
%     reg_time=(Xaxisrange(1):1/365.25:Xaxisrange(end)); % regular time basis for further fitting
%     sarkupw_reg=interp1(sar.time/365.25,Y1,reg_time,'pchip'); % interpolating the psarku power into the new regular timming
%     f=fit((1:length(sarkupw_reg))',sarkupw_reg','power2'); % fitting as power model
%     reg_time_extended=(reg_time(1):1/365.25:reg_time(1)+ 20); % extended time axis for extrapolation: 20 years
%     pw_model_extended=f.a*(1:365.25*20)'.^f.b+f.c; % evaluation of fitting in extended time axis
%     [~,CX0]=min(abs(pw_model_extended-pw_model_extended(1)+3.5)); time4ice=reg_time_extended(CX0); 
%     topXY=min(length(reg_time_extended),length(pw_model_extended));
%     figure; hold all;
%     plot(reg_time,sarkupw_reg,'linewidth',4,'DisplayName','CAL1 SAR Ku Integrated Power');
%     plot(reg_time_extended(1:topXY),pw_model_extended(1:topXY),'linewidth',2,'DisplayName','Power Model');
%     title(['S3B SRAL CAL1 SAR Ku Integrated Power. Fitting by Power Model. Drop of 3.5dB at ' num2str(floor(time4ice)) '/' num2str(round((time4ice-floor(time4ice))*12)) '.']);
%     plot(reg_time_extended(CX0),pw_model_extended(CX0),'o','markersize',10,'markerfacecolor','r','markeredgecolor','r','DisplayName','Requirement of -3.5 dB for sea-ice & ice sheets');
%     legend('show');
%     xlabel('Date');
%     ylabel('CAL1 SAR Ku Integrated Power [dB]');  
%     saveas(gcf,[folderin namef],'jpg');
%     saveas(gcf,[folderin namef],'fig');
%     close(gcf);
%     clear('reg_time_extended','pw_model_extended','topXY','Xaxisrange','sarkupw_reg','reg_time','f','CX0');
   
    % C
    namef=['S3MPC_S3B_CAL1_SAR_C_Power_' timemeas_name];
    Xaxisrange = sar.time(sar.ptrpowc_f==0)/365.25;
    Y1 = sar.ptrpowc(sar.ptrpowc_f==0);
    Y2 = sar.ptrmaxpowc(sar.ptrpowc_f==0);
    [p0,~] = polyfit(Xaxisrange,Y1,1);
    [p1,~] = polyfit(Xaxisrange,Y2,1);
    stdev0=std(Y1);
    stdev1=std(Y2);
    meanv0=mean(Y1);
    meanv1=mean(Y2);
    figure;
    [AX,H1,H2] = plotyy(sar.time(sar.ptrpowc_f==0),Y1,sar.time(sar.ptrpowc_f==0),Y2);
    set(get(AX(1),'Ylabel'),'String','Integrated Power [dB]');
    set(get(AX(2),'Ylabel'),'String','Maximum Power [dB]');
    set(H1,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Integrated Power. Mean=' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
    set(H2,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Maximum Power. Mean=' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    title(['S3B SRAL CAL1 SAR C Integrated & Max Power. Trends from ' datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Date');
    legend('show');
    set(AX(2),'XTickLabel',[]); % keep only one X axis
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick(AX(1),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        datetick(AX(2),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim(AX(1),[sar.time(1)-1 sar.time(end)+1]);
        xlim(AX(2),[sar.time(1)-1 sar.time(end)+1]);
    else
        datetick(AX(1),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        datetick(AX(2),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim(AX(1),[sar.time(1)-5 sar.time(end)+5]);
        xlim(AX(2),[sar.time(1)-5 sar.time(end)+5]);
    end
    P2Pplus=max([peak2peak(Y1) peak2peak(Y2)]);
    set(AX(1), 'YLim', [min(Y1) max(Y1)+0.5*P2Pplus]);
    set(AX(2), 'YLim', [min(Y2)-0.5*P2Pplus max(Y2)]);
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    clear('Y1','Y2','Xaxisrange');
    close(gcf);
    
    %%%%%%%% CAL1 Width ------------------------------------------------
    
    % Ku
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_Ku/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_Ku; % Hz to m
    Xaxisrange = sar.time(sar.ptrwidthku_f==0)/365.25;
    namef=['S3MPC_S3B_CAL1_SAR_Ku_Width_'  timemeas_name];
    [p,~] = polyfit(Xaxisrange,sar.ptrwidthku(sar.ptrwidthku_f==0)*freq2dist,1);
    stdev0=std(sar.ptrwidthku(sar.ptrwidthku_f==0)*freq2dist);
    mean0=mean(sar.ptrwidthku(sar.ptrwidthku_f==0)*freq2dist);
    figure; plot(sar.time(sar.ptrwidthku_f==0),sar.ptrwidthku(sar.ptrwidthku_f==0)*freq2dist, 'o-','linewidth',3,'markersize',6,'DisplayName',['Mean=' num2str(mean0) ' m. Slope = ' num2str(p(1)*1000) ' mm/yr. Stdev = ' num2str(stdev0*1000) ' mm.']);
    title(['S3B SRAL CAL1 SAR Ku PTR Width. Trend from '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Date'); ylabel('PTR Width [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_C/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_C; % Hz to m
    Xaxisrange = sar.time(sar.ptrwidthc_f==0)/365.25;
    namef=['S3MPC_S3B_CAL1_SAR_C_Width_'  timemeas_name];
    [p,~] = polyfit(Xaxisrange,sar.ptrwidthc(sar.ptrwidthc_f==0)*freq2dist,1);
    stdev0=std(sar.ptrwidthc(sar.ptrwidthc_f==0)*freq2dist);
    mean0=mean(sar.ptrwidthc(sar.ptrwidthc_f==0)*freq2dist);
    figure; plot(sar.time(sar.ptrwidthc_f==0),sar.ptrwidthc(sar.ptrwidthc_f==0)*freq2dist, 'o-','linewidth',3,'markersize',6,'DisplayName',['Mean=' num2str(mean0) ' m. Slope = ' num2str(p(1)*1000) ' mm/yr. Stdev = ' num2str(stdev0*1000) ' mm.']);
    title(['S3B SRAL CAL1 SAR C PTR Width. Trend from '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Date'); ylabel('PTR Width [m]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL1 sec lobes ------------------------------------------------
    
    % Ku distribution
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_Ku/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_Ku; % Hz to m
    namef=['S3MPC_S3B_CAL1_SAR_Ku_Sec_Lobes_Distrib'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 SAR Ku Secondary Lobes Distribution. From '  datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.']);
    for i_rec=1:recnum
        plot(sar.ptrsecfrqku(i_rec,sar.ptrsecku_f(i_rec,:)==0)*freq2dist,sar.ptrsecpowku(i_rec,sar.ptrsecku_f(i_rec,:)==0) - sar.ptrmaxpowku(i_rec), 'o-','linewidth',3,'markersize',6);
    end
    xlabel('Range [m]'); ylabel('Power wrt Main Lobe [dB]');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C distribution
    freq2dist=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_C/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_C; % Hz to m
    namef=['S3MPC_S3B_CAL1_SAR_C_Sec_Lobes_Distrib'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 SAR C Secondary Lobes Distribution. From '  datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.']);
    for i_rec=1:recnum
        plot(sar.ptrsecfrqc(i_rec,sar.ptrsecc_f(i_rec,:)==0)*freq2dist,sar.ptrsecpowc(i_rec,sar.ptrsecc_f(i_rec,:)==0) - sar.ptrmaxpowc(i_rec), 'o-','linewidth',3,'markersize',6);
    end
    xlabel('Range [m]'); ylabel('Power wrt Main Lobe [dB]');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % Ku Power Trend
    namef=['S3MPC_S3B_CAL1_SAR_Ku_Sec_Lobes_Power_Trend'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 SAR Ku Secondary Lobes Power Trend from '  datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.']);
    seclobepow_annualslope_ku(1:20)=0;
    seclobepow_std_ku(1:20)=0;
    for i_lobe=1:20
        if i_lobe<11
            num_lobe=i_lobe-11;
        else
            num_lobe=i_lobe-10;
        end
        Xaxisrange = sar.time(sar.ptrsecku_f(:,i_lobe)==0)'/365.25; 
        [p,~] = polyfit(Xaxisrange,sar.ptrsecpowku(sar.ptrsecku_f(:,i_lobe)==0,i_lobe) - sar.ptrmaxpowku(recnum),1);
        seclobepow_annualslope_ku(i_lobe)=p(1);
        seclobepow_std_ku(i_lobe)=std(sar.ptrsecpowku(sar.ptrsecku_f(:,i_lobe)==0,i_lobe));
        plot(sar.time(sar.ptrsecku_f(:,i_lobe)==0),sar.ptrsecpowku(sar.ptrsecku_f(:,i_lobe)==0,i_lobe) - sar.ptrmaxpowku(recnum), 'o-','linewidth',3,'markersize',6,'DisplayName',['SL ' num2str(num_lobe)]);
    end
    xlabel('Date'); ylabel('Power wrt Main Lobe [dB]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C Power Trend
    namef=['S3MPC_S3B_CAL1_SAR_C_Sec_Lobes_Power_Trend'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 SAR C Secondary Lobes Power Trend from '  datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.']);
    seclobepow_annualslope_c(1:20)=0;
    seclobepow_std_c(1:20)=0;
    for i_lobe=1:20
        if i_lobe<11
            num_lobe=i_lobe-11;
        else
            num_lobe=i_lobe-10;
        end
        Xaxisrange = sar.time(sar.ptrsecku_f(:,i_lobe)==0)'/365.25; 
        [p,~] = polyfit(Xaxisrange,sar.ptrsecpowc(sar.ptrsecc_f(:,i_lobe)==0,i_lobe) - sar.ptrmaxpowc(recnum),1);
        seclobepow_annualslope_c(i_lobe)=p(1);
        seclobepow_std_c(i_lobe)=std(sar.ptrsecpowc(sar.ptrsecc_f(:,i_lobe)==0,i_lobe));
%         disp(['Power --> SL C' num2str(num_lobe) ': yearslope=' num2str(p(1))]);
        plot(sar.time(sar.ptrsecku_f(:,i_lobe)==0),sar.ptrsecpowc(sar.ptrsecc_f(:,i_lobe)==0,i_lobe) - sar.ptrmaxpowc(recnum), 'o-','linewidth',3,'markersize',6,'DisplayName',['SL ' num2str(num_lobe)]);
    end
    xlabel('Date'); ylabel('Power wrt Main Lobe [dB]');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % Ku & C sec lobes POWER std & slope
    namef=['S3MPC_S3B_CAL1_SAR_Sec_Lobes_Power_std_slope'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL1 SAR Secondary Lobes Power Stdev & Inter-Annual Slope from '  datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.']);
    plot([(-10:-1) (1:10)], seclobepow_annualslope_ku, 'o-','linewidth',3,'markersize',6,'DisplayName','seclobe annual slope Ku');
    plot([(-10:-1) (1:10)], seclobepow_std_ku*100, 'o-','linewidth',3,'markersize',6,'DisplayName','seclobe std Ku');
    plot([(-10:-1) (1:10)], seclobepow_annualslope_c, 'o-','linewidth',3,'markersize',6,'DisplayName','seclobe annual slope C');
    plot([(-10:-1) (1:10)], seclobepow_std_c*100, 'o-','linewidth',3,'markersize',6,'DisplayName','seclobe std C');
    xlabel('Secondary Lobe Index'); ylabel('Stdev [dBx10^-^2] & Inter-Annual Slope [dB/year]');
    legend('show');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% intra-burst arrays  ------------------------------------------------
    
    % 2Dcolored burst power
    namef=['S3MPC_S3B_CAL1_SAR_Burst_Power_2D_'  timemeas_name];
    indexburstpower=~(max(sar.burstpower,[],2)==0 & min(sar.burstpower,[],2)==0);
    figure; imagesc(1:64,sar.time(indexburstpower),sar.burstpower(indexburstpower,:));
    title(['S3B SRAL CAL1 SAR Ku Burst Power [dB]. From '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Pulse number in the Burst'); ylabel('Date');
    datetick('y',20,'keeplimits','keepticks') % mmmyydd for short time series
    colormap jet;
    colorbar
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % 2Dcolored burst phase
    namef=['S3MPC_S3B_CAL1_SAR_Burst_Phase_2D_'  timemeas_name];
    indexburstphase=~(max(sar.burstphase,[],2)==0 & min(sar.burstphase,[],2)==0);
    figure; imagesc(1:64,sar.time(indexburstphase),sar.burstphase(indexburstphase,:));
    title(['S3B SRAL CAL1 SAR Ku Burst Phase [rad]. From '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Pulse number in the Burst'); ylabel('Date');
    datetick('y',20,'keeplimits','keepticks') % mmmyydd for short time series
    colormap jet;
    colorbar
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % averaged burst power
    namef=['S3MPC_S3B_CAL1_SAR_Burst_Power_'  timemeas_name];
    indexburstpower=~(max(sar.burstpower,[],2)==0 & min(sar.burstpower,[],2)==0);
    [p,~] = polyfit((1:64),mean(sar.burstpower(indexburstpower,:),1),1);
    stdev0=std(mean(sar.burstpower(indexburstpower,:),1));
    figure; plot(mean(sar.burstpower(indexburstpower,:),1), 'o-','linewidth',3,'markersize',6,'DisplayName',['Burst Power Correction: Slope = ' num2str(p(1)) ' dB/pulse. Stdev = ' num2str(stdev0) ' dB.']);
    title(['S3B SRAL CAL1 SAR Ku Burst Power. Averaged from '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Pulse number in the Burst'); ylabel('Burst Power [dB]');
    legend('show', 'Location', 'Northwest');
    xlim([0 65]);
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % averaged burst phase
    namef=['S3MPC_S3B_CAL1_SAR_Burst_Phase_'  timemeas_name];
    indexburstphase=~(max(sar.burstphase,[],2)==0 & min(sar.burstphase,[],2)==0);
    [p,~] = polyfit((1:64),mean(sar.burstphase(indexburstphase,:),1),1);
    stdev0=std(mean(sar.burstphase(indexburstphase,:),1));
    figure; plot(mean(sar.burstphase(indexburstphase,:),1), 'o-','linewidth',3,'markersize',6,'DisplayName',['Burst Phase Correction: Slope = ' num2str(p(1)) ' rad/pulse. Stdev = ' num2str(stdev0) ' rad.']);
    title(['S3B SRAL CAL1 SAR Ku Burst Phase. Averaged from '  datestr(sar.time(1),1) ' to ' datestr(sar.time(end),1) '.']);
    xlabel('Pulse number in the Burst'); ylabel('Burst Phase [rad]');
    legend('show');
    xlim([0 65]);
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % standard deviation of burst power and phase arrays along the products for each pulse
    namef=['S3MPC_S3B_CAL1_SAR_Burst_Arrays_std_' timemeas_name];
    figure;
    [AX,H1,H2] = plotyy((1:64),std(sar.burstpower(indexburstpower,:),0,1)*1e3, ... % mdB
        (1:64),std(sar.burstphase(indexburstphase,:),0,1)*180/pi); % deg
    set(get(AX(1),'Ylabel'),'String','Burst Power Stdev [mdB]');
    set(get(AX(2),'Ylabel'),'String','Burst Phase Stdev [deg]');
    set(H1,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Burst Power Stdev. Mean=' num2str(mean(std(sar.burstpower(indexburstpower,:),0,1))*1e3) ' mdB.']);
    set(H2,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Burst Phase Stdev. Mean=' num2str(mean(std(sar.burstphase(indexburstphase,:),0,1))*180/pi) ' deg.']);
    title(['S3B SRAL CAL1 SAR Burst Power & Phase pulse by pulse standard deviations from ' datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.']);
    xlabel('Pulse number in the Burst');
    legend('show');
    set(AX,'Xlim',[0 65]);
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % slope of burst power and phase arrays along the time
    recnumburst=length(find(indexburstpower==1));
    pbph(1:recnumburst,1:2)=NaN;
    pbpw(1:recnumburst,1:2)=NaN;
    for i_rec=1:recnumburst
        if indexburstphase(i_rec)
            [pbph(i_rec,1:2),~] = polyfit((1:64),sar.burstphase(i_rec,:),1);
        end
        if indexburstpower(i_rec)
            [pbpw(i_rec,1:2),~] = polyfit((1:64),sar.burstpower(i_rec,:),1);
        end
    end
    namef=['S3MPC_S3B_CAL1_SAR_Burst_Arrays_slopes_' timemeas_name];
    figure;
    Xaxisrangeph=sar.time(indexburstphase)/365.25;
    Xaxisrangepw=sar.time(indexburstpower)/365.25;
    [AX,H1,H2] = plotyy(sar.time(indexburstphase),pbph(indexburstphase,1)*180/pi, ... % deg
                        sar.time(indexburstpower),pbpw(indexburstpower,1)*1e3); % mdB
    set(get(AX(1),'Ylabel'),'String','Burst Phase Slope [deg/pulse]');
    set(get(AX(2),'Ylabel'),'String','Burst Power Slope [mdB/pulse]');
    set(H1,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Burst Phase Slope. Mean=' num2str(mean(pbph(indexburstphase,1)*180/pi)) ' deg/pulse.']);
    set(H2,'linewidth',3,'marker', 'o','markersize',6,'DisplayName',['Burst Power Slope. Mean=' num2str(mean(pbpw(indexburstpower,1)*1e3)) ' mdB/pulse.']);
    title(['S3B SRAL CAL1 SAR Burst Power & Phase arrays slopes from ' datestr(sartime(1),1) ' to ' datestr(sartime(end),1) '.']);
    xlabel('Date');
    legend('show');
    if (Xaxisrange(end)-Xaxisrange(1))*365.25<90
        datetick(AX(1),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        datetick(AX(2),'x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim(AX(1),[sar.time(1)-1 sar.time(end)+1]);
        xlim(AX(2),[sar.time(1)-1 sar.time(end)+1]);
    else
        datetick(AX(1),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        datetick(AX(2),'x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim(AX(1),[sar.time(1)-5 sar.time(end)+5]);
        xlim(AX(2),[sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
end


%% CAL1 SAR Auto

recnum=0;
for i_file=1:length(dirmat)
    
    load([folderin dirmat(i_file).name]);
    
    if isempty(cal.data.time_l1b_cal1_sar_auto) % no data from CAL2 SAR Ku
        
        continue
        
    else
        
        recnum=recnum+1;
        
        timein = double(cal.data.time_l1b_cal1_sar_auto)/86400+1;                   % date in days from 2000
        [tv.Y, tv.M, tv.D, tv.H, tv.MN, tv.S]=datevec(timein); tv.Y=tv.Y+2000;      % time vector-- correct pivot year
        sarauto.time(recnum)=datenum(tv.Y,tv.M,tv.D,tv.H,tv.MN,tv.S);           % new time from year 0
        sarauto.lat(recnum) = double(cal.data.lat_l1b_cal1_sar_auto) * 1e-6;    % latitude
        sarauto.lon(recnum) = double(cal.data.lon_l1b_cal1_sar_auto) * 1e-6;    % longitude
        
        sarauto.ref_ku(1:63,recnum) = double(cal.data.agc_ref_ku_l1b_cal1_sar_auto) * 1e-2;          % AGC reference values Ku band [dB]
        sarauto.cor_ku(1:63,recnum) = double(cal.data.agc_cor_ku_l1b_cal1_sar_auto) * 1e-2;          % AGC corrected values Ku band [dB]
        sarauto.ref_c(1:63,recnum) = double(cal.data.agc_ref_c_l1b_cal1_sar_auto) * 1e-2;            % AGC reference values C band [dB]
        sarauto.cor_c(1:63,recnum) = double(cal.data.agc_cor_c_l1b_cal1_sar_auto) * 1e-2;            % AGC corrected values C band [dB]
        sarauto.switch_ku(recnum) = double(cal.data.switch_att_ku_l1b_cal1_sar_auto) * 1e-2;         % ku band ku/c switch attenuation [dB]
        sarauto.switch_c(recnum) = double(cal.data.switch_att_c_l1b_cal1_sar_auto) * 1e-2;           % ku band ku/c switch attenuation [dB]
        
    end
    
end

if recnum>0
    save([folderin 'S3MPC_S3B_CAL1_SAR_Auto_data_' timemeas_name],'sarauto');   % save the data into a mat file
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Plotting CAL1 AutoCal data');
    
    %%%%%%%% CAL1 SAR Auto calibration areas plot ------------------------------------------------
    sarautotime=sarauto.time(sarauto.time~=0);
    namef=['S3MPC_S3B_CAL1_SAR_Auto_Calibration_Areas_'  timemeas_name];
    figure
    title(['S3B SRAL CAL1 Auto-Calibration Areas from '  datestr(sarautotime(1),1) ' to ' datestr(sarautotime(end),1)  '.'],'fontsize', 14);
    axesm eckert4;
    framem; gridm;
    axis off
    geoshow('landareas.shp', 'FaceColor', 'black');
    geoshow(sarauto.lat(sarauto.time~=0),sarauto.lon(sarauto.time~=0),'displaytype','point','Marker','o','markerfacecolor','r','markersize',8)
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%% CAL1 SAR Auto Corr-Ref Ku&C ------------------------------------------------
    namef=['S3MPC_S3B_CAL1_SAR_AutoCal_corr&ref_Ku&C_'  timemeas_name];
    figure; hold all;
    plot(mean(sarauto.cor_ku(:,sarauto.time~=0) - sarauto.ref_ku(:,sarauto.time~=0),2),'linewidth',3,'marker', 'o','markersize',6,'DisplayName','Ku AutoCal');
    plot(mean(sarauto.cor_c(:,sarauto.time~=0) - sarauto.ref_c(:,sarauto.time~=0),2),  'linewidth',3,'marker', 'o','markersize',6,'DisplayName','C AutoCal');
    xlabel('Attenuation step'); ylabel('Power [dB]');
    title(['S3B SRAL AutoCAL Averaged (Corrected - Reference). From '  datestr(sarautotime(1),1) ' to ' datestr(sarautotime(end),1) '.']);
    legend('show');
    xlim([0 65]);
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%% CAL1 SAR Auto progression Ku&C ------------------------------------------------
    for i=1:63
        sarauto_prog_ku(i,:) = sarauto.cor_ku(i,:)-sarauto.cor_ku(i,1);
        sarauto_prog_c(i,:)  = sarauto.cor_c(i,:)-sarauto.cor_c(i,1);
    end
    namef=['S3MPC_S3B_CAL1_SAR_AutoCal_corr_progression_Ku_'  timemeas_name];
    figure; imagesc(sarauto.time,1:63,sarauto_prog_ku);
    title(['S3B SRAL AutoCAL Ku-band Attenuation progression [dB]. From '  datestr(sarauto.time(1),1) ' to ' datestr(sarauto.time(end),1) '.']);
    ylabel('Attenuation Step [0:63]'); xlabel('Date');
    colormap bone;
    colorbar
    if (sarauto.time(end)-sarauto.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    namef=['S3MPC_S3B_CAL1_SAR_AutoCal_corr_progression_C_'  timemeas_name];
    figure; imagesc(sarauto.time,1:63,sarauto_prog_c);
    title(['S3B SRAL AutoCAL C-band Attenuation progression [dB]. From '  datestr(sarauto.time(1),1) ' to ' datestr(sarauto.time(end),1) '.']);
    ylabel('Attenuation Step [0:63]'); xlabel('Date');
    colormap bone;
    colorbar
    if (sarauto.time(end)-sarauto.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
end


%% CAL2 SAR Ku & C

recnum=0;
for i_file=1:length(dirmat)
    
    load([folderin dirmat(i_file).name]);
    
    if isempty(cal.data.time_l1b_cal2_sar_ku) % no data from CAL2 SAR Ku
        
        continue
        
    else
        
        recnum=recnum+1;
        
        timein = double(cal.data.time_l1b_cal2_sar_ku)/86400+1;  % date in days from 2000
        [tv.Y, tv.M, tv.D, tv.H, tv.MN, tv.S]=datevec(timein); tv.Y=tv.Y+2000; % time vector-- correct pivot year
        cal2sarku.time(recnum)=datenum(tv.Y,tv.M,tv.D,tv.H,tv.MN,tv.S); % new time from year 0
        cal2sarku.lat(recnum) = double(cal.data.lat_l1b_cal2_sar_ku) * 1e-6;    % latitude
        cal2sarku.lon(recnum) = double(cal.data.lon_l1b_cal2_sar_ku) * 1e-6;    % longitude
        
        
        cal2sarku.mean_left(recnum) = double(cal.data.gprw_mean_left_l1b_cal2_sar_ku) * 1e-2;                   % GPRW left side mean [dB]
        cal2sarku.mean_right(recnum) = double(cal.data.gprw_mean_right_l1b_cal2_sar_ku) * 1e-2;                 % GPRW right side mean [dB]
        cal2sarku.stdev_left(recnum) = double(cal.data.gprw_stdev_left_l1b_cal2_sar_ku) * 1e-2;                 % GPRW left side Stdev [dB]
        cal2sarku.stdev_right(recnum) = double(cal.data.gprw_stdev_right_l1b_cal2_sar_ku) * 1e-2;               % GPRW right side Stdev [dB]
        cal2sarku.diff_left(recnum) = double(cal.data.gprw_diff_left_l1b_cal2_sar_ku) * 1e-2;                   % GPRW left side peak2peak [dB]
        cal2sarku.diff_right(recnum) = double(cal.data.gprw_diff_right_l1b_cal2_sar_ku) * 1e-2;                 % GPRW right side peak2peak [dB]
        cal2sarku.slope_left(recnum) = double(cal.data.gprw_slope_left_l1b_cal2_sar_ku) * 1e-8;                 % GPRW left  side linreg Slope [dB/Hz]
        cal2sarku.slope_right(recnum) = double(cal.data.gprw_slope_right_l1b_cal2_sar_ku) * 1e-8;               % GPRW right side linreg Slope [dB/Hz]
        cal2sarku.stdev_slope_left(recnum) = double(cal.data.gprw_stdev_slope_left_l1b_cal2_sar_ku) * 1e-2;     % GPRW left  side desloped Stdev [dB]
        cal2sarku.stdev_slope_right(recnum) = double(cal.data.gprw_stdev_slope_right_l1b_cal2_sar_ku) * 1e-2;   % GPRW right side desloped Stdev [dB]
        cal2sarku.max_loc(recnum) = double(cal.data.gprw_max_loc_l1b_cal2_sar_ku);                              % GPRW maximum [Hz]
        cal2sarku.meas(recnum,1:128) = double(cal.data.gprw_meas_l1b_cal2_sar_ku) * 1e-4;                       % normalized GPRW [FFT p.u.]
        
        cal2sarku.flag_mean_left(recnum) = double(cal.data.flag_gprw_mean_left_l1b_cal2_sar_ku);
        cal2sarku.flag_mean_right(recnum) = double(cal.data.flag_gprw_mean_right_l1b_cal2_sar_ku);
        cal2sarku.flag_stdev_left(recnum) = double(cal.data.flag_gprw_stdev_left_l1b_cal2_sar_ku);
        cal2sarku.flag_stdev_right(recnum) = double(cal.data.flag_gprw_stdev_right_l1b_cal2_sar_ku);
        cal2sarku.flag_diff_left(recnum) = double(cal.data.flag_gprw_diff_left_l1b_cal2_sar_ku);
        cal2sarku.flag_diff_right(recnum) = double(cal.data.flag_gprw_diff_right_l1b_cal2_sar_ku);
        cal2sarku.flag_slope_left(recnum) = double(cal.data.flag_gprw_slope_left_l1b_cal2_sar_ku);
        cal2sarku.flag_slope_right(recnum) = double(cal.data.flag_gprw_slope_right_l1b_cal2_sar_ku);
        cal2sarku.flag_stdev_slope_left(recnum) = double(cal.data.flag_gprw_stdev_slope_left_l1b_cal2_sar_ku);
        cal2sarku.flag_stdev_slope_right(recnum) = double(cal.data.flag_gprw_stdev_slope_right_l1b_cal2_sar_ku);
        
    end
    
end

save([folderin 'S3MPC_S3B_CAL2_SAR_Ku_data_' timemeas_name],'cal2sarku');   % save the data into a mat file

recnum=0;
for i_file=1:length(dirmat)
    
    load([folderin dirmat(i_file).name]);
    
    if isempty(cal.data.time_l1b_cal2_sar_c) % no data from CAL2 SAR C
        
        continue
        
    else
        
        recnum=recnum+1;
        
        timein = double(cal.data.time_l1b_cal2_sar_c)/86400+1;  % date in days from 2000
        [tv.Y, tv.M, tv.D, tv.H, tv.MN, tv.S]=datevec(timein); tv.Y=tv.Y+2000; % time vector-- correct pivot year
        cal2sarc.time(recnum)=datenum(tv.Y,tv.M,tv.D,tv.H,tv.MN,tv.S); % new time from year 0
        cal2sarc.lat(recnum) = double(cal.data.lat_l1b_cal2_sar_c) * 1e-6;    % latitude
        cal2sarc.lon(recnum) = double(cal.data.lon_l1b_cal2_sar_c) * 1e-6;    % longitude
        
        cal2sarc.mean_left(recnum) = double(cal.data.gprw_mean_left_l1b_cal2_sar_c) * 1e-2;                   % GPRW left side mean [dB]
        cal2sarc.mean_right(recnum) = double(cal.data.gprw_mean_right_l1b_cal2_sar_c) * 1e-2;                 % GPRW right side mean [dB]
        cal2sarc.stdev_left(recnum) = double(cal.data.gprw_stdev_left_l1b_cal2_sar_c) * 1e-2;                 % GPRW left side Stdev [dB]
        cal2sarc.stdev_right(recnum) = double(cal.data.gprw_stdev_right_l1b_cal2_sar_c) * 1e-2;               % GPRW right side Stdev [dB]
        cal2sarc.diff_left(recnum) = double(cal.data.gprw_diff_left_l1b_cal2_sar_c) * 1e-2;                   % GPRW left side peak2peak [dB]
        cal2sarc.diff_right(recnum) = double(cal.data.gprw_diff_right_l1b_cal2_sar_c) * 1e-2;                 % GPRW right side peak2peak [dB]
        cal2sarc.slope_left(recnum) = double(cal.data.gprw_slope_left_l1b_cal2_sar_c) * 1e-8;                 % GPRW left  side linreg Slope [dB/Hz]
        cal2sarc.slope_right(recnum) = double(cal.data.gprw_slope_right_l1b_cal2_sar_c) * 1e-8;               % GPRW right side linreg Slope [dB/Hz]
        cal2sarc.stdev_slope_left(recnum) = double(cal.data.gprw_stdev_slope_left_l1b_cal2_sar_c) * 1e-2;     % GPRW left  side desloped Stdev [dB]
        cal2sarc.stdev_slope_right(recnum) = double(cal.data.gprw_stdev_slope_right_l1b_cal2_sar_c) * 1e-2;   % GPRW right side desloped Stdev [dB]
        cal2sarc.max_loc(recnum) = double(cal.data.gprw_max_loc_l1b_cal2_sar_c);                              % GPRW maximum [Hz]
        cal2sarc.meas(recnum,1:128) = double(cal.data.gprw_meas_l1b_cal2_sar_c) * 1e-4;                       % normalized GPRW [FFT p.u.]
        
        cal2sarc.flag_mean_left(recnum) = double(cal.data.flag_gprw_mean_left_l1b_cal2_sar_c);
        cal2sarc.flag_mean_right(recnum) = double(cal.data.flag_gprw_mean_right_l1b_cal2_sar_c);
        cal2sarc.flag_stdev_left(recnum) = double(cal.data.flag_gprw_stdev_left_l1b_cal2_sar_c);
        cal2sarc.flag_stdev_right(recnum) = double(cal.data.flag_gprw_stdev_right_l1b_cal2_sar_c);
        cal2sarc.flag_diff_left(recnum) = double(cal.data.flag_gprw_diff_left_l1b_cal2_sar_c);
        cal2sarc.flag_diff_right(recnum) = double(cal.data.flag_gprw_diff_right_l1b_cal2_sar_c);
        cal2sarc.flag_slope_left(recnum) = double(cal.data.flag_gprw_slope_left_l1b_cal2_sar_c);
        cal2sarc.flag_slope_right(recnum) = double(cal.data.flag_gprw_slope_right_l1b_cal2_sar_c);
        cal2sarc.flag_stdev_slope_left(recnum) = double(cal.data.flag_gprw_stdev_slope_left_l1b_cal2_sar_c);
        cal2sarc.flag_stdev_slope_right(recnum) = double(cal.data.flag_gprw_stdev_slope_right_l1b_cal2_sar_c);
        
    end
end

if recnum>0
    save([folderin 'S3MPC_S3B_CAL2_SAR_C_data_' timemeas_name],'cal2sarc');   % save the data into a mat file
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Plotting CAL2 data');
    cal2kutime=cal2sarku.time(cal2sarku.time~=0);
    cal2ctime=cal2sarc.time(cal2sarc.time~=0);

    %%%%%%%% CAL2 SAR calibration areas plot -------------------------------
    namef=['S3MPC_S3B_CAL2_SAR_Calibration_Areas_'  timemeas_name];
    figure
    title(['S3B SRAL CAL2 SAR Calibration Areas from '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1)  '.'],'fontsize', 14);
    axesm eckert4;
    framem; gridm;
    axis off
    geoshow('landareas.shp', 'FaceColor', 'black');
    geoshow(cal2sarc.lat,cal2sarc.lon,'displaytype','point','Marker','o','markerfacecolor','r','markersize',8)
    geoshow(cal2sarku.lat,cal2sarku.lon,'displaytype','point','Marker','o','markerfacecolor','r','markersize',8)
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL2 waveform  ------------------------------------------------
    
    % Ku waveforms Trend
    namef=['S3MPC_S3B_CAL2_SAR_Ku_WFMs'  timemeas_name];
    figure;
    surf((1:128),cal2sarku.time',cal2sarku.meas);
    xlabel('Range bin'); ylabel('Date'); zlabel('GPRW Power. [FFT p.u.]');
    title(['S3B SRAL CAL2 SAR Ku GPRW Waveforms. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('y',20,'keeplimits','keepticks') % mmmyydd for short time series
        ylim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('y',12,'keeplimits','keepticks') % mmmyy for long time series
        ylim([sar.time(1)-5 sar.time(end)+5]);
    end
    colormap jet;
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C waveforms Trend
    namef=['S3MPC_S3B_CAL2_SAR_C_WFMs'  timemeas_name];
    figure;
    surf((1:128),cal2sarc.time',cal2sarc.meas);
    xlabel('Range bin'); ylabel('Date'); zlabel('GPRW Power. [FFT p.u.]');
    title(['S3B SRAL CAL2 SAR C GPRW Waveforms. From '  datestr(cal2ctime(1),1) ' to ' datestr(cal2ctime(end),1) '.']);
    if (cal2sarc.time(end)-cal2sarc.time(1))<90
        datetick('y',20,'keeplimits','keepticks') % mmmyydd for short time series
        ylim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('y',12,'keeplimits','keepticks') % mmmyy for long time series
        ylim([sar.time(1)-5 sar.time(end)+5]);
    end
    colormap jet;
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % 2D Ku waveforms Ripples Trend
    namef=['S3MPC_S3B_CAL2_SAR_Ku_2DWFMsripples'  timemeas_name];
    figure;
    imagesc((11:118),cal2sarku.time',cal2sarku.meas(:,11:118));
    xlabel('Range bin'); ylabel('Date'); zlabel('GPRW Power. [FFT p.u.]');
    title(['S3B SRAL CAL2 SAR Ku GPRW Waveforms Ripples (samples 11 to 118). From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('y',20,'keeplimits','keepticks') % mmmyydd for short time series
        ylim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('y',12,'keeplimits','keepticks') % mmmyy for long time series
        ylim([sar.time(1)-5 sar.time(end)+5]);
    end
    colormap jet; colorbar('eastoutside');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % 2D C waveforms Ripples Trend
    namef=['S3MPC_S3B_CAL2_SAR_C_2DWFMsripples'  timemeas_name];
    figure;
    imagesc((11:118),cal2sarc.time',cal2sarc.meas(:,11:118));
    xlabel('Range bin'); ylabel('Date'); zlabel('GPRW Power. [FFT p.u.]');
    title(['S3B SRAL CAL2 SAR C GPRW Waveforms Ripples (samples 11 to 118). From '  datestr(cal2ctime(1),1) ' to ' datestr(cal2ctime(end),1) '.']);
    if (cal2sarc.time(end)-cal2sarc.time(1))<90
        datetick('y',20,'keeplimits','keepticks') % mmmyydd for short time series
        ylim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('y',12,'keeplimits','keepticks') % mmmyy for long time series
        ylim([sar.time(1)-5 sar.time(end)+5]);
    end
    colormap jet; colorbar('eastoutside');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % Ku and C avg waveforms
    namef=['S3MPC_S3B_CAL2_SAR_Ku&C_AvgWfm'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL2 SAR Ku (blue) and C (green) GPRW averaged wfm. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    plot(mean(cal2sarku.meas,1), 'bo-','linewidth',3,'markersize',4);
    plot(mean(cal2sarc.meas,1),  'go-','linewidth',3,'markersize',4);
    xlabel('Range Bin'); ylabel('GPRW Power. [FFT p.u.]');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL2 Mean  ------------------------------------------------
    
    % Ku Right & Left Mean
    namef=['S3MPC_S3B_CAL2_SAR_Ku_Mean'  timemeas_name];
    Xaxisrangeright=cal2sarku.time(cal2sarku.flag_mean_right==1)/365.25;
    Xaxisrangeleft=cal2sarku.time(cal2sarku.flag_mean_left==1)/365.25;
    [p0,~] = polyfit(Xaxisrangeright,cal2sarku.mean_right(cal2sarku.flag_mean_right==1),1);
    [p1,~] = polyfit(Xaxisrangeleft, cal2sarku.mean_left( cal2sarku.flag_mean_left==1),1);
    stdev0=std(cal2sarku.mean_right(cal2sarku.flag_mean_right==1));
    stdev1=std(cal2sarku.mean_left( cal2sarku.flag_mean_left==1));
    meanv0=mean(cal2sarku.mean_right(cal2sarku.flag_mean_right==1));
    meanv1=mean(cal2sarku.mean_left( cal2sarku.flag_mean_left==1));
    figure; hold all;
    title(['S3B SRAL CAL2 SAR Ku. Right (blue) and Left (red) GPRW wfm sides Mean. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    plot(cal2sarku.time(cal2sarku.flag_mean_right==1),cal2sarku.mean_right(cal2sarku.flag_mean_right==1), 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku mean right. Mean=' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
    plot(cal2sarku.time(cal2sarku.flag_mean_left==1), cal2sarku.mean_left(cal2sarku.flag_mean_left==1),  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku mean left. Mean=' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    xlabel('Date'); ylabel('GPRW waveform mean. [dB]');
    legend('show');
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C Right & Left Mean
    namef=['S3MPC_S3B_CAL2_SAR_C_Mean'  timemeas_name];
    Xaxisrangeright=cal2sarku.time(cal2sarku.flag_mean_right==1)/365.25;
    Xaxisrangeleft=cal2sarku.time(cal2sarku.flag_mean_left==1)/365.25;
    [p0,~] = polyfit(Xaxisrangeright,cal2sarc.mean_right(cal2sarc.flag_mean_right==1),1);
    [p1,~] = polyfit(Xaxisrangeleft, cal2sarc.mean_left( cal2sarc.flag_mean_left==1),1);
    stdev0=std(cal2sarc.mean_right(cal2sarc.flag_mean_right==1));
    stdev1=std(cal2sarc.mean_left( cal2sarc.flag_mean_left==1));
    meanv0=mean(cal2sarc.mean_right(cal2sarc.flag_mean_right==1));
    meanv1=mean(cal2sarc.mean_left( cal2sarc.flag_mean_left==1));
    figure; hold all;
    title(['S3B SRAL CAL2 SAR C. Right (blue) and Left (red) GPRW wfm sides Mean. From '  datestr(cal2ctime(1),1) ' to ' datestr(cal2ctime(end),1) '.']);
    plot(cal2sarku.time(cal2sarku.flag_mean_right==1),cal2sarc.mean_right(cal2sarc.flag_mean_right==1), 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C mean right. Mean=' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
    plot(cal2sarku.time(cal2sarku.flag_mean_left==1), cal2sarc.mean_left( cal2sarc.flag_mean_left==1),  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C mean left. Mean=' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    xlabel('Date'); ylabel('GPRW waveform mean. [dB]');
    legend('show');
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL2 stdev  ------------------------------------------------
    
    % Ku Right & Left stdev
    namef=['S3MPC_S3B_CAL2_SAR_Ku_Stdev'  timemeas_name];
    Xaxisrangeright=cal2sarku.time(cal2sarku.flag_stdev_right==0)/365.25;
    Xaxisrangeleft =cal2sarku.time(cal2sarku.flag_stdev_left==0)/365.25;
    figure; hold all;
    title(['S3B SRAL CAL2 SAR Ku. Right (blue) and Left (red) GPRW wfm sides Stdev. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    plot(cal2sarku.time(cal2sarku.flag_stdev_right==0),cal2sarku.stdev_right(cal2sarku.flag_stdev_right==0), 'bo-','linewidth',3,'markersize',6);
    plot(cal2sarku.time(cal2sarku.flag_stdev_left==0), cal2sarku.stdev_left( cal2sarku.flag_stdev_left==0),  'ro-','linewidth',3,'markersize',6);
    xlabel('Date'); ylabel('GPRW waveform Stdev. [dB]');
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C Right & Left stdev
    namef=['S3MPC_S3B_CAL2_SAR_C_Stdev'  timemeas_name];
    Xaxisrangeright=cal2sarc.time(cal2sarc.flag_stdev_right==0)/365.25;
    Xaxisrangeleft =cal2sarc.time(cal2sarc.flag_stdev_left==0)/365.25;
    figure; hold all;
    title(['S3B SRAL CAL2 SAR C. Right (blue) and Left (red) GPRW wfm sides Stdev. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    plot(cal2sarc.time(cal2sarc.flag_stdev_right==0),cal2sarc.stdev_right(cal2sarc.flag_stdev_right==0), 'bo-','linewidth',3,'markersize',6);
    plot(cal2sarc.time(cal2sarc.flag_stdev_left==0), cal2sarc.stdev_left( cal2sarc.flag_stdev_left==0),  'ro-','linewidth',3,'markersize',6);
    xlabel('Date'); ylabel('GPRW waveform Stdev. [dB]');
    if (cal2sarc.time(end)-cal2sarc.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL2 Peak to Peak  ------------------------------------------------
    
    % Ku Right & Left stdev
    namef=['S3MPC_S3B_CAL2_SAR_Ku_Peak2peak'  timemeas_name];
    Xaxisrangeright=cal2sarku.time(cal2sarku.flag_diff_right==0)/365.25;
    Xaxisrangeleft =cal2sarku.time(cal2sarku.flag_diff_left==0)/365.25;
    figure; hold all;
    title(['S3B SRAL CAL2 SAR Ku. Right (blue) and Left (red) GPRW wfm sides Peak to Peak. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    plot(cal2sarku.time(cal2sarku.flag_diff_right==0),cal2sarku.diff_right(cal2sarku.flag_diff_right==0), 'bo-','linewidth',3,'markersize',6);
    plot(cal2sarku.time(cal2sarku.flag_diff_left==0), cal2sarku.diff_left( cal2sarku.flag_diff_left==0),  'ro-','linewidth',3,'markersize',6);
    xlabel('Date'); ylabel('GPRW waveform Peak to Peak. [dB]');
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C Right & Left stdev
    namef=['S3MPC_S3B_CAL2_SAR_C_Peak2peak'  timemeas_name];
    Xaxisrangeright=cal2sarc.time(cal2sarc.flag_diff_right==0)/365.25;
    Xaxisrangeleft =cal2sarc.time(cal2sarc.flag_diff_left==0)/365.25;
    figure; hold all;
    title(['S3B SRAL CAL2 SAR C. Right (blue) and Left (red) GPRW wfm sides Peak to Peak. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    plot(cal2sarc.time(cal2sarc.flag_diff_right==0),cal2sarc.diff_right(cal2sarc.flag_diff_right==0), 'bo-','linewidth',3,'markersize',6);
    plot(cal2sarc.time(cal2sarc.flag_diff_left==0), cal2sarc.diff_left( cal2sarc.flag_diff_left==0),  'ro-','linewidth',3,'markersize',6);
    xlabel('Date'); ylabel('GPRW waveform Peak to Peak. [dB]');
    if (cal2sarc.time(end)-cal2sarc.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL2 Slope  ------------------------------------------------
    
    % at 2018.22 there is a change in flag data due to CHD update
    % before everything was unvalid, now it changes
    % TODO: when flag is stable change from ==1 to ==0 (0 is valid) -> done
    % done at 05/06/2018
    
    % Ku Right & Left Slope
    namef=['S3MPC_S3B_CAL2_SAR_Ku_Slope'  timemeas_name];
%     [p0,~] = polyfit(cal2sarku.time(cal2sarku.flag_slope_right==0)/365.25,cal2sarku.slope_right(cal2sarku.flag_slope_right==0),1);
%     [p1,~] = polyfit(cal2sarku.time(cal2sarku.flag_slope_left==0)/365.25, cal2sarku.slope_left( cal2sarku.flag_slope_left==0),1);
%     stdev0=std(cal2sarku.slope_right(cal2sarku.flag_slope_right==0));
%     stdev1=std(cal2sarku.slope_left( cal2sarku.flag_slope_left==0));
%     meanv0=mean(cal2sarku.slope_right(cal2sarku.flag_slope_right==0));
%     meanv1=mean(cal2sarku.slope_left( cal2sarku.flag_slope_left==0));
    [p0,~] = polyfit(cal2sarku.time/365.25,cal2sarku.slope_right,1);
    [p1,~] = polyfit(cal2sarku.time/365.25, cal2sarku.slope_left,1);
    stdev0=std(cal2sarku.slope_right);
    stdev1=std(cal2sarku.slope_left);
    meanv0=mean(cal2sarku.slope_right);
    meanv1=mean(cal2sarku.slope_left);
    figure; hold all;
    title(['S3B SRAL CAL2 SAR Ku. Right (blue) and Left (red) GPRW wfm sides linear regression Slopes. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
%     plot(cal2sarku.time(cal2sarku.flag_slope_right==0)/365.25,cal2sarku.slope_right(cal2sarku.flag_slope_right==0), 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku slope right. Mean=' num2str(meanv0) ' dB/Hz. Slope=' num2str(p0(1)) ' dB/Hz/yr. Stdev=' num2str(stdev0) ' dB/Hz']);
%     plot(cal2sarku.time(cal2sarku.flag_slope_left==0)/365.25, cal2sarku.slope_left( cal2sarku.flag_slope_left==0),  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku slope left. Mean=' num2str(meanv1) ' dB/Hz. Slope=' num2str(p1(1)) ' dB/Hz/yr. Stdev=' num2str(stdev1) ' dB/Hz']);
    plot(cal2sarku.time,cal2sarku.slope_right, 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku slope right. Mean=' num2str(meanv0) ' dB/Hz. Slope=' num2str(p0(1)) ' dB/Hz/yr. Stdev=' num2str(stdev0) ' dB/Hz']);
    plot(cal2sarku.time,cal2sarku.slope_left,  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku slope left. Mean=' num2str(meanv1) ' dB/Hz. Slope=' num2str(p1(1)) ' dB/Hz/yr. Stdev=' num2str(stdev1) ' dB/Hz']);
    xlabel('Date'); ylabel('GPRW waveform slope. [dB/Hz]');
    legend('show');
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C Right & Left Slope
    namef=['S3MPC_S3B_CAL2_SAR_C_Slope'  timemeas_name];
%     [p0,~] = polyfit(cal2sarc.time(cal2sarc.flag_slope_right==0)/365.25,cal2sarc.slope_right(cal2sarc.flag_slope_right==0),1);
%     [p1,~] = polyfit(cal2sarc.time(cal2sarc.flag_slope_left==0)/365.25, cal2sarc.slope_left( cal2sarc.flag_slope_left==0),1);
%     stdev0=std(cal2sarc.slope_right(cal2sarc.flag_slope_right==0));
%     stdev1=std(cal2sarc.slope_left( cal2sarc.flag_slope_left==0));
%     meanv0=mean(cal2sarc.slope_right(cal2sarc.flag_slope_right==0));
%     meanv1=mean(cal2sarc.slope_left( cal2sarc.flag_slope_left==0));
    [p0,~] = polyfit(cal2sarc.time/365.25,cal2sarc.slope_right,0);
    [p1,~] = polyfit(cal2sarc.time/365.25, cal2sarc.slope_left,0);
    stdev0=std(cal2sarc.slope_right);
    stdev1=std(cal2sarc.slope_left);
    meanv0=mean(cal2sarc.slope_right);
    meanv1=mean(cal2sarc.slope_left);
    figure; hold all;
    title(['S3B SRAL CAL2 SAR C. Right (blue) and Left (red) GPRW wfm sides linear regression Slopes. From '  datestr(cal2ctime(1),1) ' to ' datestr(cal2ctime(end),1) '.']);
%     plot(cal2sarc.time(cal2sarc.flag_slope_right==0)/365.25,cal2sarc.slope_right(cal2sarc.flag_slope_right==0), 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C slope right. Mean=' num2str(meanv0) ' dB/Hz. Slope=' num2str(p0(1)) ' dB/Hz/yr. Stdev=' num2str(stdev0) ' dB/Hz']);
%     plot(cal2sarc.time(cal2sarc.flag_slope_left==0)/365.25, cal2sarc.slope_left( cal2sarc.flag_slope_left==0),  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C slope left. Mean=' num2str(meanv1) ' dB/Hz. Slope=' num2str(p1(1)) ' dB/Hz/yr. Stdev=' num2str(stdev1) ' dB/Hz']);
    plot(cal2sarc.time,cal2sarc.slope_right, 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C slope right. Mean=' num2str(meanv0) ' dB/Hz. Slope=' num2str(p0(1)) ' dB/Hz/yr. Stdev=' num2str(stdev0) ' dB/Hz']);
    plot(cal2sarc.time, cal2sarc.slope_left,  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C slope left. Mean=' num2str(meanv1) ' dB/Hz. Slope=' num2str(p1(1)) ' dB/Hz/yr. Stdev=' num2str(stdev1) ' dB/Hz']);
    xlabel('Date'); ylabel('GPRW waveform slope. [dB/Hz]');
    legend('show');
    if (cal2sarc.time(end)-cal2sarc.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    %%%%%%%% CAL2 Desloped Stdev  ---------------------------------------------
    
    % Ku Right & Left Desloped Stdev
    namef=['S3MPC_S3B_CAL2_SAR_Ku_DeslopedStdev'  timemeas_name];
%     [p0,~] = polyfit(cal2sarku.time(cal2sarku.flag_stdev_slope_right==0)/365.25,cal2sarku.stdev_slope_right(cal2sarku.flag_stdev_slope_right==0),1);
%     [p1,~] = polyfit(cal2sarku.time(cal2sarku.flag_stdev_slope_left==0)/365.25, cal2sarku.stdev_slope_left( cal2sarku.flag_stdev_slope_left==0),1);
%     stdev0=std(cal2sarku.stdev_slope_right(cal2sarku.flag_stdev_slope_right==0));
%     stdev1=std(cal2sarku.stdev_slope_left( cal2sarku.flag_stdev_slope_left==0));
%     meanv0=mean(cal2sarku.stdev_slope_right(cal2sarku.flag_stdev_slope_right==0));
%     meanv1=mean(cal2sarku.stdev_slope_left( cal2sarku.flag_stdev_slope_left==0));
    [p0,~] = polyfit(cal2sarku.time/365.25,cal2sarku.stdev_slope_right,1);
    [p1,~] = polyfit(cal2sarku.time/365.25, cal2sarku.stdev_slope_left,1);
    stdev0=std(cal2sarku.stdev_slope_right);
    stdev1=std(cal2sarku.stdev_slope_left);
    meanv0=mean(cal2sarku.stdev_slope_right);
    meanv1=mean(cal2sarku.stdev_slope_left);
    figure; hold all;
    title(['S3B SRAL CAL2 SAR Ku. Right (blue) and Left (red) GPRW wfm sides desloped Stdev. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
%     plot(cal2sarku.time(cal2sarku.flag_stdev_slope_right==0)/365.25,cal2sarku.stdev_slope_right(cal2sarku.flag_stdev_slope_right==0), 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku desloped Stdev right. Mean=' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
%     plot(cal2sarku.time(cal2sarku.flag_stdev_slope_left==0)/365.25, cal2sarku.stdev_slope_left( cal2sarku.flag_stdev_slope_left==0),  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku desloped Stdev left. Mean=' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    plot(cal2sarku.time,cal2sarku.stdev_slope_right, 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku desloped Stdev right. Mean=' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
    plot(cal2sarku.time, cal2sarku.stdev_slope_left,  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 Ku desloped Stdev left. Mean=' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    xlabel('Date'); ylabel('GPRW desloped wfm side Stdev. [dB]');
    legend('show');
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
    % C Right & Left Desloped Stdev
    namef=['S3MPC_S3B_CAL2_SAR_C_DeslopedStdev'  timemeas_name];
%     [p0,~] = polyfit(cal2sarc.time(cal2sarc.flag_stdev_slope_right==0)/365.25,cal2sarc.stdev_slope_right(cal2sarc.flag_stdev_slope_right==0),1);
%     [p1,~] = polyfit(cal2sarc.time(cal2sarc.flag_stdev_slope_left==0)/365.25, cal2sarc.stdev_slope_left( cal2sarc.flag_stdev_slope_left==0),1);
%     stdev0=std(cal2sarc.stdev_slope_right(cal2sarc.flag_stdev_slope_right==0));
%     stdev1=std(cal2sarc.stdev_slope_left( cal2sarc.flag_stdev_slope_left==0));
%     meanv0=mean(cal2sarc.stdev_slope_right(cal2sarc.flag_stdev_slope_right==0));
%     meanv1=mean(cal2sarc.stdev_slope_left( cal2sarc.flag_stdev_slope_left==0));
    [p0,~] = polyfit(cal2sarc.time/365.25,cal2sarc.stdev_slope_right,1);
    [p1,~] = polyfit(cal2sarc.time/365.25, cal2sarc.stdev_slope_left,1);
    stdev0=std(cal2sarc.stdev_slope_right);
    stdev1=std(cal2sarc.stdev_slope_left);
    meanv0=mean(cal2sarc.stdev_slope_right);
    meanv1=mean(cal2sarc.stdev_slope_left);
    figure; hold all;
    title(['S3B SRAL CAL2 SAR C. Right (blue) and Left (red) GPRW wfm sides desloped Stdev. From '  datestr(cal2ctime(1),1) ' to ' datestr(cal2ctime(end),1) '.']);
%     plot(cal2sarc.time(cal2sarc.flag_stdev_slope_right==0)/365.25,cal2sarc.stdev_slope_right(cal2sarc.flag_stdev_slope_right==0), 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C desloped Stdev right. Mean=' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
%     plot(cal2sarc.time(cal2sarc.flag_stdev_slope_left==0)/365.25, cal2sarc.stdev_slope_left( cal2sarc.flag_stdev_slope_left==0),  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C desloped Stdev left. Mean=' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    plot(cal2sarc.time,cal2sarc.stdev_slope_right, 'bo-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C desloped Stdev right. Mean=' num2str(meanv0) ' dB. Slope=' num2str(p0(1)) ' dB/yr. Stdev=' num2str(stdev0) ' dB']);
    plot(cal2sarc.time, cal2sarc.stdev_slope_left,  'ro-','linewidth',3,'markersize',6,'DisplayName',['CAL2 C desloped Stdev left. Mean=' num2str(meanv1) ' dB. Slope=' num2str(p1(1)) ' dB/yr. Stdev=' num2str(stdev1) ' dB']);
    xlabel('Date'); ylabel('GPRW desloped wfm side Stdev. [dB]');
    if (cal2sarc.time(end)-cal2sarc.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    legend('show');
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);

    %%%%%%%% CAL2 Max location  -----------------------------------------------
    
    % Ku & C Max location
    namef=['S3MPC_S3B_CAL2_SAR_Ku_MaxLocation'  timemeas_name];
    figure; hold all;
    title(['S3B SRAL CAL2 SAR Ku (blue) and C (green) GPRW wfm Maximum Location. From '  datestr(cal2kutime(1),1) ' to ' datestr(cal2kutime(end),1) '.']);
    plot(cal2sarku.time(cal2sarku.time~=0),cal2sarku.max_loc(cal2sarku.time~=0)/365.25, 'bo-','linewidth',3,'markersize',6);
    plot(cal2sarc.time(cal2sarc.time~=0),  cal2sarc.max_loc(cal2sarc.time~=0)/365.25,   'go-','linewidth',3,'markersize',6);
    xlabel('Date'); ylabel('GPRW Maximum Location. [Hz]');
    if (cal2sarku.time(end)-cal2sarku.time(1))<90
        datetick('x',20,'keeplimits','keepticks') % mmmyydd for short time series
        xlim([sar.time(1)-1 sar.time(end)+1]);
    else
        datetick('x',12,'keeplimits','keepticks') % mmmyy for long time series
        xlim([sar.time(1)-5 sar.time(end)+5]);
    end
    saveas(gcf,[folderin namef],'jpg');
    saveas(gcf,[folderin namef],'fig');
    close(gcf);
    
end



%% Excel writing

tabxls=[timemeas1(28:35) '_' timemeas2(44:51)];
xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {'Reporting Table'}, tabxls, 'A1');
xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {['SRAL Calibration data period: ' timemeas_name '.']}, tabxls, 'A2');
xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {'mean', 'annual slope', 'standard deviation', 'mean', 'annual slope', 'standard deviation'}, tabxls, 'C4');
xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {['BOM to C' cycle_number]}, tabxls, 'A12');

freq2distku=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_Ku/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_Ku; % Hz to m
freq2distc=cst.light_speed/2*chd.SRAL_PTR_FFT_Features.FFT_Step_Time_C/chd.SRAL_PTR_FFT_Features.FFT_Step_Freq_C; % Hz to m

[p_lrmdelku,~] = polyfit(lrmiq.time(lrmiq.ptrdelku_f==0)/365.25,lrmiq.ptrdelku(lrmiq.ptrdelku_f==0),1);
[p_lrmpowku,~] = polyfit(lrmiq.time(lrmiq.ptrpowku_f==0)/365.25,lrmiq.ptrpowku(lrmiq.ptrpowku_f==0),1);
[p_lrmwidku,~] = polyfit(lrmiq.time(lrmiq.ptrwidthku_f==0)/365.25,lrmiq.ptrwidthku(lrmiq.ptrwidthku_f==0)*freq2distku,1);
[p_lrmdelc,~] = polyfit(lrmiq.time(lrmiq.ptrdelc_f==0)/365.25,lrmiq.ptrdelc(lrmiq.ptrdelc_f==0),1);
[p_lrmpowc,~] = polyfit(lrmiq.time(lrmiq.ptrpowc_f==0)/365.25,lrmiq.ptrpowc(lrmiq.ptrpowc_f==0),1);
[p_lrmwidc,~] = polyfit(lrmiq.time(lrmiq.ptrwidthc_f==0)/365.25,lrmiq.ptrwidthc(lrmiq.ptrwidthc_f==0)*freq2distc,1);

xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {'LRM CAL1 Time Delay',mean(lrmiq.ptrdelku(lrmiq.ptrdelku_f==0)),p_lrmdelku(1)*1000,std(lrmiq.ptrdelku(lrmiq.ptrdelku_f==0))*1000,mean(lrmiq.ptrdelc(lrmiq.ptrdelc_f==0)),p_lrmdelc(1)*1000,std(lrmiq.ptrdelc(lrmiq.ptrdelc_f==0))*1000}, tabxls, 'B5');
xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {'LRM CAL1 Power',mean(lrmiq.ptrpowku(lrmiq.ptrpowku_f==0)),p_lrmpowku(1),std(lrmiq.ptrpowku(lrmiq.ptrpowku_f==0)),mean(lrmiq.ptrpowc(lrmiq.ptrpowc_f==0)),p_lrmpowc(1),std(lrmiq.ptrpowc(lrmiq.ptrpowc_f==0))}, tabxls, 'B7');
xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {'LRM CAL1 PTR Width',mean(lrmiq.ptrwidthku(lrmiq.ptrwidthku_f==0))*freq2distku,p_lrmwidku(1)*1000,std(lrmiq.ptrwidthku(lrmiq.ptrwidthku_f==0)*1000*freq2distku),mean(lrmiq.ptrwidthc(lrmiq.ptrwidthc_f==0))*freq2distc,p_lrmwidc(1)*1000,std(lrmiq.ptrwidthc(lrmiq.ptrwidthc_f==0)*1000*freq2distc)}, tabxls, 'B9');

[p_sardelku,~] = polyfit(sar.time(sar.ptrdelku_f==0)/365.25,sar.ptrdelku(sar.ptrdelku_f==0),1);
[p_sarpowku,~] = polyfit(sar.time(sar.ptrpowku_f==0)/365.25,sar.ptrpowku(sar.ptrpowku_f==0),1);
[p_sarwidku,~] = polyfit(sar.time(sar.ptrwidthku_f==0)/365.25,sar.ptrwidthku(sar.ptrwidthku_f==0)*freq2distku,1);
[p_sardelc,~] = polyfit(sar.time(sar.ptrdelc_f==0)/365.25,sar.ptrdelc(sar.ptrdelc_f==0),1);
[p_sarpowc,~] = polyfit(sar.time(sar.ptrpowc_f==0)/365.25,sar.ptrpowc(sar.ptrpowc_f==0),1);
[p_sarwidc,~] = polyfit(sar.time(sar.ptrwidthc_f==0)/365.25,sar.ptrwidthc(sar.ptrwidthc_f==0)*freq2distc,1);

xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {'SAR CAL1 Time Delay',mean(sar.ptrdelku(sar.ptrdelku_f==0)),p_sardelku(1)*1000,std(sar.ptrdelku(sar.ptrdelku_f==0))*1000,mean(sar.ptrdelc(sar.ptrdelc_f==0)),p_sardelc(1)*1000,std(sar.ptrdelc(sar.ptrdelc_f==0))*1000}, tabxls, 'B6');
xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {'SAR CAL1 Power',mean(sar.ptrpowku(sar.ptrpowku_f==0)),p_sarpowku(1),std(sar.ptrpowku(sar.ptrpowku_f==0)),mean(sar.ptrpowc(sar.ptrpowc_f==0)),p_sarpowc(1),std(sar.ptrpowc(sar.ptrpowc_f==0))}, tabxls, 'B8');
xlswrite('C:\Users\Pablo\Desktop\S3MPC\S3B\code\SRALreporting.xlsx', {'SAR CAL1 PTR Width',mean(sar.ptrwidthku(sar.ptrwidthku_f==0))*freq2distku,p_sarwidku(1)*1000,std(sar.ptrwidthku(sar.ptrwidthku_f==0)*1000*freq2distku),mean(sar.ptrwidthc(sar.ptrwidthc_f==0))*freq2distc,p_sarwidc(1)*1000,std(sar.ptrwidthc(sar.ptrwidthc_f==0)*1000*freq2distc)}, tabxls, 'B10');
        

end