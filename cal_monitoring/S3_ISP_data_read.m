function L0 = S3_ISP_data_read(filename)
%S3_ISP_read S3 ISP reading routine
%   This function reads a S3 SRAL ISP product
%   made by isardSAT 2016

disp(['reading S3 SRAL ISP product: ' filename]);

fid = fopen(filename,'r','b'); % open the input file for reading
dirf=dir(filename); filesize=dirf.bytes;
frewind(fid); % go to beg of file

tic

i_rec=0;
while ftell(fid)<filesize
    
    i_rec=i_rec+1;
    
    if mod(i_rec,10000)==0
        disp(['reading record ' num2str(i_rec) '  -------  Timing: ' datestr(now)]);
    end
    
    %% read ISP Packet Header
    % Packet ID
    L0.packet_header.version_number(i_rec) = fread(fid,1,'ubit3'); % byte 1
    L0.packet_header.type(i_rec) = fread(fid,1,'ubit1'); % byte 1
    L0.packet_header.data_field_header_flag(i_rec) = fread(fid,1,'ubit1'); % byte 1
    L0.packet_header.process_identifier(i_rec) = fread(fid,1,'ubit7'); % byte 1-2
    L0.packet_header.packet_category(i_rec) = fread(fid,1,'ubit4'); % byte 2
    % Packet Sequence Control
    L0.packet_header.grouping_flag(i_rec) = fread(fid,1,'ubit2'); % byte 3
    L0.packet_header.sequence_count(i_rec) = fread(fid,1,'ubit14'); % byte 3-4
    % Packet Length
    L0.packet_header.packet_length(i_rec) = fread(fid,1,'uint16'); % bytes 5-6
    
    %% read ISP Data Field Header
    % Data Field Header
    fread(fid,1,'ubit1'); % spare bit, byte 7
    L0.data_field_header.pus_version_number(i_rec) = fread(fid,1,'ubit3'); % byte 7
    fread(fid,1,'ubit4'); % spare bits, byte 7
    L0.data_field_header.service_type(i_rec) = fread(fid,1,'uint8'); % byte 8
    L0.data_field_header.service_subtype(i_rec) = fread(fid,1,'uint8'); % byte 9
    L0.data_field_header.destination_identifier(i_rec) = fread(fid,1,'uint8'); % byte 10
    L0.data_field_header.coarse_time(i_rec) = fread(fid,1,'uint32'); % bytes 11-14
    L0.data_field_header.fine_time(i_rec) = fread(fid,1,'ubit24'); % bytes 15-17
    fread(fid,1,'ubit7'); % spare bits, byte 18
    L0.data_field_header.time_status(i_rec) = fread(fid,1,'ubit1'); % byte 18
    
    %% read ISP Application Data
    switch L0.packet_header.packet_category(i_rec)
        case 0 % TM_ACQ  ----------------------------------------------------------------------------------------------
            %             disp(['Product Category: TM_ACQ, record ' num2str(i_rec)]);
            fseek(fid,458-18-2,'cof');
            
        case 1 % TM_ST_SAR  ----------------------------------------------------------------------------------------------
            disp(['Product Category: TM_ST_SAR, no reading record ' num2str(i_rec)]);
            fseek(fid,16926-18-2,'cof');
            
        case 2 % TM_CAL1_LRM_I2Q2  ----------------------------------------------------------------------------------------------
%             disp('Product Category: TM_CAL1_LRM_I2Q2');
%             fseek(fid,566-18-2,'cof');
            L0.cal1lrmi2q2.val(i_rec)=1; % ground flag: use it for identifying the mode validity
            L0.cal1lrmi2q2.fine_datation(i_rec) = fread(fid,1,'uint32'); % bytes 19-22 -- SRAL Fine datation (LSB = 137.5 ns) 
            fread(fid,1,'ubit7'); % byte 23, spare bits 0-6
            L0.cal1lrmi2q2.cal1_conf(i_rec) = fread(fid,1,'ubit1'); % byte 23, bit 7 -- CAL1 Configuration. 0:attenuator enabled, 1:attenuator disabled
            L0.cal1lrmi2q2.distance_instruction(i_rec) = fread(fid,1,'uint8'); % byte 24 -- Distance instruction (from 0 to 70)
            L0.cal1lrmi2q2.ncycle(i_rec) = fread(fid,1,'uint16'); % byte 25-26 -- NCYCLE
            L0.cal1lrmi2q2.cycle_number(i_rec) = fread(fid,1,'uint16'); % byte 27-28 -- Number of cycle for a given distance instruction (from 1 to NCYCLE)
            L0.cal1lrmi2q2.att1_ku(i_rec) = fread(fid,1,'uint8'); % byte 29 -- Attenuation ATT1 for the Ku band
            L0.cal1lrmi2q2.att2_ku(i_rec) = fread(fid,1,'uint8'); % byte 30 -- Attenuation ATT2 for the Ku band
            L0.cal1lrmi2q2.att1_c(i_rec) = fread(fid,1,'uint8'); % byte 31 -- Attenuation ATT1 for the C band
            L0.cal1lrmi2q2.att2_c(i_rec) = fread(fid,1,'uint8'); % byte 32 -- Attenuation ATT2 for the C band
            fseek(fid,1,'cof'); % byte 33 reserved
            fread(fid,1,'ubit4'); % byte 34, spare bits 0-3 -- Synthesizer frequency.
            L0.cal1lrmi2q2.synthesizer_freq(i_rec) = fread(fid,1,'ubit4'); % byte 34, bits 4-7 -- Synthesizer frequency. 0b1001:649MHz, 0b1010:650MHz, 0b1011:651MHz
            L0.cal1lrmi2q2.therm.TH_KU_SSPA_THR1(i_rec) = fread(fid,1,'uint8'); % byte 35 -- Thermistances telemetry. THR1 : TH_KU_SSPA_THR1
            L0.cal1lrmi2q2.therm.TH_C_SSPA_THR2(i_rec) = fread(fid,1,'uint8'); % byte 36 -- Thermistances telemetry. THR2 : TH_C_SSPA_THR2
            L0.cal1lrmi2q2.therm.TH_KU_RX_THR3(i_rec) = fread(fid,1,'uint8'); % byte 37 -- Thermistances telemetry. THR3 : TH_KU_RX_THR3
            L0.cal1lrmi2q2.therm.TH_C_RX_THR4(i_rec) = fread(fid,1,'uint8'); % byte 38 -- Thermistances telemetry. THR4 : TH_C_RX_THR4
            L0.cal1lrmi2q2.therm.TH_FI_CNG_THR5(i_rec) = fread(fid,1,'uint8'); % byte 39 -- Thermistances telemetry. THR5 : TH_FI_CNG_THR5
            L0.cal1lrmi2q2.therm.TH_X16_THR6(i_rec) = fread(fid,1,'uint8'); % byte 40 -- Thermistances telemetry. THR6 : TH_X16_THR6
            L0.cal1lrmi2q2.therm.TH_DCDC_RFU_THR7(i_rec) = fread(fid,1,'uint8'); % byte 41 -- Thermistances telemetry. THR7 : TH_DCDC_RFU_THR7
            L0.cal1lrmi2q2.therm.ref_can1(i_rec) = fread(fid,1,'uint8'); % byte 42 -- Thermistances telemetry. THR8 : Reference CAN 1
            L0.cal1lrmi2q2.therm.TH_PLL_THR8(i_rec) = fread(fid,1,'uint8'); % byte 43 -- Thermistances telemetry. THR9 : TH_PLL_THR8
            L0.cal1lrmi2q2.therm.TH_CHIRP_THR9(i_rec) = fread(fid,1,'uint8'); % byte 44 -- Thermistances telemetry. THR10 : TH_CHIRP_THR9
            L0.cal1lrmi2q2.therm.TH_FI_THR10(i_rec) = fread(fid,1,'uint8'); % byte 45 -- Thermistances telemetry. THR11 : TH_FI_THR10
            L0.cal1lrmi2q2.therm.TH_ASIC_DFFT_THR11(i_rec) = fread(fid,1,'uint8'); % byte 46 -- Thermistances telemetry. THR12 : TH_ASIC_DFFT_THR11
            L0.cal1lrmi2q2.therm.TH_DCDC_THR12(i_rec) = fread(fid,1,'uint8'); % byte 47 -- Thermistances telemetry. THR13 : TH_DCDC_THR12
            fseek(fid,1,'cof'); % byte 48 -- Thermistances telemetry. THR14 : reserved
            fseek(fid,1,'cof'); % byte 49 -- Thermistances telemetry. THR15 : reserved
            L0.cal1lrmi2q2.therm.ref_can2(i_rec) = fread(fid,1,'uint8'); % byte 50 -- Thermistances telemetry. THR16 : Reference CAN 2
            fseek(fid,1,'cof'); % byte 51 reserved
            fseek(fid,1,'cof'); % byte 52 reserved
            for i_samp=1:128
                L0.cal1lrmi2q2.wfm_ku(i_rec,i_samp) = fread(fid,1,'uint16'); % bytes 53-308 -- 128 (I2+Q2) samples for the Ku band (averaged echo): each sample on 2 bytes (LSB = 1)
            end
            for i_samp=1:128
                L0.cal1lrmi2q2.wfm_c(i_rec,i_samp) = fread(fid,1,'uint16'); % bytes 309-564 -- 128 (I2+Q2) samples for the C band (averaged echo): each sample on 2 bytes (LSB = 1)
            end
            
            
        case 3 % TM_CAL1_LRM_IQ  ----------------------------------------------------------------------------------------------
%             disp('Product Category: TM_CAL1_LRM_IQ');
%             fseek(fid,25142-18-2,'cof');
            
            L0.cal1lrmiq.val(i_rec)=1; % ground flag: use it for identifying the mode validity
            L0.cal1lrmiq.fine_datation(i_rec) = fread(fid,1,'uint32'); % bytes 19-22 -- SRAL Fine datation (LSB = 137.5 ns) 
            fread(fid,1,'ubit6'); % byte 23, spare bits 0-5
            L0.cal1lrmiq.cal1_conf_att(i_rec) = fread(fid,1,'ubit1'); % byte 23, bit 6 -- CAL1 Configuration. 0:attenuator enabled, 1:attenuator disabled
            L0.cal1lrmiq.cal1_conf_instrmode(i_rec) = fread(fid,1,'ubit1'); % byte 23, bit 7 -- CAL1 Configuration. 0:Acquisition mode, 1:LRM mode
            L0.cal1lrmiq.distance_instruction(i_rec) = fread(fid,1,'uint8'); % byte 24 -- Distance instruction (from 0 to 70)
            L0.cal1lrmiq.ncycle(i_rec) = fread(fid,1,'uint16'); % byte 25-26 -- NCYCLE
            L0.cal1lrmiq.cycle_number(i_rec) = fread(fid,1,'uint16'); % byte 27-28 -- Number of cycle for a given distance instruction (from 1 to NCYCLE)
            L0.cal1lrmiq.att1_ku(i_rec) = fread(fid,1,'uint8'); % byte 29 -- Attenuation ATT1 for the Ku band
            L0.cal1lrmiq.att2_ku(i_rec) = fread(fid,1,'uint8'); % byte 30 -- Attenuation ATT2 for the Ku band
            L0.cal1lrmiq.att1_c(i_rec) = fread(fid,1,'uint8'); % byte 31 -- Attenuation ATT1 for the C band
            L0.cal1lrmiq.att2_c(i_rec) = fread(fid,1,'uint8'); % byte 32 -- Attenuation ATT2 for the C band
            fseek(fid,1,'cof'); % byte 33 reserved
            fread(fid,1,'ubit4'); % byte 34, spare bits 0-3 -- Synthesizer frequency.
            L0.cal1lrmiq.synthesizer_freq(i_rec) = fread(fid,1,'ubit4'); % byte 34, bits 4-7 -- Synthesizer frequency. 0b1001:649MHz, 0b1010:650MHz, 0b1011:651MHz
            L0.cal1lrmiq.therm.TH_KU_SSPA_THR1(i_rec) = fread(fid,1,'uint8'); % byte 35 -- Thermistances telemetry. THR1 : TH_KU_SSPA_THR1
            L0.cal1lrmiq.therm.TH_C_SSPA_THR2(i_rec) = fread(fid,1,'uint8'); % byte 36 -- Thermistances telemetry. THR2 : TH_C_SSPA_THR2
            L0.cal1lrmiq.therm.TH_KU_RX_THR3(i_rec) = fread(fid,1,'uint8'); % byte 37 -- Thermistances telemetry. THR3 : TH_KU_RX_THR3
            L0.cal1lrmiq.therm.TH_C_RX_THR4(i_rec) = fread(fid,1,'uint8'); % byte 38 -- Thermistances telemetry. THR4 : TH_C_RX_THR4
            L0.cal1lrmiq.therm.TH_FI_CNG_THR5(i_rec) = fread(fid,1,'uint8'); % byte 39 -- Thermistances telemetry. THR5 : TH_FI_CNG_THR5
            L0.cal1lrmiq.therm.TH_X16_THR6(i_rec) = fread(fid,1,'uint8'); % byte 40 -- Thermistances telemetry. THR6 : TH_X16_THR6
            L0.cal1lrmiq.therm.TH_DCDC_RFU_THR7(i_rec) = fread(fid,1,'uint8'); % byte 41 -- Thermistances telemetry. THR7 : TH_DCDC_RFU_THR7
            L0.cal1lrmiq.therm.ref_can1(i_rec) = fread(fid,1,'uint8'); % byte 42 -- Thermistances telemetry. THR8 : Reference CAN 1
            L0.cal1lrmiq.therm.TH_PLL_THR8(i_rec) = fread(fid,1,'uint8'); % byte 43 -- Thermistances telemetry. THR9 : TH_PLL_THR8
            L0.cal1lrmiq.therm.TH_CHIRP_THR9(i_rec) = fread(fid,1,'uint8'); % byte 44 -- Thermistances telemetry. THR10 : TH_CHIRP_THR9
            L0.cal1lrmiq.therm.TH_FI_THR10(i_rec) = fread(fid,1,'uint8'); % byte 45 -- Thermistances telemetry. THR11 : TH_FI_THR10
            L0.cal1lrmiq.therm.TH_ASIC_DFFT_THR11(i_rec) = fread(fid,1,'uint8'); % byte 46 -- Thermistances telemetry. THR12 : TH_ASIC_DFFT_THR11
            L0.cal1lrmiq.therm.TH_DCDC_THR12(i_rec) = fread(fid,1,'uint8'); % byte 47 -- Thermistances telemetry. THR13 : TH_DCDC_THR12
            fseek(fid,1,'cof'); % byte 48 -- Thermistances telemetry. THR14 : reserved
            fseek(fid,1,'cof'); % byte 49 -- Thermistances telemetry. THR15 : reserved
            L0.cal1lrmiq.therm.ref_can2(i_rec) = fread(fid,1,'uint8'); % byte 50 -- Thermistances telemetry. THR16 : Reference CAN 2
            fseek(fid,1,'cof'); % byte 51 reserved
            fseek(fid,1,'cof'); % byte 52 reserved
            % bytes 53-25140 -- 98x128 (I,Q) samples, Each (I,Q) sample on 2 consecutive bytes (I : LSB and Q : MSB) with LSB=1.
            % LRM mode : 3Ku1C3Ku, 3Ku1C3Ku... (total of 84 Ku pulses and 14 C pulses). Acquisition mode : 98 Ku pulses
            if L0.cal1lrmiq.cal1_conf_instrmode(i_rec)==0 % acquisition mode
                for i_ku_echo=1:98
                    for i_samp=1:128
                        L0.cal1lrmiq.wfmacq_ku_q(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                        L0.cal1lrmiq.wfmacq_ku_i(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                    end
                end
            else % LRM mode
                i_ku_echo=0; i_c_echo=0;
                for i_rep1=1:14
                    for i_rep2=1:3 % 3 Ku
                        i_ku_echo=i_ku_echo+1;
                        for i_samp=1:128
                            L0.cal1lrmiq.wfm_ku_q(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                            L0.cal1lrmiq.wfm_ku_i(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                        end
                    end
                    i_c_echo=i_c_echo+1;
                    for i_samp=1:128 % 1 C
                        L0.cal1lrmiq.wfm_c_q(i_rec,i_c_echo,i_samp) = fread(fid,1,'int8');
                        L0.cal1lrmiq.wfm_c_i(i_rec,i_c_echo,i_samp) = fread(fid,1,'int8');
                    end
                    for i_rep2=1:3 % 3 Ku
                        i_ku_echo=i_ku_echo+1;
                        for i_samp=1:128
                            L0.cal1lrmiq.wfm_ku_q(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                            L0.cal1lrmiq.wfm_ku_i(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                        end
                    end
                end
            end
                        
            
        case 4 % TM_CAL1_SAR  ----------------------------------------------------------------------------------------------
%             disp('Product Category: TM_CAL1_SAR');
%             fseek(fid,16950-18-2,'cof');
            
            L0.cal1sar.val(i_rec)=1; % ground flag: use it for identifying the mode validity
            L0.cal1sar.fine_datation(i_rec) = fread(fid,1,'uint32'); % bytes 19-22 -- SRAL Fine datation (LSB = 137.5 ns) 
            fread(fid,1,'ubit6'); % byte 23, spare bits 0-5
            L0.cal1sar.calibration_config(i_rec) = fread(fid,1,'ubit1'); % byte 23, bit 6 -- CAL1 Configuration. 0:Normal calibration, 1:Automatic calibration
            L0.cal1sar.attenuation_config(i_rec) = fread(fid,1,'ubit1'); % byte 23, bit 7 -- CAL1 Configuration. 0:Attenuator enabled, 1:Attenuator disabled
            L0.cal1sar.number_atten_autocal(i_rec) = fread(fid,1,'uint8'); % byte 24 -- Number of attenuation for automatic calibration (from 1 to C_RDB_CAL1_SAR_ATT_NB + 1)
            L0.cal1sar.ncycle(i_rec) = fread(fid,1,'uint16'); % byte 25-26 -- NCYCLE
            L0.cal1sar.cycle_number(i_rec) = fread(fid,1,'uint16'); % byte 27-28 -- Number of cycle for a given distance instruction (from 1 to NCYCLE)
            L0.cal1sar.att1_ku(i_rec) = fread(fid,1,'uint8'); % byte 29 -- Attenuation ATT1 for the Ku band
            L0.cal1sar.att2_ku(i_rec) = fread(fid,1,'uint8'); % byte 30 -- Attenuation ATT2 for the Ku band
            L0.cal1sar.att1_c(i_rec) = fread(fid,1,'uint8'); % byte 31 -- Attenuation ATT1 for the C band
            L0.cal1sar.att2_c(i_rec) = fread(fid,1,'uint8'); % byte 32 -- Attenuation ATT2 for the C band
            fseek(fid,1,'cof'); % byte 33 reserved
            fread(fid,1,'ubit4'); % byte 34, spare bits 0-3 -- Synthesizer frequency.
            L0.cal1sar.synthesizer_freq(i_rec) = fread(fid,1,'ubit4'); % byte 34, bits 4-7 -- Synthesizer frequency. 0b1001:649MHz, 0b1010:650MHz, 0b1011:651MHz
            L0.cal1sar.therm.TH_KU_SSPA_THR1(i_rec) = fread(fid,1,'uint8'); % byte 35 -- Thermistances telemetry. THR1 : TH_KU_SSPA_THR1
            L0.cal1sar.therm.TH_C_SSPA_THR2(i_rec) = fread(fid,1,'uint8'); % byte 36 -- Thermistances telemetry. THR2 : TH_C_SSPA_THR2
            L0.cal1sar.therm.TH_KU_RX_THR3(i_rec) = fread(fid,1,'uint8'); % byte 37 -- Thermistances telemetry. THR3 : TH_KU_RX_THR3
            L0.cal1sar.therm.TH_C_RX_THR4(i_rec) = fread(fid,1,'uint8'); % byte 38 -- Thermistances telemetry. THR4 : TH_C_RX_THR4
            L0.cal1sar.therm.TH_FI_CNG_THR5(i_rec) = fread(fid,1,'uint8'); % byte 39 -- Thermistances telemetry. THR5 : TH_FI_CNG_THR5
            L0.cal1sar.therm.TH_X16_THR6(i_rec) = fread(fid,1,'uint8'); % byte 40 -- Thermistances telemetry. THR6 : TH_X16_THR6
            L0.cal1sar.therm.TH_DCDC_RFU_THR7(i_rec) = fread(fid,1,'uint8'); % byte 41 -- Thermistances telemetry. THR7 : TH_DCDC_RFU_THR7
            L0.cal1sar.therm.ref_can1(i_rec) = fread(fid,1,'uint8'); % byte 42 -- Thermistances telemetry. THR8 : Reference CAN 1
            L0.cal1sar.therm.TH_PLL_THR8(i_rec) = fread(fid,1,'uint8'); % byte 43 -- Thermistances telemetry. THR9 : TH_PLL_THR8
            L0.cal1sar.therm.TH_CHIRP_THR9(i_rec) = fread(fid,1,'uint8'); % byte 44 -- Thermistances telemetry. THR10 : TH_CHIRP_THR9
            L0.cal1sar.therm.TH_FI_THR10(i_rec) = fread(fid,1,'uint8'); % byte 45 -- Thermistances telemetry. THR11 : TH_FI_THR10
            L0.cal1sar.therm.TH_ASIC_DFFT_THR11(i_rec) = fread(fid,1,'uint8'); % byte 46 -- Thermistances telemetry. THR12 : TH_ASIC_DFFT_THR11
            L0.cal1sar.therm.TH_DCDC_THR12(i_rec) = fread(fid,1,'uint8'); % byte 47 -- Thermistances telemetry. THR13 : TH_DCDC_THR12
            fseek(fid,1,'cof'); % byte 48 -- Thermistances telemetry. THR14 : reserved
            fseek(fid,1,'cof'); % byte 49 -- Thermistances telemetry. THR15 : reserved
            L0.cal1sar.therm.ref_can2(i_rec) = fread(fid,1,'uint8'); % byte 50 -- Thermistances telemetry. THR16 : Reference CAN 2
            fseek(fid,1,'cof'); % byte 51 reserved
            L0.cal1sar.burst_number(i_rec) = fread(fid,1,'uint8'); % byte 52 -- Burst number (from 1 to 4)
            % 66x128 (I,Q) samples : 1 C pulse + 64 Ku pulses + 1 C pulse. Each (I,Q) sample on 2 consecutive bytes (I : LSB and Q : MSB) with LSB = 1
            for i_samp=1:128 % 1 C
                L0.cal1sar.wfm_c_q(i_rec,1,i_samp) = fread(fid,1,'int8');
                L0.cal1sar.wfm_c_i(i_rec,1,i_samp) = fread(fid,1,'int8');
            end
            for i_ku_echo=1:64 % 64 Ku
                for i_samp=1:128
                    L0.cal1sar.wfm_ku_q(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                    L0.cal1sar.wfm_ku_i(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                end
            end
            for i_samp=1:128 % 1 C
                L0.cal1sar.wfm_c_q(i_rec,2,i_samp) = fread(fid,1,'int8');
                L0.cal1sar.wfm_c_i(i_rec,2,i_samp) = fread(fid,1,'int8');
            end


        case 5 % TM_CAL2_SAR  ----------------------------------------------------------------------------------------------
%             disp('Product Category: TM_CAL2_SAR');
%             fseek(fid,16946-18-2,'cof');
            
            L0.cal2.val(i_rec)=1; % ground flag: use it for identifying the mode validity
            L0.cal2.fine_datation(i_rec) = fread(fid,1,'uint32'); % bytes 19-22 -- SRAL Fine datation (LSB = 137.5 ns) 
            L0.cal2.ncycle(i_rec) = fread(fid,1,'uint16'); % byte 23-24 -- NCYCLE
            L0.cal2.cycle_number(i_rec) = fread(fid,1,'uint16'); % byte 25-26 -- Number of cycle for a given distance instruction (from 1 to NCYCLE)
            L0.cal2.att1_ku(i_rec) = fread(fid,1,'uint8'); % byte 27 -- Attenuation ATT1 for the Ku band
            L0.cal2.att2_ku(i_rec) = fread(fid,1,'uint8'); % byte 28 -- Attenuation ATT2 for the Ku band
            L0.cal2.att1_c(i_rec) = fread(fid,1,'uint8'); % byte 29 -- Attenuation ATT1 for the C band
            L0.cal2.att2_c(i_rec) = fread(fid,1,'uint8'); % byte 30 -- Attenuation ATT2 for the C band
            L0.cal2.therm.TH_KU_SSPA_THR1(i_rec) = fread(fid,1,'uint8'); % byte 31 -- Thermistances telemetry. THR1 : TH_KU_SSPA_THR1
            L0.cal2.therm.TH_C_SSPA_THR2(i_rec) = fread(fid,1,'uint8'); % byte 32 -- Thermistances telemetry. THR2 : TH_C_SSPA_THR2
            L0.cal2.therm.TH_KU_RX_THR3(i_rec) = fread(fid,1,'uint8'); % byte 33 -- Thermistances telemetry. THR3 : TH_KU_RX_THR3
            L0.cal2.therm.TH_C_RX_THR4(i_rec) = fread(fid,1,'uint8'); % byte 34 -- Thermistances telemetry. THR4 : TH_C_RX_THR4
            L0.cal2.therm.TH_FI_CNG_THR5(i_rec) = fread(fid,1,'uint8'); % byte 35 -- Thermistances telemetry. THR5 : TH_FI_CNG_THR5
            L0.cal2.therm.TH_X16_THR6(i_rec) = fread(fid,1,'uint8'); % byte 36 -- Thermistances telemetry. THR6 : TH_X16_THR6
            L0.cal2.therm.TH_DCDC_RFU_THR7(i_rec) = fread(fid,1,'uint8'); % byte 37 -- Thermistances telemetry. THR7 : TH_DCDC_RFU_THR7
            L0.cal2.therm.ref_can1(i_rec) = fread(fid,1,'uint8'); % byte 38 -- Thermistances telemetry. THR8 : Reference CAN 1
            L0.cal2.therm.TH_PLL_THR8(i_rec) = fread(fid,1,'uint8'); % byte 39 -- Thermistances telemetry. THR9 : TH_PLL_THR8
            L0.cal2.therm.TH_CHIRP_THR9(i_rec) = fread(fid,1,'uint8'); % byte 40 -- Thermistances telemetry. THR10 : TH_CHIRP_THR9
            L0.cal2.therm.TH_FI_THR10(i_rec) = fread(fid,1,'uint8'); % byte 41 -- Thermistances telemetry. THR11 : TH_FI_THR10
            L0.cal2.therm.TH_ASIC_DFFT_THR11(i_rec) = fread(fid,1,'uint8'); % byte 42 -- Thermistances telemetry. THR12 : TH_ASIC_DFFT_THR11
            L0.cal2.therm.TH_DCDC_THR12(i_rec) = fread(fid,1,'uint8'); % byte 43 -- Thermistances telemetry. THR13 : TH_DCDC_THR12
            fseek(fid,1,'cof'); % byte 44 -- Thermistances telemetry. THR14 : reserved
            fseek(fid,1,'cof'); % byte 45 -- Thermistances telemetry. THR15 : reserved
            L0.cal2.therm.ref_can2(i_rec) = fread(fid,1,'uint8'); % byte 46 -- Thermistances telemetry. THR16 : Reference CAN 2
            fseek(fid,1,'cof'); % byte 47 reserved
            L0.cal2.burst_number(i_rec) = fread(fid,1,'uint8'); % byte 48 -- Burst number (from 1 to 4)
            % 66x128 (I,Q) samples : 1 C pulse + 64 Ku pulses + 1 C pulse. Each (I,Q) sample on 2 consecutive bytes (I : LSB and Q : MSB) with LSB = 1
            for i_samp=1:128 % 1 C
                L0.cal2.wfm_c_q(i_rec,1,i_samp) = fread(fid,1,'int8');
                L0.cal2.wfm_c_i(i_rec,1,i_samp) = fread(fid,1,'int8');
            end
            for i_ku_echo=1:64 % 64 Ku
                for i_samp=1:128
                    L0.cal2.wfm_ku_q(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                    L0.cal2.wfm_ku_i(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                end
            end
            for i_samp=1:128 % 1 C
                L0.cal2.wfm_c_q(i_rec,2,i_samp) = fread(fid,1,'int8');
                L0.cal2.wfm_c_i(i_rec,2,i_samp) = fread(fid,1,'int8');
            end
            
            
        case 6 % TM_ECHO_LRM  ----------------------------------------------------------------------------------------------
            %              disp(['Product Category: TM_ECHO_LRM, record ' num2str(i_rec)]);
            %              fseek(fid,588-18-2,'cof'); % jump reading this data if needed
            
            % initialise only wfms field: 10000 more records as NaN.
            if mod(i_rec,10000)==0 || i_rec==1
                disp('initialising LRM, 10000 more records');
                L0.lrm.wfm_ku(i_rec:i_rec+10000,1:128) = NaN;
            end
            
            L0.lrm.val(i_rec)=1; % ground flag: use it for identifying the mode validity
            L0.lrm.fine_datation(i_rec) = fread(fid,1,'uint32'); % bytes 19-22 -- SRAL Fine datation (LSB = 137.5 ns) 
            fseek(fid,1,'cof'); % byte 23 reserved
            L0.lrm.mode_id(i_rec) = fread(fid,1,'uint8'); % byte 24 -- LRM CL or LRM OL
            fread(fid,1,'ubit5'); % byte 25, spare bits 0-4
            L0.lrm.tracking_conf.closed_loop_gain(i_rec) = fread(fid,1,'ubit1'); % byte 25, bit 5 -- 0:nominal value, 1:nominal value with back-off
            L0.lrm.tracking_conf.acquisition(i_rec) = fread(fid,1,'ubit1'); % byte 25, bit 6 -- 0:no acquisition, 1:acquisition
            L0.lrm.tracking_conf.dem_eeprom_read_access(i_rec) = fread(fid,1,'ubit1'); % byte 25, bit 7 -- DEM EEPROM read access. 0:enabled, 1:disabled
            fread(fid,1,'ubit7'); % byte 26, spare bits 0-6
            L0.lrm.altimeter_conf.weighting_function(i_rec) = fread(fid,1,'ubit1'); % byte 26, bit 7 -- 0:enabled, 1:disabled
            L0.lrm.ho_nav(i_rec) = fread(fid,1,'uint32'); % byte 27-30 -- Distance command H0 computed with navigation information and DEM (LSB = 3.125/64 ns)
            L0.lrm.cor2_nav(i_rec) = fread(fid,1,'int16'); % byte 31-32 -- Distance command COR2 computed with navigation information and DEM (signed integer) (LSB = 3.125/1024 ns)
            L0.lrm.ho(i_rec) = fread(fid,1,'uint32'); % byte 33-36 -- Applied distance command H0 (LSB = 3.125/64 ns)
            L0.lrm.cor2(i_rec) = fread(fid,1,'int16'); % byte 37-38 -- Applied distance command COR2 (signed integer)  (LSB = 3.125/1024 ns)
            fseek(fid,1,'cof'); % byte 39 reserved
            fread(fid,1,'ubit7'); % byte 40, spare bits 0-6
            L0.lrm.loss_track(i_rec) = fread(fid,1,'ubit1'); % byte 40, bit 7 -- Loss of track criterion computed on the echo of the cycle (N-2) in OL mode (current cycle = N). 0:normal, 1:loss of track
            L0.lrm.dist_error_ol(i_rec) = fread(fid,1,'int32'); % byte 41-44 -- Distance error computed on the echo of the cycle (N-2) in OL mode (current cycle = N)  (signed integer)  (LSB = 3.125/64 ns)
            L0.lrm.att_code_ku(i_rec) = fread(fid,1,'uint8'); % byte 45 -- Ku band ATTCODE (LSB = 1 dB)
            L0.lrm.att_code_c(i_rec) = fread(fid,1,'uint8'); % byte 46 -- C band ATTCODE (LSB = 1 dB)
            fread(fid,1,'ubit7'); % byte 47, spare bits 0-6
            L0.lrm.nav_boletin_status(i_rec) = fread(fid,1,'ubit1'); % byte 47, bit 7 -- Navigation bulletin status. 0:Bulletin OK, 1:Bulletin KO
            fread(fid,1,'ubit6'); % byte 48, spare bits 0-5
            L0.lrm.nav_boletin_source(i_rec) = fread(fid,1,'ubit2'); % byte 48, bits 6-7 -- Navigation bulletin source identifier. 0b01:GNSS, 0b10:DORIS
            L0.lrm.nav_coarse_time(i_rec) = fread(fid,1,'uint32'); % bytes 49-52 -- Time of the current navigation bulletin: Navigation bulletin coarse time on 4 bytes (LSB = 1 s)
            L0.lrm.nav_fine_time(i_rec) = fread(fid,1,'ubit24'); % bytes 53-55 -- Time of the current navigation bulletin: Navigation bulletin fine time on 3 bytes (LSB = 2^-24s)
            fseek(fid,1,'cof'); % byte 56 reserved
            L0.lrm.therm.TH_KU_SSPA_THR1(i_rec) = fread(fid,1,'uint8'); % byte 57 -- Thermistances telemetry. THR1 : TH_KU_SSPA_THR1
            L0.lrm.therm.TH_C_SSPA_THR2(i_rec) = fread(fid,1,'uint8'); % byte 58 -- Thermistances telemetry. THR2 : TH_C_SSPA_THR2
            L0.lrm.therm.TH_KU_RX_THR3(i_rec) = fread(fid,1,'uint8'); % byte 59 -- Thermistances telemetry. THR3 : TH_KU_RX_THR3
            L0.lrm.therm.TH_C_RX_THR4(i_rec) = fread(fid,1,'uint8'); % byte 60 -- Thermistances telemetry. THR4 : TH_C_RX_THR4
            L0.lrm.therm.TH_FI_CNG_THR5(i_rec) = fread(fid,1,'uint8'); % byte 61 -- Thermistances telemetry. THR5 : TH_FI_CNG_THR5
            L0.lrm.therm.TH_X16_THR6(i_rec) = fread(fid,1,'uint8'); % byte 62 -- Thermistances telemetry. THR6 : TH_X16_THR6
            L0.lrm.therm.TH_DCDC_RFU_THR7(i_rec) = fread(fid,1,'uint8'); % byte 63 -- Thermistances telemetry. THR7 : TH_DCDC_RFU_THR7
            L0.lrm.therm.ref_can1(i_rec) = fread(fid,1,'uint8'); % byte 64 -- Thermistances telemetry. THR8 : Reference CAN 1
            L0.lrm.therm.TH_PLL_THR8(i_rec) = fread(fid,1,'uint8'); % byte 65 -- Thermistances telemetry. THR9 : TH_PLL_THR8
            L0.lrm.therm.TH_CHIRP_THR9(i_rec) = fread(fid,1,'uint8'); % byte 66 -- Thermistances telemetry. THR10 : TH_CHIRP_THR9
            L0.lrm.therm.TH_FI_THR10(i_rec) = fread(fid,1,'uint8'); % byte 67 -- Thermistances telemetry. THR11 : TH_FI_THR10
            L0.lrm.therm.TH_ASIC_DFFT_THR11(i_rec) = fread(fid,1,'uint8'); % byte 68 -- Thermistances telemetry. THR12 : TH_ASIC_DFFT_THR11
            L0.lrm.therm.TH_DCDC_THR12(i_rec) = fread(fid,1,'uint8'); % byte 69 -- Thermistances telemetry. THR13 : TH_DCDC_THR12
            fseek(fid,1,'cof'); % byte 70 -- Thermistances telemetry. THR14 : reserved
            fseek(fid,1,'cof'); % byte 71 -- Thermistances telemetry. THR15 : reserved
            L0.lrm.therm.ref_can2(i_rec) = fread(fid,1,'uint8'); % byte 72 -- Thermistances telemetry. THR16 : Reference CAN 2
            fseek(fid,2,'cof'); % byte 73-74 reserved
            for i_samp=1:128
                L0.lrm.wfm_ku(i_rec,i_samp) = fread(fid,1,'uint16'); % bytes 75-330 -- 128 (I2+Q2) samples for the Ku band (averaged echo): each sample on 2 bytes (LSB = 1)
            end
            for i_samp=1:128
                L0.lrm.wfm_c(i_rec,i_samp) = fread(fid,1,'uint16'); % bytes 331-586 -- 128 (I2+Q2) samples for the C band (averaged echo): each sample on 2 bytes (LSB = 1)
            end
            
        case 7 % TM_ECHO_SAR ----------------------------------------------------------------------------------------------
            %             disp(['Product Category: TM_ECHO_SAR, record ' num2str(i_rec)]);
            %             fseek(fid,17484-18-2,'cof');
            
            % initialise only wfms field: 10000 more records as NaN.
            if mod(i_rec,10000)==0 || i_rec==1
                disp('initialising SAR, 10000 more records');
                L0.lrm.wfm_ku(i_rec:i_rec+10000,1:128) = NaN;
                L0.sar.wfm_c_q(i_rec:i_rec+10000,2,1:128) = NaN;
                L0.sar.wfm_c_i(i_rec:i_rec+10000,2,1:128) = NaN;
                L0.sar.wfm_ku_q(i_rec:i_rec+10000,64,1:128) = NaN;
                L0.sar.wfm_ku_i(i_rec:i_rec+10000,64,1:128) = NaN;
                L0.sar.avgwfm_ku(i_rec:i_rec+10000,1:128) = NaN;
                L0.sar.avgwfm_c(i_rec:i_rec+10000,1:128) = NaN;
            end
            
            L0.sar.val(i_rec)=1; % ground flag: use it for identifying the mode validity
            L0.sar.fine_datation(i_rec) = fread(fid,1,'uint32'); % bytes 19-22 -- SRAL Fine datation (LSB = 137.5 ns) 
            fseek(fid,1,'cof'); % byte 23 reserved
            L0.sar.mode_id(i_rec) = fread(fid,1,'uint8'); % byte 24 -- SAR CL or SAR OL
            fread(fid,1,'ubit5'); % byte 25, spare bits 0-4
            L0.sar.tracking_conf.closed_loop_gain(i_rec) = fread(fid,1,'ubit1'); % byte 25, bit 5 -- 0:nominal value, 1:nominal value with back-off
            L0.sar.tracking_conf.acquisition(i_rec) = fread(fid,1,'ubit1'); % byte 25, bit 6 -- 0:no acquisition, 1:acquisition
            L0.sar.tracking_conf.dem_eeprom_read_access(i_rec) = fread(fid,1,'ubit1'); % byte 25, bit 7 -- DEM EEPROM read access. 0:enabled, 1:disabled
            fread(fid,1,'ubit7'); % byte 26, spare bits 0-6
            L0.sar.altimeter_conf.weighting_function(i_rec) = fread(fid,1,'ubit1'); % byte 26, bit 7 -- 0:enabled, 1:disabled
            L0.sar.ho_nav(i_rec) = fread(fid,1,'uint32'); % byte 27-30 -- Distance command H0 computed with navigation information and DEM (LSB = 3.125/64 ns)
            L0.sar.cor2_nav(i_rec) = fread(fid,1,'int16'); % byte 31-32 -- Distance command COR2 computed with navigation information and DEM (signed integer) (LSB = 3.125/1024 ns)
            L0.sar.ho(i_rec) = fread(fid,1,'uint32'); % byte 33-36 -- Applied distance command H0 (LSB = 3.125/64 ns)
            L0.sar.cor2(i_rec) = fread(fid,1,'int16'); % byte 37-38 -- Applied distance command COR2 (signed integer)  (LSB = 3.125/1024 ns)
            fseek(fid,1,'cof'); % byte 39 reserved
            fread(fid,1,'ubit7'); % byte 40, spare bits 0-6
            L0.sar.loss_track(i_rec) = fread(fid,1,'ubit1'); % byte 40, bit 7 -- Loss of track criterion computed on the echo of the cycle (N-2) in OL mode (current cycle = N). 0:normal, 1:loss of track
            L0.sar.dist_error_ol(i_rec) = fread(fid,1,'int32'); % byte 41-44 -- Distance error computed on the echo of the cycle (N-2) in OL mode (current cycle = N)  (signed integer)  (LSB = 3.125/64 ns)
            L0.sar.att_code_ku(i_rec) = fread(fid,1,'uint8'); % byte 45 -- Ku band ATTCODE (LSB = 1 dB)
            L0.sar.att_code_c(i_rec) = fread(fid,1,'uint8'); % byte 46 -- C band ATTCODE (LSB = 1 dB)
            fread(fid,1,'ubit7'); % byte 47, spare bits 0-6
            L0.sar.nav_boletin_status(i_rec) = fread(fid,1,'ubit1'); % byte 47, bit 7 -- Navigation bulletin status. 0:Bulletin OK, 1:Bulletin KO
            fread(fid,1,'ubit6'); % byte 48, spare bits 0-5
            L0.sar.nav_boletin_source(i_rec) = fread(fid,1,'ubit2'); % byte 48, bits 6-7 -- Navigation bulletin source identifier. 0b01:GNSS, 0b10:DORIS
            L0.sar.nav_coarse_time(i_rec) = fread(fid,1,'uint32'); % bytes 49-52 -- Time of the current navigation bulletin: Navigation bulletin coarse time on 4 bytes (LSB = 1 s)
            L0.sar.nav_fine_time(i_rec) = fread(fid,1,'ubit24'); % bytes 53-55 -- Time of the current navigation bulletin: Navigation bulletin fine time on 3 bytes (LSB = 2^-24s)
            fseek(fid,1,'cof'); % byte 56 reserved
            L0.sar.therm.TH_KU_SSPA_THR1(i_rec) = fread(fid,1,'uint8'); % byte 57 -- Thermistances telemetry. THR1 : TH_KU_SSPA_THR1
            L0.sar.therm.TH_C_SSPA_THR2(i_rec) = fread(fid,1,'uint8'); % byte 58 -- Thermistances telemetry. THR2 : TH_C_SSPA_THR2
            L0.sar.therm.TH_KU_RX_THR3(i_rec) = fread(fid,1,'uint8'); % byte 59 -- Thermistances telemetry. THR3 : TH_KU_RX_THR3
            L0.sar.therm.TH_C_RX_THR4(i_rec) = fread(fid,1,'uint8'); % byte 60 -- Thermistances telemetry. THR4 : TH_C_RX_THR4
            L0.sar.therm.TH_FI_CNG_THR5(i_rec) = fread(fid,1,'uint8'); % byte 61 -- Thermistances telemetry. THR5 : TH_FI_CNG_THR5
            L0.sar.therm.TH_X16_THR6(i_rec) = fread(fid,1,'uint8'); % byte 62 -- Thermistances telemetry. THR6 : TH_X16_THR6
            L0.sar.therm.TH_DCDC_RFU_THR7(i_rec) = fread(fid,1,'uint8'); % byte 63 -- Thermistances telemetry. THR7 : TH_DCDC_RFU_THR7
            L0.sar.therm.ref_can1(i_rec) = fread(fid,1,'uint8'); % byte 64 -- Thermistances telemetry. THR8 : Reference CAN 1
            L0.sar.therm.TH_PLL_THR8(i_rec) = fread(fid,1,'uint8'); % byte 65 -- Thermistances telemetry. THR9 : TH_PLL_THR8
            L0.sar.therm.TH_CHIRP_THR9(i_rec) = fread(fid,1,'uint8'); % byte 66 -- Thermistances telemetry. THR10 : TH_CHIRP_THR9
            L0.sar.therm.TH_FI_THR10(i_rec) = fread(fid,1,'uint8'); % byte 67 -- Thermistances telemetry. THR11 : TH_FI_THR10
            L0.sar.therm.TH_ASIC_DFFT_THR11(i_rec) = fread(fid,1,'uint8'); % byte 68 -- Thermistances telemetry. THR12 : TH_ASIC_DFFT_THR11
            L0.sar.therm.TH_DCDC_THR12(i_rec) = fread(fid,1,'uint8'); % byte 69 -- Thermistances telemetry. THR13 : TH_DCDC_THR12
            fseek(fid,1,'cof'); % byte 70 -- Thermistances telemetry. THR14 : reserved
            fseek(fid,1,'cof'); % byte 71 -- Thermistances telemetry. THR15 : reserved
            L0.sar.therm.ref_can2(i_rec) = fread(fid,1,'uint8'); % byte 72 -- Thermistances telemetry. THR16 : Reference CAN 2
            fseek(fid,1,'cof'); % byte 73 reserved
            L0.sar.burst_number(i_rec) = fread(fid,1,'uint8'); % byte 74 -- Burst number (from 1 to 4)
            i_ku_echo=0; i_c_echo=0;
            for i_echo=1:66 % bytes 75-16970 -- 66 x 128 (I,Q) samples : 1 C pulse + 64 Ku pulses + 1 C pulse. Each (I,Q) sample on 2 bytes (I : LSB and Q : MSB) with LSB = 1
                if i_echo==1 || i_echo==66
                    i_c_echo=i_c_echo+1;
                    for i_samp=1:128
                        L0.sar.wfm_c_q(i_rec,i_c_echo,i_samp) = fread(fid,1,'int8');
                        L0.sar.wfm_c_i(i_rec,i_c_echo,i_samp) = fread(fid,1,'int8');
                    end
                else
                    i_ku_echo=i_ku_echo+1;
                    for i_samp=1:128
                        L0.sar.wfm_ku_q(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                        L0.sar.wfm_ku_i(i_rec,i_ku_echo,i_samp) = fread(fid,1,'int8');
                    end
                end
            end
            for i_samp=1:128 % bytes 16791-17226 -- 128 Ku samples (averaged echo) : each sample on 2 bytes (LSB = 1) (Available only in the TM corresponding to the fourth burst. For the TM corresponding to the 3 first bursts, these fields are filled with the value 0.)
                L0.sar.avgwfm_ku(i_rec,i_samp) = fread(fid,1,'uint16');
            end
            for i_samp=1:128 % bytes 17227-17482 -- 128 C samples (averaged echo) : each sample on 2 bytes (LSB = 1) (Available only in the TM corresponding to the fourth burst. For the TM corresponding to the 3 first bursts, these fields are filled with the value 0. )
                L0.sar.avgwfm_c(i_rec,i_samp) = fread(fid,1,'uint16');
            end
            
        otherwise
            disp(['PCAT = ' num2str(L0.packet_header.packet_category(i_rec)) '. No Valid Packet Category. ABORT reading routine']);
            return
            
    end
    
    %% read ISP packet error control (last 2 bytes)
    L0.packet_error_control(i_rec) = fread(fid,1,'uint16'); % Packet Error Control (last field)
    
end

fclose (fid);

save([fileparts(filename) '.mat'],'L0');


toc
end