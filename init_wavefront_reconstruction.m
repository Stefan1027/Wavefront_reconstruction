% @file    init_wavefront_reconstruction.m
% @author  Stefan-Tudor Ilie
% @date    03/09/18
% @version 1.2 - Release 


% Clearing Variables and Environment

close all;
clear all;
clc;

% 0 - Pixel geometry, 1 - Micron geometry

geometry_select = 1;

% Checks if all the files required are in Matlab directory and loads them, else errors
% Pixel geometry, deviation oriented on centre of lenslet.

if geometry_select == 1
    if exist('WFS_deviation_file_x_processed.csv', 'file')
        WFSdeviationfilex = load('WFS_deviation_file_x_processed.csv');
    elseif exist('WFS_deviation_file_x.csv', 'file')
        WFSdeviationfilex = load('WFS_deviation_file_x.csv');
    else
      % File does not exist.
      warningMessage = sprintf('Warning: files do not exist in current location:\n%s\n%s', 'WFS_deviation_file_x.csv','WFS_deviation_file_x_processed.csv');
      uiwait(msgbox(warningMessage));
    end

    if exist('WFS_deviation_file_y_processed.csv', 'file')
        WFSdeviationfiley = load('WFS_deviation_file_y_processed.csv');
    elseif exist('WFS_deviation_file_y.csv', 'file')
        WFSdeviationfiley = load('WFS_deviation_file_y.csv');
    else
      % File does not exist.
      warningMessage = sprintf('Warning: files do not exist in current location:\n%s\n%s', 'WFS_deviation_file_y.csv', 'WFS_deviation_file_y_processed.csv');
      uiwait(msgbox(warningMessage));
    end

    if exist('WFS_field_intensity_processed.csv', 'file')
        WFSfieldintensity = load('WFS_field_intensity_processed.csv');
    elseif exist('WFS_field_intensity.csv', 'file')
        WFSfieldintensity = load('WFS_field_intensity.csv');
    else
      % File does not exist.
      warningMessage = sprintf('Warning: files do not exist in current location:\n%s\n%s', 'WFS_field_intensity.csv', 'WFS_field_intensity_processed.csv');
      uiwait(msgbox(warningMessage));
    end
else
    
% Checks if all the files required are in Matlab directory and loads them, else errors
% Micron geometry, deviation oriented in between two adjecent lenslets.

    if exist('WFS_deviation_file_x_avg.csv', 'file')
        WFSdeviationfilex = load('WFS_deviation_file_x.csv');
    elseif exist('WFS_deviation_file_x_avg_processed.csv', 'file')
        WFSdeviationfilex = load('WFS_deviation_file_x_avg_processed.csv');
    else
      % File does not exist.
      warningMessage = sprintf('Warning: files do not exist in current location:\n%s\n%s', 'WFS_deviation_file_x_avg.csv','WFS_deviation_file_x_avg_processed.csv');
      uiwait(msgbox(warningMessage));
    end

    if exist('WFS_deviation_file_y_avg.csv', 'file')
        WFSdeviationfiley = load('WFS_deviation_file_y_avg.csv');
    elseif exist('WFS_deviation_file_y_avg_processed.csv', 'file')
        WFSdeviationfiley = load('WFS_deviation_file_y_avg_processed.csv');
    else
      % File does not exist.
      warningMessage = sprintf('Warning: files do not exist in current location:\n%s\n%s', 'WFS_deviation_file_y_avg.csv', 'WFS_deviation_file_y_avg_processed.csv');
      uiwait(msgbox(warningMessage));
    end

    if exist('WFS_field_intensity_avg.csv', 'file')
        WFSfieldintensity = load('WFS_field_intensity_avg.csv');
    elseif exist('WFS_field_intensity_avg_processed.csv', 'file')
        WFSfieldintensity = load('WFS_field_intensity_avg_processed.csv');
    else
      % File does not exist.
      warningMessage = sprintf('Warning: files do not exist in current location:\n%s\n%s', 'WFS_field_intensity_avg.csv', 'WFS_field_intensity_avg_processed.csv');
      uiwait(msgbox(warningMessage));
    end

end

% Reconstructed zernike wavefront from WFS loading

if exist('WFS_wavefront_zernike_fit.csv', 'file')
    WFSzernikefit = load('WFS_wavefront_zernike_fit.csv');
end

% Algorithm supports only 2^N by 2^N input datasets, N - integer
% Checking if dataset corresponds in size, else cuts it to nearest
% available unit, by cutting the first a-31 or a-15 elements (rows and
% columns), where a is size of the input matrix

[a,b]=size(WFSdeviationfilex);

if a > 32
    WFSdeviationfilex = WFSdeviationfilex(((a-31):end),((a-31):end));
    WFSdeviationfiley = WFSdeviationfiley(((a-31):end),((a-31):end));
    WFSfieldintensity = WFSfieldintensity(((a-31):end),((a-31):end));
    WFSzernikefit = WFSzernikefit(((a-31):end),((a-31):end));
elseif (a > 17 && a < 32)
    WFSdeviationfilex = WFSdeviationfilex(((a-15):end),((a-15):end));
    WFSdeviationfiley = WFSdeviationfiley(((a-15):end),((a-15):end));
    WFSfieldintensity = WFSfieldintensity(((a-15):end),((a-15):end));
    WFSzernikefit = WFSzernikefit(((a-15):end),((a-15):end));
elseif a < 16
    warningMessage = sprintf('Warning: Input array is not big enough, make sure the input has more than 16x16 elements');
    uiwait(msgbox(warningMessage));
else
    sprintf('The input datasets are either 32 by 32 or 16 by 16');
end

% input_type = 0 if random generated data is required, else turbulance
% based data set is used for the initial phase and ampltiude
% Note: Turbulance phase/amplitude dataset needs to generated using Waveprop

input_type = 0;

if input_type == 0
    phase2pi=2*pi+2.*rand(450,450);
    amplitude2pi=( 3*10^(-3))+2.*rand(450,450);
else
    load('amp_n_phase_0.567_1.mat');
    phase2pi=phase2pi+pi;
end

if exist('WFS_focal.csv', 'file')
    f=load('WFS_focal.csv');
end

SNR0=60; %Signal to Noise Ratio
NN=2^6; %Number of spots per whole pupil of the sensor
NoiseFlag=0;  %1 adds noise. 0 does NOT add noise

% sensor_type = 0 if Model of Wavefront is intended to be used else 1, if
% the sensor data will be used in calculation of the x,y slopes.

sensor_type = 1;

if sensor_type == 0
    [Fx,Fy,Magnitudes,SNR] = slope_WgtAvg(phase2pi,amplitude2pi,NN,SNR0,NoiseFlag);
else
    [NN_sensor,n]=size(WFSdeviationfilex);
    [Fx,Fy,Magnitudes,SNR,NN] = wavefront_sensor(WFSdeviationfilex,WFSdeviationfiley,NN_sensor,SNR0,NoiseFlag,WFSfieldintensity);
end

VLQx=1./Magnitudes;
VLQy=VLQx;
if sensor_type == 0
    dx=size(phase2pi,1)/NN;
else
    dx=size(WFSdeviationfilex,1)/NN;
end

%WFSzernikefit=flip(WFSzernikefit,1);
%WFSzernikefit=flip(WFSzernikefit,2);  %matrix flip to correct error in GUI output
for i = 1:32
    for j = 1:32
        WFSzernikefit(i,j)= -WFSzernikefit(i,j);
    end
end
[E,phase_CER,VLQ]=NVWCER4HSD(Fx,Fy,VLQx,VLQy);      %NVWCER Phase Reconstruction
phase_zonal=zonal_2(Fx,Fy);                        %Zonal Phase Reconstruction
[Zernike_Phase,Terms]=ZernikeModal(Fx*dx,Fy*dx);       %Modal Phase Reconstruction

% Streth ratio calculation and figure of comparison will only work when the
% wavefront model is in place. As the streth ratio weighted measurement of
% the difference between original phase and reconstructed phase at each point.

% The phase data presented is wrapped, meaning it is on the domain of zero to 2?.
% The interest of reconstruction is generating an unwrapped phase, meaning it exist on the domain of  -? to ?

%Max
[t,f]=maxk(phase_CER,1,1);
[t1,f1]=maxk(phase_CER,1,2);
max_horiz = mean(maxk(t,5));
max_vert = mean(maxk(t1,5));

%Min
[t2,f2]=mink(phase_CER,1,1);
[t3,f4]=mink(phase_CER,1,2);
min_horiz = mean(mink(t2,5));
min_vert = mean(mink(t3,5));

str = sprintf('Maximum average H = %f, V = %f, Minimum average H = %f, V= %f', max_horiz, max_vert, min_horiz, min_vert);
dim = [0 0 0.1 0.2];


if sensor_type == 0
    % Strehl Ratio Calculations;

    ims=size(phase2pi,1);
    dd=ims/NN;  %Spacing for each subsample distance 

    undersamp_southwell=amplitude2pi(dd/2:dd:ims-dd/2,dd/2:dd:ims-dd/2).*exp(i*phase2pi(dd/2:dd:ims-dd/2,dd/2:dd:ims-dd/2));
    undersamp_fried=amplitude2pi(1:dd:ims,1:dd:ims).*exp(i*phase2pi(1:dd:ims,1:dd:ims));

    reff_sw=(sum(sum(abs(undersamp_southwell))))^2;
    reff_fr=(sum(sum(abs(undersamp_fried))))^2;

    comp_zonal=exp(-1i*phase_zonal);
    comp_modal=exp(-1i*Zernike_Phase); 

    comp_CER_unwrap=exp(-1i*phase_CER(1:64,1:64));
    comp_CER_wrap=conj(E(1:64,1:64))./abs(E(1:64,1:64));

    strehl_zonal=(abs(sum(sum(undersamp_southwell.*comp_zonal))))^2/reff_sw; 
    strehl_zernike=(abs(sum(sum(undersamp_southwell.*comp_modal))))^2/reff_sw;   

    strehl_CER_unwrap=(abs(sum(sum(undersamp_fried.*comp_CER_unwrap))))^2/reff_fr; 
    strehl_CER_wrap=(abs(sum(sum(undersamp_fried.*comp_CER_wrap))))^2/reff_fr; 


    % Wrapped and Unwrapped plots of all the reconstructors when using WFS Model 



    figure;
    surf(phase_zonal);title(['Zonal Unwrapped ' 'Strehl=' num2str(strehl_zonal) ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap gray; 
    saveas(gcf,'Zonal_Unwrapped.fig','fig'); 

    figure;
    surf(Zernike_Phase);title(['Zernike Fit Unwrapped ' 'Strehl=' num2str(strehl_zernike) ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap gray; 
    saveas(gcf,'Zernike_Fit_Unwrapped.fig','fig'); 

    figure;
    imagesc(mod(phase_CER,2*pi));title(['NVWCER Wrapped ' 'Strehl=' num2str(strehl_CER_unwrap) ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap gray; 
    saveas(gcf,'NVWCER_Wrapped.fig','fig'); 

    figure;
    imagesc(mod(phase_zonal,2*pi));title(['Zonal Wrapped ' 'Strehl=' num2str(strehl_zonal) ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap gray; 
    saveas(gcf,'Zonal_Wrapped.fig','fig'); 

    figure;
    imagesc(mod(Zernike_Phase,2*pi));title(['Zernike Fit Wrapped ' 'Strehl=' num2str(strehl_zernike) ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap gray; 
    saveas(gcf,'Zernike_Fit_Wrapped.fig','fig');
    
    figure;
    surf(phase_CER);title(['NVWCER Unwrapped ' 'Strehl=' num2str(strehl_CER_unwrap) ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap gray;
    saveas(gcf,'NVWCER_Unwrapped.fig','fig'); 
    
else
    
    % Wrapped and Unwrapped plots of all the reconstructors when using data from WFS

    figure;
    surf(phase_zonal); rotate3d on; title(['Zonal Unwrapped ' ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap jet; 
    saveas(gcf,'Zonal_Unwrapped_sensor.fig','fig'); 

    figure;
    surf(Zernike_Phase);rotate3d on;title(['Zernike Fit Unwrapped ' ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap jet; 
    saveas(gcf,'Zernike_Fit_Unwrapped.fig','fig'); 

    figure;
    imagesc(mod(phase_CER,2*pi));rotate3d on;title(['NVWCER Wrapped ' ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap jet; 
    saveas(gcf,'NVWCER_Wrapped_sensor.fig','fig'); 

    figure;
    imagesc(mod(phase_zonal,2*pi));rotate3d on;title(['Zonal Wrapped ' ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap jet; 
    saveas(gcf,'Zonal_Wrapped_sensor.fig','fig'); 

    figure;
    imagesc(mod(Zernike_Phase,2*pi));rotate3d on;title(['Zernike Fit Wrapped ' ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap jet; 
    saveas(gcf,'Zernike_Fit_Wrapped.fig','fig');
    
    figure;
    surf(WFSzernikefit);rotate3d on;title(['Zernike Fit based on WFS ' 'N=' num2str(a)]);colorbar; colormap jet; 
    saveas(gcf,'Zernike_Fit_WFS.fig','fig'); 
    
    figure;
    surf(phase_CER); rotate3d on; title(['NVWCER Unwrapped ' ' SNR=' num2str(SNR0) 'N=' num2str(NN)]);colorbar; colormap jet;annotation('textbox',dim,'String',str,'FitBoxToText','on');
    saveas(gcf,'NVWCER_Unwrapped_sensor.fig','fig'); 
    
end
