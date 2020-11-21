function [settings] = init_settings()
%init_settings 
%   
settings.num_of_Beams = 7;
settings.num_of_Antenna = 7;
settings.selected_Users_per_Beam = 4;
settings.users_per_Beam = settings.selected_Users_per_Beam * (settings.selected_Users_per_Beam + 1) / 2;
settings.frequency = 28e9;  %% Hz
settings.Boltzmann_Constant = 1.38e-23;
settings.user_Link_Bandwidth = 500e6;  %% Hz
settings.max_Transmit_Gain = 26;  %% In dB
settings.three_dB_Angle = 0.4; %% Degree
settings.receive_Gain_over_Temperature = 34; %% dB/K
settings.rain_Fading_Mean = -2.6; %% dB
settings.rain_Fading_Standard_Deviation = 1.63; %% dB
settings.phase_Error_Standard_Deviation = 5;
settings.power_per_Antenna = 250 ;%% Watt
settings.cell_Radius = 250e3; %% meter
settings.satellite_Height = 35800e3; %% meter
settings.light_Speed = 2.99792458e8; %% m/s
settings.earth_Radius = 6371e3; %% meter
settings.noise_Power = 10 ^(0);  %% Normalized noise power
settings.SINR_Threshold = [0.46,0.46,0.46,0.46,0.46,0.46,0.46]-0.35; %%SINR requirement for each user group
settings.outage_Probability = 0.1;
settings.retry_count = 10000;

end

