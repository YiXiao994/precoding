function [channel_Matrix] = init_channel_matrix(users , settings)
%init_channel_Matrix
%   
channel_Matrix = [];
path_Loss_dB = 20 * log10((settings.light_Speed)./(4 * pi * settings.frequency * users.distance_to_Satellite));
large_Scale_Fading_dB = path_Loss_dB + settings.receive_Gain_over_Temperature + ...
    - 10 * log10(settings.Boltzmann_Constant * settings.user_Link_Bandwidth);

for k = 1:settings.num_of_Beams
    rain_attenuation_dB = normrnd(settings.rain_Fading_Mean , settings.rain_Fading_Standard_Deviation);
    for q = 1:settings.users_per_Beam
       phase_Vector = exp(1j * 2 * pi * (repelem(unifrnd(0,360) , settings.num_of_Antenna) / 360 ))';
       u = 2.07123 * sind(users.angle_to_Beam_Centrals{q,k}) / sind(settings.three_dB_Angle);
       beam_Gain = 10 * log10(settings.max_Transmit_Gain) + 10 * log10( (besselj(1,u)./(2 * u) ...
                                                                + 36 * besselj(3,u)./(u .^ 3)).^2 );
       channel_Amplitude_dB = (rain_attenuation_dB + large_Scale_Fading_dB(q,k) + beam_Gain)/2;
       channel_Amplitude = 10 .^ (channel_Amplitude_dB / 10);
       channel_Vector = channel_Amplitude .* phase_Vector;
       channel_Matrix = [channel_Matrix , channel_Vector];
    end
end


end

