function [users] =init_user_positions(settings , beam_Centrals)
%init_user_positions 
%   

users.positions = [];
for k = 1:settings.num_of_Beams
    user_Distribution = settings.cell_Radius * unifrnd(0, 0.6, settings.users_per_Beam,1).* exp(1j*unifrnd(0,360,settings.users_per_Beam,1));
    user_Coordinates = beam_Centrals(k,:) + [real(user_Distribution), imag(user_Distribution)];
    users.positions = [users.positions , user_Coordinates];
    
end

end

