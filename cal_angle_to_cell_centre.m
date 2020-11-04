function [users] = cal_angle_to_cell_centre(users,settings, beam_Centrals)
%UNTITLED2 
%   
angle_to_Beam_Centrals = cell(settings.users_per_Beam,settings.num_of_Beams);
beam_Central_to_Original = sqrt(beam_Centrals(:,1) .^2 + beam_Centrals(:,2) .^2);

phi = beam_Central_to_Original/settings.earth_Radius;
beam_Central_to_Satellite = sqrt(settings.earth_Radius ^ 2 + ...
                          (settings.earth_Radius + settings.satellite_Height)^2 ... 
                           - 2 * settings.earth_Radius * (settings.earth_Radius + settings.satellite_Height) * cos(phi));

for k = 1:settings.num_of_Beams
    position = users.positions(:,2*k-1 : 2*k);
    for q = 1:settings.users_per_Beam
        user_Coordinate = position(q,:);
        user_to_Beam_Central = beam_Centrals - user_Coordinate;
        arc_Distance_to_Beam_Centrals = sqrt(user_to_Beam_Central(:,1) .^2 + user_to_Beam_Central(:,2) .^2);
        theta = arc_Distance_to_Beam_Centrals/settings.earth_Radius;
        chord_Distance_to_Beam_Centrals = 2 * settings.earth_Radius * sin(theta/2);
        angle_to_Beam_Centrals{q,k} = acosd((beam_Central_to_Satellite .^2 + users.distance_to_Satellite(q,k) .^2 ...
            -chord_Distance_to_Beam_Centrals.^2) ./(2 * users.distance_to_Satellite(q,k) * beam_Central_to_Satellite));
        
    end
end
users.angle_to_Beam_Centrals = angle_to_Beam_Centrals;
    
end

