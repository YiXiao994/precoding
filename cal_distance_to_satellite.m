function [users] = cal_distance_to_satellite(users,settings)
%cal_distance_to_satellite
%   
distance_to_Original = users.positions(:,1:2:2*settings.num_of_Beams - 1) .^ 2 + users.positions(:,2:2:2*settings.num_of_Beams);
distance_to_Original = sqrt(distance_to_Original);
phi = distance_to_Original/settings.earth_Radius;
users.distance_to_Satellite = sqrt(settings.earth_Radius ^ 2 + ...
                          (settings.earth_Radius + settings.satellite_Height)^2 ... 
                           - 2 * settings.earth_Radius * (settings.earth_Radius + settings.satellite_Height) * cos(phi));
end

