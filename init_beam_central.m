function [central_Positions] = init_beam_central(settings)
%init_beam_central 
%  
allowed_Cell_Num = [1,3,4,7];
if ismember(settings.num_of_Antenna , allowed_Cell_Num)
    if settings.num_of_Antenna == 7
        half_Central_Distance = sqrt(3) * settings.cell_Radius/2;
        central_Positions = zeros(7,2);
        central_Positions(1,:) = [0 ,0];
        central_Positions(2,:) = [-2 * half_Central_Distance,0];
        central_Positions(3,:) = [-half_Central_Distance, 1.5 * settings.cell_Radius];
        central_Positions(4,:) = [half_Central_Distance, 1.5 * settings.cell_Radius];
        central_Positions(5,:) = [2 * half_Central_Distance,0];
        central_Positions(6,:) = [half_Central_Distance, -1.5 * settings.cell_Radius];
        central_Positions(7,:) = [-half_Central_Distance, -1.5 * settings.cell_Radius];
    end
end

end

