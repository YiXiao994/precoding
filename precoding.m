function [ratio,W_opt] = precoding(channel_Matrix , settings, num_of_Users)
%Precoding
% 
cvx_begin quiet sdp
  variable W(settings.num_of_Antenna, settings.num_of_Antenna,settings.num_of_Beams) complex semidefinite;
  power_Matrix = 0;
  for k = 1:settings.num_of_Beams
     power_Matrix = power_Matrix + W(:,:,k);
  end
  variable p;
  minimize p;
  subject to
    p <= 1;
    0 <= p;
    for n = 1:settings.num_of_Antenna
     real (power_Matrix(n,n)) <= p * settings.power_per_Antenna;
     real (power_Matrix(n,n)) >= 0;
    end
    for k = 1:settings.num_of_Beams
       for q = 1:num_of_Users 
           h = channel_Matrix(:,(k-1)*num_of_Users + q);
           R = h * h';
           signal_Power = trace(R*W(:,:,k));
           interference_Power = 0;
           for l = 1:settings.num_of_Beams
              if l ~= k
                 interference_Power = interference_Power + trace(R * W(:,:,l)); 
              end
           end
           signal_Power - settings.SINR_Threshold(k) * (interference_Power + settings.noise_Power) >= 0;
       end

    end
    
cvx_end
ratio = p;
W_opt = W;
end

