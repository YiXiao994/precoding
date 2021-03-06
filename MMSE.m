function [result] = MMSE(channel_Matrix, settings, num_of_Users)
%MMSE
% 
H_i = zeros(settings.num_of_Beams, settings.num_of_Antenna ,num_of_Users);
for q = 1:num_of_Users
   for k = 1:settings.num_of_Beams

     H_i(k,:,q) = channel_Matrix(:,(k-1)*num_of_Users + q)';
       
   end
   
end
H_hat = mean(H_i , 3);
W_mmse = inv(H_hat' * H_hat + (1 / settings.power_per_Antenna) * eye(settings.num_of_Antenna) ) * H_hat';
beta_mmse = sqrt(settings.power_per_Antenna * settings.num_of_Antenna / trace(W_mmse' * W_mmse));
W = beta_mmse * W_mmse;

SINR_real = zeros( num_of_Users , settings.num_of_Beams);
data_rate_real = zeros(1 , settings.num_of_Beams);
utility_real = zeros(num_of_Users, settings.num_of_Beams);
sum_weighted_SINR = 0;
outage_count = zeros(num_of_Users, settings.num_of_Beams);
for count = 1:settings.retry_count
  for k = 1:settings.num_of_Beams
    W_k = W(:,k) * W(:,k)';
    W_l = 0;
    for l = 1:settings.num_of_Beams
        if l ~= k
            W_l = W_l + W(:,l) * W(:,l)';
        end
    end
        
    for q = 1:num_of_Users
       h_estimated = channel_Matrix(:,(k-1)*num_of_Users + q);
       phase_error = normrnd(0 , settings.phase_Error_Standard_Deviation, settings.num_of_Antenna, 1);
       h = h_estimated .* exp(1j * 2 * pi * (phase_error/360) );
       signal_power = real(h' * W_k * h);
       interference_power = real(h' * W_l * h);
       SINR_temp = (signal_power) / (interference_power + 1);
       SINR_real(q,k) = SINR_real(q,k) +SINR_temp;  
       utility_real(q,k) = utility_real(q,k) + settings.SINR_Threshold(k) * log2(SINR_temp);
       sum_weighted_SINR = sum_weighted_SINR + SINR_temp / settings.SINR_Threshold(k);
       if SINR_temp < settings.SINR_Threshold(k)
          outage_count(q,k) = outage_count(q,k) + 1; 
       end
       
    end    
  end
     
end
P_antenna = diag(W * W');
P_cluster = diag(W' * W);
      SINR_real = SINR_real / settings.retry_count;
      utility_real = utility_real / settings.retry_count;
      sum_weighted_SINR = sum_weighted_SINR / settings.retry_count;
      outage_prob = outage_count / settings.retry_count;

      
      individual_utility = sum(utility_real(:)) / (settings.num_of_Beams * num_of_Users);
      if num_of_Users == 1
        min_SINR = SINR_real;
      else
        min_SINR = min(SINR_real);
      end
      P_antenna = real(diag(W_k + W_l)');
      jain_Index = square( sum(SINR_real(:) ) ) / (settings.num_of_Beams * num_of_Users * sum(square(SINR_real(:))));
      individual_rate_real = log2(1 + min_SINR) ;
      sum_rate_real = num_of_Users * sum(individual_rate_real); 
      %individual_rate_real = sum_rate_real / (settings.num_of_Beams * num_of_Users);
    result.precoding_Matrix = W;
    result.beam_power = P_cluster;
%    result.obj = obj;
    result.outage_Prob = outage_prob;
    result.SINR_optimized = SINR_real;
    result.min_SINR = min_SINR;
    result.sum_Rate = sum_rate_real;
    result.data_Rate = individual_rate_real;
    result.utility = individual_utility;
    result.jain_Index = jain_Index;
    result.antenna_Power = P_antenna;
    result.sum_weighted_SINR = sum_weighted_SINR;
end