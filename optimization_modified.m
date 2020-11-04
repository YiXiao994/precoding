function [result] = optimization_modified(channel_Matrix,settings)
%optimization
% 
[p,W_k] = precoding(channel_Matrix, settings);


M = 1;
R_sum_temp = 0;
SINR_temp = [];
P_n_temp = [];
W_temp = [];
data_Rate_temp = [];
for m = 1:M
    V = [];

    for k = 1:settings.num_of_Beams
    [U,S] = schur(W_k(:,:,k));
    v_l = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
    w_k = U * sqrtm(S) * v_l;
    v_k = w_k / norm(w_k);
    V = [V,v_k];
    end
    W = [];
    SINR = [];
    cvx_begin quiet gp
       variable P(settings.num_of_Beams,1);
       min_SINR_inversed = [];
       for k = 1:settings.num_of_Beams
           Beam_SINR = [];
           for q = 1:settings.users_per_Beam                
               h = channel_Matrix(:,(k-1)*settings.users_per_Beam + q);
               signal_Power = P(k) * square(norm(h' * V(:,k)));
               interference_Power = [];
               for l = 1:settings.num_of_Beams
                  if l ~= k
                     interference_Power =  [interference_Power;P(l) * square(norm(h' * V(:,l)))]; 
                  end
               end
               %user_SINR_changed = (signal_Power + sum(interference_Power) + settings.noise_Power )/ (sum(interference_Power) + settings.noise_Power);
               user_SINR_inversed = (sum(interference_Power) + settings.noise_Power)/ signal_Power ;
               Beam_SINR = [Beam_SINR ; user_SINR_inversed];                           
           end
           %min_SINR = [min_SINR , min(Beam_SINR)];
           min_SINR_inversed = [min_SINR_inversed , max(Beam_SINR)];
          % min_SINR = max(Beam_SINR);
           
           SINR = [SINR , Beam_SINR];
       end
       R_sum_approximate = sum(settings.users_per_Beam * settings.SINR_Threshold .* log((1)./min_SINR_inversed));
       maximize R_sum_approximate;
              
       subject to
         for k = 1:settings.num_of_Beams
            0 <= P(k); 
            %t * settings.SINR_Threshold < min(SINR);
         end
         P_Antenna = [];
         for n = 1:settings.num_of_Antenna
             P_n = 0;
            % P_Antenna = [];
             for k = 1:settings.num_of_Beams
                %P_Antenna = [P_Antenna , P(k) * square(norm(V(n,k))) ];
                P_n = P_n +  P(k) * square(norm(V(n,k)));
             end
             %sum(P_Antenna)  <= settings.power_per_Antenna;
             P_n  <= settings.power_per_Antenna;
             P_Antenna = [P_Antenna, P_n];
         end
         
         for k = 1:settings.num_of_Beams
           %Beam_SINR = [];
           for q = 1:settings.users_per_Beam                
               h = channel_Matrix(:,(k-1)*settings.users_per_Beam + q);
               signal_Power = P(k) * square(norm(h' * V(:,k)));
               interference_Power = [];
               for l = 1:settings.num_of_Beams
                  if l ~= k
                     interference_Power =  [interference_Power,P(l) * square(norm(h' * V(:,l)))]; 
                  end
               end
               signal_Power >=  (sum(interference_Power) + settings.noise_Power) * settings.SINR_Threshold(k);                        
           end
         end 
    cvx_end
    min_SINR = 1./ min_SINR_inversed;
    min_SINR_inversed;
    
    
    W = V .* sqrt(P(k))
    sum_Rate = settings.num_of_Beams * sum(log2(1+min_SINR))
    data_Rate = log2(1+min_SINR);
   % sum_rate = settings.num_of_Beams * sum(1og2(1 + min_SINR));
         if sum_Rate > R_sum_temp
           R_sum_temp = sum_Rate;
           SINR_temp = 1./SINR;
           P_n_temp = P_Antenna;
           W_temp = W;
           data_Rate_temp = data_Rate
           stp = m;
        end
     end
%
 result.sum_Rate = R_sum_temp;
 result.antenna_Power = P_n_temp;
 result.precoding_Matrix = W_temp;
 result.SINR_optimized = SINR_temp;
 result.data_Rate = data_Rate_temp;
 stp
end

