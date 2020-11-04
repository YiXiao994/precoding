function [precoding_Matrix , SINR] = optimization(channel_Matrix,settings)
%optimization
% 
[p,W_k] = precoding(channel_Matrix, settings);
V = [];
for k = 1:settings.num_of_Beams
    [U,S] = schur(W_k(:,:,k));
    v_l = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
    w_k = U * sqrtm(S) * v_l;
    v_k = w_k / norm(w_k);
    V = [V,v_k];
end
M = 20;
for m = 1:M
    W = [];
    SINR = []
    cvx_begin quiet
       variable P(1,settings.num_of_Beams);
       R_sum = 0;
       power_Matrix = 0;
       for k = 1:settings.num_of_Beams
           power_Matrix = power_Matrix + P(k) * V(:,k) * V(:,k)';
           Beam_SINR = [];
           for q = 1:settings.users_per_Beam                
               h = channel_Matrix(:,(k-1)*settings.users_per_Beam + q);
               R = h * h';
               signal_Power = real(trace(R * (P(k) * V(:,k) * V(:,k)')));
               interference_Power = 0;
               for l = 1:settings.num_of_Beams
                  if l ~= k
                     interference_Power = interference_Power + real(trace(R *(P(l) * V(:,l) * V(:,l)'))); 
                  end
               end
               %user_SINR = signal_Power / (interference_Power + settings.noise_Power);
               user_SINR_inversed = (interference_Power + settings.noise_Power)/ signal_Power ;
               Beam_SINR = [Beam_SINR , user_SINR_inversed];                           
           end
          % min_SINR = min(Beam_SINR);
           min_SINR = max(Beam_SINR);
           R_sum = R_sum + settings.users_per_Beam * log2(1+1/min_SINR);
           SINR = [SINR , Beam_SINR];
       end                     
       maximize R_sum;
       

        subject to
         for k = 1:settings.num_of_Beams
            P(k) >= 0; 
         end
         for n = 1:settings.num_of_Antenna
            power_Matrix(n,n)  <= settings.power_per_Antenna;      
         end
         
         for k = 1:settings.num_of_Beams
           for q = 1:settings.users_per_Beam 
               h = channel_Matrix(:,(k-1)*settings.users_per_Beam + q);
               R = h * h';
               signal_Power = real(trace(R * (P(k) * V(:,k) * V(:,k)')));
               interference_Power = 0;
               for l = 1:settings.num_of_Beams
                  if l ~= k
                     interference_Power = interference_Power + real(trace(R *(P(l) * V(:,l) * V(:,l)'))); 
                  end
               end
               signal_Power - settings.SINR_Threshold * (interference_Power + settings.noise_Power) >= 0;
           end

    end
    cvx_end
    W = V .* sqrt(P);
end
r
P
precoding_Matrix = W
end

