function [precoding_Matrix , SINR] = optimization_changed(channel_Matrix,settings)
%optimization
% 
[p,W_k] = precoding(channel_Matrix, settings);
M = 2;


for m = 1:M
    V = [];
      for k = 1:settings.num_of_Beams
        [U,S] = schur(W_k(:,:,k));

        v_l = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
        w_k = U * sqrtm(S) * v_l;
        v_k = w_k / norm(w_k);
        V = [V,v_k];  
      end
    SINR = [];
    cvx_begin quiet gp
       variable P(1 ,settings.num_of_Beams);
       variable t;
       R = [];
       maximize t;       
        subject to
         t >= 1;
         for k = 1:settings.num_of_Beams
            0 <= P(k); 
            %t * settings.SINR_Threshold < min(SINR);
         end
         for n = 1:settings.num_of_Antenna
             P_tot = 0
            % P_Antenna = [];
             for k = 1:settings.num_of_Beams
                %P_Antenna = [P_Antenna , P(k) * square(norm(V(n,k))) ];
                P_tot = P_tot +  P(k) * square(norm(V(n,k)));
             end
             %sum(P_Antenna)  <= settings.power_per_Antenna;
             P_tot  <= settings.power_per_Antenna;
         end
%          W_sum = 0;
%          for k = 1:settings.num_of_Antenna
%             W_sum = W_sum + P(k) * real(V(:,k) * V(:,k)'); 
%          end
         
         
         for k = 1:settings.num_of_Beams
           Beam_SINR = [];
           for q = 1:settings.users_per_Beam                
               h = channel_Matrix(:,(k-1)*settings.users_per_Beam + q);
               signal_Power = P(k) * square(norm(h' * V(:,k)));
               interference_Power = [];
               for l = 1:settings.num_of_Beams
                  if l ~= k
                     interference_Power =  [interference_Power,P(l) * square(norm(h' * V(:,l)))]; 
                  end
               end
               signal_Power >= t *(sum(interference_Power) + settings.noise_Power) * settings.SINR_Threshold;                      
           end
       end 
         


    cvx_end
    W = V .* sqrt(P);
end
P
t
precoding_Matrix = W


end

