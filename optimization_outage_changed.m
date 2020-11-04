function [result] = optimization_outage_changed(channel_Matrix,settings,num_of_Users)
%optimization
% 
[p,W_k] = precoding_outage(channel_Matrix, settings , num_of_Users);


M = 1;
R_sum_temp = 0;
SINR_temp = [];
P_n_temp = [];
W_temp = [];
data_Rate_temp = [];

A = [];
a = [1];
G2 = [];
g2 = [exp( -(2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
  
for n = 2:settings.num_of_Antenna
  a = [a , exp( -(2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
  g2 = [g2 , exp( -2 * (2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
end
  
for n = 1:settings.num_of_Antenna
  A = [A ; a];
  a = circshift(a ,1);
  G2 = [G2 ; g2];
  g2 = circshift(g2, 1);
end
G1 = A;
G = [];
for n1 = 1:settings.num_of_Antenna
  G_temp = [];
  for n2 = 1:settings.num_of_Antenna
     if n1 == n2
        G_temp = [G_temp , G1]; 
     else
        G_temp = [G_temp , G2];
        
     end
         
  end
  G = [G ; G_temp];
end

for m = 1:M
    V = [];

    for k = 1:settings.num_of_Beams
    [U,S] = schur(W_k(:,:,k));
    v_s = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
    w_k = U * sqrtm(S) * v_s;
    v_k = w_k / norm(w_k);
    V = [V,v_k];
    end
    V
    W = [];
    SINR = [];
    cvx_begin quiet gp
       variable P(settings.num_of_Beams,1);
       variable SINR_t(settings.num_of_Beams,1);
       min_SINR_inversed = [];

       for k = 1:settings.num_of_Beams
           Beam_SINR = [];
           for q = 1:num_of_Users               
               h = channel_Matrix(:,(k-1)*num_of_Users + q);
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
       R_sum_approximate = sum(num_of_Users * settings.SINR_Threshold .* log(SINR_t'));
       maximize R_sum_approximate;
              
       subject to
         for k = 1:settings.num_of_Beams
            0 <= P(k); 
           % P(k) <= 50;
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
           v_k = V(:,k);
           W_sum = 0;
           for q = 1:num_of_Users                
             h = channel_Matrix(:,(k-1)*num_of_Users+ q);                          
             m_k = SINR_t(k) * settings.noise_Power;
             n_k = sqrt(2) * erfinv(1 - 2 * settings.outage_Probability);             
             b_k = [];
             d_k = [];
             e_k = m_k / (sqrt(n_k^2 + 1));
             for l = 1:settings.num_of_Beams
                if l ~= k
                    v_l = V(:,l);
                    d = SINR_t(k)*sqrt(n_k^2 + 1) * trace(real(diag(h)*v_l*v_l'*diag(h'))*A) / (n_k);
                    d_k = [d_k ; d];
                    b = SINR_t(k)*norm(sqrtm(G) * vec((diag(h)*v_l*v_l'*diag(h'))'));
                    b_k = [b_k;b];
                else
                    d_k = [d_k ; 0];
                    c_k = sqrt(n_k^2 + 1) * trace(real(diag(h)*v_k*v_k'*diag(h'))*A) / (n_k);
                    b_k = [b_k;0];
                    a_k = norm(sqrtm(G) * vec(diag(h)*v_k*v_k'*diag(h')));
                end
             end
            c_k - a_k
   
             sum ((d_k - b_k).* P) + e_k <= (c_k - a_k) * P(k);
             %end
             %norm(sqrtm(G) * vec(C)) <= (1/b) * (sqrt(b^2 + 1) * real(mu) - a/(sqrt(b^2 + 1)) );                       
           end
           
         end 
    cvx_end
    min_SINR = 1./ min_SINR_inversed;
    min_SINR_inversed;
    
    
    W = V .* sqrt(P(k))
    sum_Rate = settings.num_of_Beams * sum(log2(1+min_SINR))
    data_Rate = log2(1+min_SINR);
   % sum_rate = settings.num_of_Beams * sum(1og2(1 + min_SINR));
        
           R_sum_temp = sum_Rate;
           SINR_temp = 1./SINR;
           P_n_temp = P_Antenna;
           W_temp = W;
           data_Rate_temp = data_Rate
     
     end
%
 result.beam_power = P';
 result.fairness = R_sum_approximate;
 result.sum_Rate = R_sum_temp;
 result.antenna_Power = P_n_temp;
 result.precoding_Matrix = W_temp;
 result.SINR_optimized = SINR_temp;
 result.data_Rate = data_Rate_temp;
end

