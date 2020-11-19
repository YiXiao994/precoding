function [P_cluster_next , SINR_next] = CVX_optimization(V,SINR ,settings , num_of_Users, channel_Matrix)
%CVX_optimization
%  
  A = [];
  a_m = [1];

  G2 = [];
  g2 = [exp( -( 2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
  
  for n = 2:settings.num_of_Antenna
     a_m = [a_m , exp( -(2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
     g2 = [g2 , exp( -2 * (2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
  end
  
  for n = 1:settings.num_of_Antenna
     A = [A ; a_m];
     a_m = circshift(a_m ,1);
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
  cvx_begin
     variable P(settings.num_of_Beams , 1)
     expression S(settings.num_of_Beams, num_of_Users)
     expression I(settings.num_of_Beams, num_of_Users)
     obj = 0;
     for k = 1:settings.num_of_Beams
          v_k = V(:,k);
          W_kk =  P(k) * v_k * v_k';
          W_ll = 0;
          for l=1:settings.num_of_Beams
            if l ~=k
                v_l = V(:,l);
                W_ll = W_ll + P(l) * v_l * v_l';
            end
          end
          for q = 1:num_of_Users
            h = channel_Matrix(:,(k-1)*num_of_Users + q);
            S(k,q) = real( h' * W_kk * h);
            I(k,q) = real( h' * W_ll * h) + 1;
            obj = obj + (S(k,q) - SINR(k,q)* I(k,q));
          end
        end
        
        maximize obj
        
        subject to
        power_matrix = 0;
        for k = 1:settings.num_of_Beams
          v_k = V(:,k);
          P(k) >= 0
          W_kk =  P(k) * v_k * v_k';
          power_matrix = power_matrix + W_kk;
        end
        
        for n = 1:settings.num_of_Antenna
           real(power_matrix(n,n)) >= 0;
           real(power_matrix(n,n)) <= settings.power_per_Antenna;
           
        end
        
        for k = 1:settings.num_of_Beams
          W_kk = P(k) * v_k * v_k';
          W_ll = 0;
          for l = 1:settings.num_of_Beams
            if l ~= k
               W_ll = W_ll + P(l) * v_l * v_l';
               %Z = W(:,:,k) - settings.outage_Ratio * settings.SINR_Threshold(k) * W(:,:,l); 
            end
          end
          Z = W_kk - settings.SINR_Threshold(k) * W_ll;
          for q = 1:num_of_Users 
             h = channel_Matrix(:,(k-1)*num_of_Users + q);
             C = diag(h) * Z * diag(h');
             mu = real(trace(C * A));
           %sig_square = (vec(C.').')
             a = settings.SINR_Threshold(k) * settings.noise_Power;
             b = sqrt(2) * erfinv(1 - 2 * settings.outage_Probability);
             norm(sqrtm(G) * vec(C')) <=  (1/b) * ( sqrt(b^2 + 1) * mu -  a/(sqrt(b^2 + 1)) );
           %norm(sqrtm(G) * vec(C')) <= (a - mu)/(b)
            % S(k,q) >= settings.SINR_Threshold(k) * I(k,q);
             %S(k,q) >= 0.01 * I(k,q);
          end
          end
     cvx_end
       SINR_next = S ./ I;
       P_cluster_next = P;


