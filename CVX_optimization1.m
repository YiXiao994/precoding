function [P_cluster_next , SINR_next] = CVX_optimization1(V,SINR ,settings , num_of_Users, channel_Matrix)
%CVX_optimization
%
  cvx_begin quiet
     variable P(settings.num_of_Beams , 1)
     expression S(settings.num_of_Beams, num_of_Users)
     expression I(settings.num_of_Beams, num_of_Users)
     obj = 0;
     for k = 1:settings.num_of_Beams
          obj_temp = [];
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
            %S(k,q) = real( trace(diag(h) * A * diag(h') * W_kk ) );
            %I(k,q) = real( trace(diag(h) * A * diag(h') * W_ll) ) + 1;
            S(k,q) = real(h' * W_kk * h);
            I(k,q) = real(h' * W_ll * h) + 1;
            obj_temp = [obj_temp , S(k,q) - settings.SINR_Threshold(k) * SINR(k,q)* I(k,q) ];
            %obj = obj + (S(k,q) - settings.SINR_Threshold(k) * SINR(k,q)* I(k,q));
          end
          obj = obj + num_of_Users * min(obj_temp);
     end
        
        maximize obj
        
        subject to
        power_matrix = 0;
        for k = 1:settings.num_of_Beams
          v_k = V(:,k);
          P(k) >= 0
          W_kk =  P(k) * v_k * v_k';
          power_matrix = power_matrix + W_kk;
           for q = 1:num_of_Users
              S(k,q) >= settings.SINR_Threshold(k) * I(k,q); 
           end
        end
        
        for n = 1:settings.num_of_Antenna
           real(power_matrix(n,n)) >= 0;
           real(power_matrix(n,n)) <= settings.power_per_Antenna;
           
        end
        
%         for k = 1:settings.num_of_Beams
%           W_kk = P(k) * v_k * v_k';
%           W_ll = 0;
%           for l = 1:settings.num_of_Beams
%             if l ~= k
%                W_ll = W_ll + P(l) * v_l * v_l';
%                %Z = W(:,:,k) - settings.outage_Ratio * settings.SINR_Threshold(k) * W(:,:,l); 
%             end
%           end
%           Z = W_kk - settings.SINR_Threshold(k) * W_ll;
%           for q = 1:num_of_Users 
%              h = channel_Matrix(:,(k-1)*num_of_Users + q);
%              C = diag(h) * Z * diag(h');
%              mu = real(trace(C * A));
%            %sig_square = (vec(C.').')
%              a = settings.SINR_Threshold(k) * settings.noise_Power;
%              b = sqrt(2) * erfinv(1 - 2 * settings.outage_Probability);
%              norm(sqrtm(G) * vec(C')) <=  (1/b) * ( sqrt(b^2 + 1) * mu -  a/(sqrt(b^2 + 1)) );
%            %norm(sqrtm(G) * vec(C')) <= (a - mu)/(b)
%             % S(k,q) >= settings.SINR_Threshold(k) * I(k,q);
%              %S(k,q) >= 0.01 * I(k,q);
%           end
%           end
     cvx_end
       SINR_next = S ./ I;
       P_cluster_next = P;
