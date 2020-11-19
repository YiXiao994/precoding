function [result] = optimization_changed(channel_Matrix,settings,num_of_Users)
%optimization
% 
[p,W_opt] = precoding_outage(channel_Matrix, settings , num_of_Users);
M = 1;

A = [];
a_m = [1];
G2 = [];
g2 = [exp( -(2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
  
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


for m = 1:M
    V = [];
    P_cluster = zeros(settings.num_of_Beams , 1);
    P_cluster_next = zeros(settings.num_of_Beams , 1);
    SINR = zeros(settings.num_of_Beams , num_of_Users);
   % SINR1 = zeros(settings.num_of_Beams , num_of_Users);
      for k = 1:settings.num_of_Beams
        [U,S] = schur(W_opt(:,:,k));

        v_l = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
        %w_k = U * sqrtm(S) * v_l;
        w_k = mvnrnd(zeros(settings.num_of_Antenna,1),W_opt(:,:,k)).';
        v_k = w_k / norm(w_k);
        V = [V,v_k];
      end
      
      W_init = V .* (sqrt( 0.3 * settings.power_per_Antenna)); 
      for k = 1:settings.num_of_Beams
        P_cluster(k) = norm( W_init(:,k))^2;
        v_k = V(:,k);
        W_k =  P_cluster(k) * v_k * v_k';
        W_l = 0;
        W_opt_l = 0;
        for l=1:settings.num_of_Beams
            if l ~=k
                v_l = V(:,l);
                W_l = W_l + P_cluster(l) * v_l * v_l';
                W_opt_l = W_opt_l +  W_opt(:,:,l) ;
            end
        end
        for q = 1:num_of_Users
          h = channel_Matrix(:,(k-1)*num_of_Users + q);
          SINR(k,q) = real((h' * W_k * h)) / (real( h' * W_l * h) + 1);
        end
      end
    count = 0;
    P_cluster
    SINR
    while true
      
       [P_cluster_next , SINR] = CVX_optimization(V,SINR,settings , num_of_Users,channel_Matrix);
       count = count + 1;
       SINR
       if norm(P_cluster_next - P_cluster) <= 0.01
           P_cluster = P_cluster_next;
           count
           break;
       end
       
       if count >= 50
          disp('Exceed the maximum iterative times.')
          P_cluster = P_cluster_next;
          count
          break;
       end
       
       P_cluster = P_cluster_next;
    end
    
    W = sqrt(P_cluster') .* V;
    cvx_begin
    variable r(settings.num_of_Beams);
    obj = 0;
    for k = 1:settings.num_of_Beams
        obj = obj + r(k) * settings.SINR_Threshold(k);
    end
    maximize obj
    
    subject to
    
    for k = 1:settings.num_of_Beams
        W_k = W(:,k) * W(:,k)';
        W_l = 0;
        for l = 1:settings.num_of_Beams
            if l ~= k
                W_l = W_l + W(:,l) * W(:,l)';
            end
        end
        
        for q = 1:num_of_Users
           h = channel_Matrix(:,(k-1)*num_of_Users + q);
           Z_k = W_k - settings.SINR_Threshold(k) * r(k) * W_l;
           C_k = diag(h) * Z_k * diag(h');
           mu = real(trace(C_k * A));
           a = r(k) * settings.SINR_Threshold(k) * settings.noise_Power;
           b = sqrt(2) * erfinv( 1 - 2 * settings.outage_Probability);
           norm(sqrtm(G) * vec(C_k')) <= (1/b) * ( sqrt(b^2 + 1) * mu -  a/(sqrt(b^2 + 1)) )
        end
    end
    
    cvx_end
    
    r;

       
    end
end