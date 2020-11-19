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
    SINR1 = zeros(settings.num_of_Beams , num_of_Users);
      for k = 1:settings.num_of_Beams
        [U,S] = schur(W_opt(:,:,k));

        v_l = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
        %w_k = U * sqrtm(S) * v_l;
        w_k = mvnrnd(zeros(settings.num_of_Antenna,1),W_opt(:,:,k)).';
        v_k = w_k / norm(w_k);
        V = [V,v_k];
        P_cluster(k) = max( eig(W_opt(:,:,k) ) );

      end
      
      for k = 1:settings.num_of_Beams
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
          SINR1(k,q) = real(h' * W_opt(:,:,k)* h) / (1 + real(h' * W_opt_l * h));
          SINR(k,q) = real((h' * W_k * h)) / (real( h' * W_l * h) + 1);
        end
      end
    count = 0;
    while true
      
       [P_cluster_next , SINR] = CVX_optimization(V,SINR1,settings , num_of_Users,channel_Matrix);
       count = count + 1;
       if norm(P_cluster_next - P_cluster) <= 0.01
           result = P_cluster_next;
           break;
       end
       
       if count >= 20
          disp('Exceed the maximum iterative times.')
          result = P_cluster_next;
          break;
       end
       
       P_cluster = P_cluster_next;
    end
       %result = P_cluster;
    end
end