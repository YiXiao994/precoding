function [result] = optimization_changed1(channel_Matrix,settings,num_of_Users)
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
        w_k = U * sqrtm(S) * v_l;
       % w_k = mvnrnd(zeros(settings.num_of_Antenna,1),W_opt(:,:,k)).';
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
    %count = 0;
%     if num_of_Users > 1
%         SINR = min(SINR')';
%         SINR = repmat(SINR , 1, num_of_Users);
%     end
    
    weighted_sum_SINR_iteration = sum(vec(settings.SINR_Threshold' .* SINR));
    for count = 1:10
      
       [P_cluster_next , SINR] = CVX_optimization1(V,SINR,settings , num_of_Users,channel_Matrix);
     %  count = count + 1;
%        if norm(P_cluster_next - P_cluster) <= 0.01
%            P_cluster = P_cluster_next;
%      %      count;
%            break;
%        end
       weighted_SINR = SINR ./ settings.SINR_Threshold'  ;
       weighted_sum = sum(vec(weighted_SINR));
       P_cluster = P_cluster_next;
       
       
       weighted_sum_SINR_iteration = [weighted_sum_SINR_iteration , weighted_sum ];
%        if num_of_Users > 1
%         SINR = min(SINR')';
%         SINR = repmat(SINR , 1, num_of_Users);
%        end
    end
    
    W = sqrt(P_cluster') .* V;
    cvx_begin quiet
    variable r(settings.num_of_Beams);
    obj = 0;
    for k = 1:settings.num_of_Beams
        obj = obj + r(k);% * settings.SINR_Threshold(k);
    end
    maximize obj
    
    subject to
    
    for k = 1:settings.num_of_Beams
        r(k) >= 1;
        W_k = W(:,k) * W(:,k)';
        W_l = 0;
        for l = 1:settings.num_of_Beams
            if l ~= k
                W_l = W_l + W(:,l) * W(:,l)';
            end
        end
        
        for q = 1:num_of_Users
           h = channel_Matrix(:,(k-1)*num_of_Users + q);
           Z_k = W_k - r(k) * settings.SINR_Threshold(k) * W_l;
           C_k = diag(h) * Z_k * diag(h');
           mu = real(trace(C_k * A));
           a = r(k)* settings.SINR_Threshold(k) * settings.noise_Power;
           b = sqrt(2) * erfinv( 1 - 2 * settings.outage_Probability);
           norm(sqrtm(G) * vec(C_k')) <= (1/b) * ( sqrt(b^2 + 1) * mu -  a/(sqrt(b^2 + 1)) )
        end
    end
    cvx_end
    outage_Threshold = r' .* settings.SINR_Threshold;
    
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
           if SINR_temp < outage_Threshold(k)
              outage_count(q,k) = outage_count(q,k) + 1; 
           end
        end    
      end
     
    end
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
    result.outage_Threshold = outage_Threshold;
    result.r = r;
    result.precoding_Matrix = W;
    result.beam_power = P_cluster;
    result.obj = obj;
    result.outage_Prob = outage_prob;
    result.SINR_optimized = SINR_real;
    result.min_SINR = min_SINR;
    result.sum_Rate = sum_rate_real;
    result.data_Rate = individual_rate_real;
    result.utility = individual_utility;
    result.jain_Index = jain_Index;
    result.antenna_Power = P_antenna;
    result.sum_weighted_SINR = sum_weighted_SINR;
    result.weighted_sum_SINR_iteration = weighted_sum_SINR_iteration;
    end
end