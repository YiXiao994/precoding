function [ratio,precoding_Result] = precoding_outage(channel_Matrix , settings , num_of_Users)
%Precoding
% 

  A = [];
  a = [1];

  G2 = [];
  g2 = [exp( -( 2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
  
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

cvx_begin quiet sdp
  variable W(settings.num_of_Antenna, settings.num_of_Antenna,settings.num_of_Beams) complex semidefinite;
  power_Matrix = 0;
  for k = 1:settings.num_of_Beams
     power_Matrix = power_Matrix + W(:,:,k);
  end
  variable p;
  minimize p;
  subject to
    p <= 1;
    0 <= p;
    for n = 1:settings.num_of_Antenna
     real (power_Matrix(n,n)) <= p * settings.power_per_Antenna;
     real (power_Matrix(n,n)) >= 0;
    end
    for k = 1:settings.num_of_Beams
       W_k = W(:,:,k);
       W_l = 0;
       for l = 1:settings.num_of_Beams
          if l ~= k
             W_l = W_l + W(:,:,l);
             %Z = W(:,:,k) - settings.outage_Ratio * settings.SINR_Threshold(k) * W(:,:,l); 
          end
       end
       Z = W_k - settings.SINR_Threshold(k) * W_l;
       for q = 1:num_of_Users 
           h = channel_Matrix(:,(k-1)*num_of_Users + q);
           C = diag(h) * Z * diag(h');
           mu = real(trace(C * A));
           %sig_square = (vec(C.').')
           a = settings.SINR_Threshold(k) * settings.noise_Power;
           b = sqrt(2) * erfinv(1 - 2 * settings.outage_Probability);
           norm(sqrtm(G) * vec(C')) <= (1/b) * ( sqrt(b^2 + 1) * mu -  a/(sqrt(b^2 + 1)) );
           %norm(sqrtm(G) * vec(C')) <= (a - mu)/(b)
       end

    end
    
cvx_end
outage_Probability = zeros(settings.num_of_Beams , num_of_Users);

    
    for k = 1:settings.num_of_Beams
       for l = 1:settings.num_of_Beams
          if l ~= k
             Z = W(:,:,k) -  settings.SINR_Threshold(k) * W(:,:,l); 
          end
       end
       
       for q = 1:num_of_Users 
           h = channel_Matrix(:,(k-1)*num_of_Users + q);
           C = diag(h) * Z * diag(h');
           mu = real(trace(C * A));
           %sig_square = (vec(C.').' )
           v = sqrt((norm(sqrtm(G) * vec(C')))^2 - mu^2);
           a = settings.SINR_Threshold(k) * settings.noise_Power;
           b = sqrt(2) * erfinv( 1 - 2 * settings.outage_Probability);
           %norm(sqrtm(G) * vec(C')) <= (1/b) * (a/(sqrt(b^2 + 1)) - sqrt(b^2 + 1) * real(mu)  );
           %prob = 0.5+ 0.5 * erf(( settings.SINR_Threshold(k)-mu)/(sqrt(2) * norm(sqrtm(G) * vec(C'))) )()^2- mu^2)
           outage_Probability(k,q) = 0.5 + 0.5 * erf((mu -  settings.SINR_Threshold(k))/(sqrt(2)*v));

       end
    end
    
    SINR = zeros(num_of_Users , settings.num_of_Beams);
    min_SINR = zeros(1, settings.num_of_Beams);
    for k = 1:settings.num_of_Beams
       W_k = W(:,:,k);
       for q = 1:num_of_Users
           h = channel_Matrix(:,(k-1)*num_of_Users + q);
           signal_power = real(trace(diag(h) * A * diag(h') * W_k));
           interference_power = 0;
           
           for l = 1:settings.num_of_Beams
              if l ~= k
                 W_l = W(:,:,l);
                 interference_power = interference_power +  real(trace(diag(h) * A * diag(h') * W_l));
                  
              end
           end
           SINR(q,k) = signal_power / (interference_power + 1);
          
       end
       min_SINR(k) = min(SINR(:,k));
    end

ratio = p;
precoding_Result = W;
end

