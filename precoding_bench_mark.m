function [W] = precoding_bench_mark(channel_Matrix, settings, num_of_Users)
% precoding_bench_mark
%
C = (2 * pi * (settings.phase_Error_Standard_Deviation / 360))^2 * eye(settings.num_of_Antenna);

cvx_begin quiet sdp;

variable W_opt(settings.num_of_Antenna, settings.num_of_Antenna,settings.num_of_Beams) complex semidefinite;
variable t(settings.num_of_Beams, num_of_Users);
expression f_X(settings.num_of_Antenna,settings.num_of_Antenna);
expression g_Y(settings.num_of_Antenna,1);
power_Matrix = 0;
for k = 1:settings.num_of_Beams
   power_Matrix = power_Matrix + W_opt(:,:,k);
end

sum_Power = trace(power_Matrix);
minimize real(sum_Power);
subject to

for n = 1:settings.num_of_Antenna
   real (power_Matrix(n,n)) <= settings.power_per_Antenna;
   real (power_Matrix(n,n)) >= 0;
end

for k = 1:settings.num_of_Beams
    W_k = W_opt(:,:,k);
    W_l = 0;
    for l = 1:settings.num_of_Beams
      if l ~= k
         W_l = W_l + W_opt(:,:,l);
      end
    end
    for q = 1:num_of_Users
       h = channel_Matrix(:,(k-1)*num_of_Users + q);
       Z = diag(h')*(W_k - settings.SINR_Threshold(k) * W_l)*diag(h);
       X = real(Z);
       Y = imag(Z);
 %      f_X = zeros(settings.num_of_Antenna, settings.num_of_Antenna);
 %      g_Y = zeros(settings.num_of_Antenna, 1);
       for m = 1:settings.num_of_Antenna
          g_Y(m) = 2 * sum(Y(m,:));
          for n = 1:settings.num_of_Antenna
              if m == n
                  f_X(m,n) = 2 * X(m,n) - sum(X(m,:));
              else
                  f_X(m,n) = X(m,n);
              end
              
          end           
       end
       Q = sqrtm(C) * f_X * sqrtm(C);
       r = 0.5 * sqrtm(C) * g_Y ;
       s = sum(Z(:)) - settings.SINR_Threshold(k);
       R = sqrt(chi2inv(1 - settings.outage_Probability,7));
       M = [Q+t(k,q)*eye(settings.num_of_Antenna) , r;r',s - t(k,q)*( R^2 )];
       M >= 0;
    end
    
    
end
cvx_end;

W = W_opt;
end

