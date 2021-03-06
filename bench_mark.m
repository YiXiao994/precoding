function [result] = bench_mark(channel_Matrix, settings, num_of_Users)
%bench_mark
%

R_sum_temp = 0;
SINR_temp = [];
P_n_temp = [];
W_temp = [];
data_Rate_temp = [];


C = (2 * pi * (settings.phase_Error_Standard_Deviation / 360))^2 * eye(settings.num_of_Antenna);
W_opt = precoding_bench_mark(channel_Matrix, settings, num_of_Users);
V = [];
P_antenna = [];
for k = 1:settings.num_of_Beams
  [U,S] = schur(W_opt(:,:,k));
  v_s = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
  w_k = U * sqrtm(S) * v_s;
  %w_k = mvnrnd(zeros(settings.num_of_Antenna,1), W_opt(:,:,1)).';
  v_k = w_k / norm(w_k);
  V = [V,v_k];
end

SINR = [];
precoding_Matrix = [];
cvx_begin sdp
variable P(1,settings.num_of_Beams);
variable t(settings.num_of_Beams, num_of_Users);
expression f_X(settings.num_of_Antenna,settings.num_of_Antenna);
expression g_Y(settings.num_of_Antenna,1);
power_Matrix = 0;
for k = 1:settings.num_of_Beams
   v_k = V(:,k);
   power_Matrix = power_Matrix + P(k) * v_k * v_k';

end
sum_power = trace(power_Matrix);
minimize real(sum_power);
subject to
for n = 1:settings.num_of_Antenna
   P_antenna = [P_antenna,power_Matrix(n,n)];
   real (power_Matrix(n,n)) <= settings.power_per_Antenna;
   real (power_Matrix(n,n)) >= 0;
end

for k = 1:settings.num_of_Beams
    v_k = V(:,k);
    W_k = P(k)* v_k * v_k';
    W_l = 0;
    for l = 1:settings.num_of_Beams
      if l ~= k
         v_l = V(:,l);
         W_l = W_l + P(l) * v_l * v_l';
      end
    end
    for q = 1:num_of_Users
       
       h = channel_Matrix(:,(k-1)*num_of_Users + q);
       %SINR = (h' * W_k * h)/(h' * W_l * h + 1);
       %real(h' * W_k * h)>= settings.SINR_Threshold(k) * (real(h' * W_l * h) + 1);
       %(h' * (W_k - settings.SINR_Threshold(k) * W_l )* h) >= settings.SINR_Threshold(k);
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
P_antenna
P
result = sqrt(P) .* V ;
end

