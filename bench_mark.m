function [result] = bench_mark(channel_Matrix, settings, num_of_Users)
%bench_mark
%

R_sum_temp = 0;
SINR_temp = [];
P_n_temp = [];
W_temp = [];
data_Rate_temp = [];


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
minimize sum_Power;
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
         W_l = W_l + W(:,:,l);
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
       delta = normrnd(0 , 1 ,  7 , 1);
       R = sqrt(chi2inv(1 - settings.outage_Probability,7));
       M = [Q+t(k,q)*eye(settings.num_of_Antenna) , r;r',s - t(k,q)*( R^2 )];
       M >= 0;
    end
    
    
end
cvx_end;

V = [];
for k = 1:settings.num_of_Beams
  [U,S] = schur(W_opt(:,:,k));
  v_s = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
  w_k = U * sqrtm(S) * v_s;
  v_k = w_k / norm(w_k);
  V = [V,v_k];
end

SINR = [];
precoding_Matrix = [];
cvx_begin quiet
variable P(1,settings.num_of_Beams);
W = P .* V;
minimize trace(W*W')
subject to


cvx_end

result = W;
end

