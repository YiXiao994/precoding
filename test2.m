%W_opt = precoding_bench_mark(H, settings, 4);
C = (2 * pi * (settings.phase_Error_Standard_Deviation / 360))^2 * eye(settings.num_of_Antenna);
V = [];
for k = 1:settings.num_of_Beams
  [U,S] = schur(W_opt(:,:,k));
  v_s = sqrt(1/2) * (randn(settings.num_of_Antenna , 1) + 1j * randn(settings.num_of_Antenna , 1));
  w_k = U * sqrtm(S) * v_s;
  v_k = w_k / norm(w_k);
  V = [V,v_k];
end

P = [100,150,200,110,210,200,170];
%W_1 = P(1) * V(:,1) * V(:,1)'
%eig(W_1)
f_X = zeros(settings.num_of_Antenna,settings.num_of_Antenna);
g_Y = zeros(settings.num_of_Antenna,1);

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
       %eig(Z)
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
       M = [Q+5*eye(settings.num_of_Antenna) , r;r',s - 5*( R^2 )];
       eig(M)
    end
    
    
end