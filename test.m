settings.phase_Error_Standard_Deviation = 10;
A = [];
a = [1];
G2 = [];
g2 = [exp( -(2 * pi * settings.phase_Error_Standard_Deviation / 360)^2)];
  
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

W = result1(1).precoding_Matrix;
h = channel_Matrix(:,1);
W_k = W(:,1) * W(:,1)';
W_l = 0;
k = 1;
t = result1(1).t;
C = diag(h) * (W_k - settings.SINR_Threshold(k)*t *W_l) * diag(h');

for l = 1:settings.num_of_Beams
   if l ~= k
    W_l = W_l + W(:,l) * W(:,l)';
    
   end
end
power = real(trace(diag(h)*A*diag(h')*W_k))
interference = real(trace(diag(h)*A*diag(h')*W_l))
SINR = power/(interference + 1)
mu = real(trace(C * A))
v =  (norm(sqrtm(G)*vec(C')))^2 - mu^2 
prob = 0.5 + 0.5 * erfinv( (mu - t * settings.SINR_Threshold(k) )/(sqrt(2 * v)) )

