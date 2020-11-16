w_k = result1(1).precoding_Matrix(:,2);
W_k = w_k * w_k';
h = channel_Matrix(:,11);
W_l = 0;
for l = 1:7
    if l ~= 2
        
    w_l = result1(1).precoding_Matrix(:,l);
    W_l = w_l * w_l' + W_l;
    end
end
count = 0;
for i = 1:10000
theta_error = normrnd(0 , settings.phase_Error_Standard_Deviation ,  settings.num_of_Antenna , 1)
h_real = h .* exp(1j * 2 * pi * theta_error ./ (360))

SINR_real = real(( h_real' * W_k * h_real))/(real(h_real' * W_l * h_real) + 1)
%SINR = real(trace(diag(h) * A * diag(h') * W_k )) / ( real(trace( diag(h) * A * diag(h') * W_l ) ) + 1 )

if SINR_real >= settings.SINR_Threshold(2)
   count = count + 1; 
end

end

prob = count / 10000