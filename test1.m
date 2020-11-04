W_k = W(:,:,1);
h = channel_Matrix(:,1);
W_l = 0;
for l = 1:7
    if l ~= 1
    W_l = W(:,:,l) + W_l;
    end
end
count = 0;
for i = 1:1000
theta_error = normrnd(0 , settings.phase_Error_Standard_Deviation ,  settings.num_of_Antenna , 1)
h_real = h .* exp(1j * 2 * pi * theta_error ./ (360))

SINR_real = real(( h_real' * W_k * h_real)/(h_real' * W_l * h_real + 1))
SINR = real(trace(diag(h) * A * diag(h') * W_k )) / ( real(trace( diag(h) * A * diag(h') * W_l ) ) + 1 )

if SINR_real >= settings.SINR_Threshold(1)
   count = count + 1; 
end

end

prob = count / 1000