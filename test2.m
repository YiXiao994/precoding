sum = 0;
for i = 1:10000
    
theta_error = normrnd(0 , 30);
v = exp(1j * 2 * pi * theta_error ./ (360));
sum = sum + v;

end

avg = sum/10000
exp(-1 * (  2 *pi* 30 / 360)^2)