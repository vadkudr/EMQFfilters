function [num,den]=emqf_co2tf(alpha,beta,alpha1);

% transfer function of EMQF filter by given coefficients

beta_a=beta(2:2:end); beta_b=beta(1:2:end);

p0=[1 alpha* (1+beta_a(1)) beta_a(1)];
for i=2:length(beta_a),
    p0=conv(p0, [1 alpha* (1+beta_a (i) ) beta_a(i)]);
end

p0=conv(p0,[1 alpha1]);
d0=p0;
p0=fliplr(p0);

p1=[1 alpha* (1+beta_b(1)) beta_b(1)];
for i=2:length(beta_b),
    p1=conv(p1,[1 alpha*(1+beta_b(i)) beta_b(i)]);
end

d1=p1;
p1=fliplr(p1);

den=conv(d0,d1);
num=(conv(p0,d1)+conv(p1,d0))/2;
