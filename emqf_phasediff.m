function [psi,f]=emqf_phasediff (alpha,beta, alpha1,freq,Fs);

% phase difference function for EMQF
beta_a=beta(2:2:end); beta_b=beta(1:2:end);
p0=[1 alpha*(1+beta_a(1)) beta_a(1)];
for i=2:length(beta_a),
    p0=conv(p0, [1 alpha *(1+beta_a(i)) beta_a(i) ]);
end
p0=conv(p0,[1 alpha1]); d0=p0;
p0=fliplr(p0);

p1=[1 alpha*(1+beta_b(1)) beta_b(1)];
for i=2:length(beta_b),
    p1=conv(p1, [1 alpha*(1+beta_b(i)) beta_b(i)]);
end
d1=p1; p1=fliplr(p1);

[h0, f ]=freqz (p0, d0, freq,Fs);
[h1, f ]=freqz (p1, d1, freq,Fs);

psi=(phase(h0)-phase(h1))/2;
