function diffbeta=beta_ellip(fa,n,Ap,Aa,F3dB,betamaxmax);

%  compute  difference  of beta
%  used  for solving betamax<fa)<0

fp=atan(tan(pi*F3dB)^2/tan(pi*fa))/pi;
omega_a=tan(pi*fa)/tan(pi*fp);
L=lmodule(omega_a,n);
ap=10*log10(1+1/L);
aa=10*log10(1+L);
% L,n,ap,aa,fp
[z,p,c]=ellip(n,ap,aa,2*fp);
betamax=max(abs(p).^2) ;

diffbeta=betamax-betamaxmax;
