function [p0,d0,p1,d1]=apellip_du(n,fp,fs);
% polyphase allpass filter realization of EMQF elliptic filter
% design is oriented to realization of filter as a sum of two allpasses
% H(z)=z^(-l)*A0(z^2)+Al(z^2)
%
% Important properties:
% - has constant alpha's in two-mult, two-delay realization
%
% Input:                                                                                                                                                                       9
% n - order of filter (only odd values accepted)
% fp - cutoff passband frequency
% fs - cutoff stopband frequency
%                                                                                                                                                                                                                                          Jjj
% Output:
% p0,d0 - numerator and denominator of first allpass
% p1,d1 - numerator and denominator of second allpass
%
% Reference:
% Lj. D. Milic, M. D. Lutovac,
% Efficient Algorithm for the Design of High-Speed Elliptic IIR Filters,
% Int. Journal of Electronics and Communications, AEU, Vol. 57, No. 4, 2003, pp. 255-262.

% debug
% n=7;                                                                                                                                                                                !
%   fp=0.2;
% fs=0.3;

% check for odd
if (rem(n,2)~=1) & (n<3)
	error('filter order should be odd greater than or equal to 3');
end

f3=atan(sqrt(tan(pi*fs)*tan(pi*fp)))/pi;

% alpha coeff
alpha=-cos(2*pi*f3);
% real pole
alpha1=-(1-tan(pi*f3))/(1+tan(pi*f3));
% beta coeffs
ksi=tan(pi*fs)/tan(pi*fp);
l=(n-1)/2;
x=ellipj(((2*(1:l)-1)/n+1)*ellipke(1/ksi^2),1/ksi^2);
lambda=2*sqrt(ksi)*tan(pi*fp)/(1+ksi*tan(pi*fp)^2);
beta=(ksi+x.^2-lambda*sqrt((1-x.^2).*(ksi^2-x.^2)))./(ksi+x.^2+lambda*sqrt((1-x.^2).*(ksi^2-x.^2)));
% distributing beta
beta=sort(beta);
p0=[alpha1 1];

for i=2:2:length(beta),
    p0=conv(p0,[beta(i) alpha*(1+beta(i)) 1]); end

d0=fliplr(p0);
p1=[1];

for i=1:2:length(beta),
    p1=conv(p1,[beta(i) alpha*(1+beta(i)) 1]); end

d1=fliplr(p1);