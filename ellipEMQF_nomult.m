function [alpha,beta,alpha1]=ellipEMQF_nomult(ord, Fp, Fa, Ap, Aa) ;
% function [a0_q,b0_q,a1_q,b1_q]=ellipEMQF_nomult(ord, Fp, Fa);
% function [a0_q,b0_q,a1_q,b1_q]=ellipEMQF_nomult(Fp, Fa, Ap, Aa);
% Multiplierless elliptic lowpass IIR filter design
% Lj. D. Milic, M. D. Lutovac,
% Design of multiplierless elliptic IIR filters with a small quantization error,
% IEEE Trans. Signal Processing, Vol. 47, no. 2, Feb. 1999, pp. 469-479,
% Fp, Ap - passband cutoff and ripple
% Fa, Aa - attenuation band cutoff and ripple

% restrictions
% 1. ord -> odd number
% for debug
% ord=9;
% Fp=0.25-0.04;
% Fa=0.25+0.04;
% Ap=0.2;     % dB
% Aa=30;     % dB
% Fp=0.13 5;
% Fa=0.2;
% Ap=0.2;     % dB
% Aa=30;     % dB
% Fp=0.115-0.0113; % Fa=0.115+0.0113;
% Ap=5;     % dB % Aa=32;     % dB
%
% ord=7;%ord=ellipord(2*Fp,2*Fa,Ap,Aa);      % B1
% ord=ceil(ord)+(1-rem(ord,2));      % B2  % close greater odd number
%%% Choice 'alpha' %%%
disp('Calculating ''alpha''...');
% find quantized 3dB frequency
F3dB=(Fp+Fa)/2;
% set of frequencies corresponding to alpha=2^i+2^k
alpha=[-1+1./2.^(8:-1:2) -1./2.^(1:8) 0 1./2.^(8:-1:1) 1-1./2.^(2:8)];
f=1/pi*atan(sqrt((1+alpha)./(1-alpha)));
[temp,ind]=min(abs(f-F3dB));
F3dB_q=f(ind) ;                                            % new optimized F3dB corresponding to alpha(ind)
alpha=alpha(ind);
%%% Choice of 'betamax' %%%

disp('Calculating maximal ''beta''...');
fa1=Fa;fp1=atan(tan(pi*F3dB_q)^2/tan(pi*fa1))/pi;   % B3
omega_a=tan(pi*fa1)/tan(pi*fp1);                                            % B4
L=lmodule(omega_a, ord) ;                             % B5
ap1=10*log10(1+1/L);                                   % B6
aa1=10*log10(1+L) ;                                       % B7
[z,p,c]=ellip(ord,ap1,aa1,2*fp1);   % B8
betamax1=max(abs(p).^2);                           % B9

L=10^(Aa/10)-1;                                              % B10
omega_a=fsolve (@lmodule_minus, 10, optimset('fsolve' ), ord,L) ;    % B11
fp2=atan(tan(pi*F3dB_q)/sqrt(omega_a))/pi;  %B12
ap2=10*log10(1+1/L);                                   % B13
aa2=Aa;                                                               % B14
[z,p,c]=ellip(ord,ap2,aa2,2*fp2);   % B15
betamax2=max(abs(p).^2);                           % B16
betamax=pof2appr(betamax1, betamax2);
% betamax=1-2^(-2)+2^(-5)+2^(-7);   % for debug

% find 'alpha1' nad 'beta's and filter response
fa=fzero(@beta_ellip,Fa,optimset ('fzero'), ord,Ap,Aa,F3dB_q, betamax);    % B20
fp=atan(tan(pi*F3dB_q)^2/tan(pi*fa))/pi;                                                                % B21
omega_a=tan(pi*fa)/tan(pi*fp);                                                                                     % B22
L=lmodule(omega_a,ord);                                                                                                   % B23
ap=10*log10(1+1/L);                                                                                                           %B24
aa=10*log10(1+L);                                                                                                                % B25
[z,p,c]=ellip(ord,ap,aa,2*fp);                                                                                     % B26
%%% prepare list of coeffs for 'a' and fb' branch %%%
beta=p(find(imag(p)>0)) ;                          % remove real pole and pole with negative imag part
beta=sort(abs(beta).^2)';                          % all betas
alpha1=-(1-tan(pi*F3dB_q))/(1+tan(pi*F3dB_q));
%%% beta's quantization %%%
disp('Calculating ''beta''s...');
Da=pi/2-acos(10^(-Aa/20));                                      % C4
for i=1:length(beta)-1,%length(beta)-1:-1:1
%   calculate  response   in  stopband
%         [num,den]=emqfco2tf(alpha,beta,alpha1);                                    %   transfer  function
%         [h,f]=freqz(num,den, [fa:(0.5-fa)/511:0.5],1) ;   h=abs (h) ;         %  magnitude response
%         psi=acos(h);                                                         %  phase-difference response   (4)(9)
    [psi,f]=emqf_phasediff(alpha,beta,alpha1, [Fa: (0.5-Fa)/511:0.5] , 1) ;
    psi=-psi;   %   if mean(psi)<0  psi=-psi;end

    % C1 find extremums
    t1=sign(diff(psi)); t1(find(t1==0))=1;      % differentiate and compare with zero
    t2=abs(diff(t1));   t2=find(t2>0)+1;       % indices of extremum points
    t3=diff(t1); t3=t3(t2-1); % max or min flag
    if t3(1)>0, t3=[-2 t3]; else t3=[2 t3]; end

    extr_freq=[Fa f(t2)];                                                                                                                                                      % C1
    extr_psi=psi([1 t2]);                                                                                                                                                     % C2
    delta_psi=extr_psi-pi/2*sign(mean(extr_psi));                                                                                                                                % C3
    if any(abs(delta_psi)>Da),
        error('Required stopband ripple conditions are too strong. Relax initial conditions.');
    end
                                            % C4
    dd=delta_psi+(sign(t3) ) *Da;        % C5
    diff_psi=(4*(alpha+cos(2*pi*extr_freq)).*sin(2*pi*extr_freq))./((1-beta(i))^2+4*beta(i)*cos(2*pi*extr_freq).^2+alpha*(1+beta(i))^2*(alpha+2*cos(2*pi*extr_freq)))/2*(1-2*rem(i,2));    % C6
    dbeta=dd./diff_psi;                                                                                                                                                                   % C7

    % find allowable boundaries
    ind=find(dbeta<=0); beta_bndry_min=max(beta(i)+dbeta(ind));
    ind=find(dbeta>=0); beta_bndry_max=min(beta(i)+dbeta(ind));

    temp=pof2appr(beta_bndry_min(1), beta_bndry_max(1)) ;
    if isempty(temp) ,
        error(['Parameter "beta(' num2str(i) ')" can not be represented even as a sum of four power-of-two terms. Relax initial conditions.']) ;
    end
    beta (i)=temp(1) ;
end

%%%% alphal quantization %%%
disp('Calculating ''alpha1''...');
[psi,f]=emqf_phasediff(alpha,beta,alpha1, [Fa: (0.5-Fa)/511:0.5],1);
psi=-psi; % if mean(psi)<0 psi=-psi;end

%  C10   find  extremums
t1=sign(diff(psi));   t1(find(t1==0))=1;             % differentiate and compare with zero
t2=abs(diff(t1));       t2=find(t2>0)+1;                 %   indices   of  extremum points
t3=diff(t1);   t3=t3(t2-1);   % max  or min  flag
if  t3(1)>0,   t3=[-2  t3]; else   t3=[2   t3]; end

extr_freq=[Fa f(t2)];                                                                                                             % C10
extr_psi=psi([1 t2]);                                                                                                             % C11
delta_psi=extr_psi-pi/2*sign(mean(extr_psi));                                                                                             % C12
                                                    % C13
dd=delta_psi+(sign(t3))*Da;        % C14
diff_psi=2*sin(2*pi*extr_freq)./(1+alpha1^2+2*alpha1*cos(2*pi*extr_freq))/2;    % C15
dalpha1=dd./diff_psi;                                                                                                                       % C16                                                                   ~

% find allowable boundaries
ind=find(dalpha1<=0);
alpha1_bndry_min=max(alpha1+dalpha1(ind));
ind=find(dalpha1>=0);
alpha1_bndry_max=min(alpha1+dalpha1(ind));

temp=pof2appr (alpha1_bndry_min, alpha1_bndry_max);
if isempty(temp),
    error (['Parameter ''alpha1'' can not be represented even as a sum of four power-of-two terms. Relax initial conditions.']);
end
alpha1=temp(1);

% alpha,   beta,   alpha1
% dec2bin(round (abs (alpha) *2^16) ) , dec2bin (round (beta*2^16)) , dec2bin(round (abs (alpha1) *2^16))
