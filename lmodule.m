function L=lmodule( omega_a, n) ;

% module of elliptic filter characteristic function
% Lj. Milic, M. Lutovac
% "Design of elliptic IIR filters with
% a reduced number of shift-and-add
% operations in multipliers"
% 3rd IEEE Conf. ICECS96, Rodos, Greece

phi=asin(1/omega_a);
if phi<=pi/4,
    t=1/2*(1-sqrt(cos(phi)) ) /(1+sqrt(cos(phi) ) ) ;
    q=t+2*t^5+15*t^9+150*t^13;
else
    t=1/2*(1-sqrt(sin(phi) ) ) /(1+sqrt(sin(phi)));
    q=t+2*t^5+15*t^9+150*t^13;
    q=exp(pi^2/log(q));
end

L=sqrt(1/16/q^n-1);
