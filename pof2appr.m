function y=pof2appr(minn,maxx)

% load sets of pof2s
if exist('pof2nums_2.mat','file')
    load pof2nums_2    %  loads   'p2nums_2' array
else
    p2nums_2=p2p2(2,16);
    save pof2nums_2 p2nums_2
end

if exist('pof2nums_3.mat','file')
    load pof2nums_3    %  loads   'p2nums_3' array
else
    p2nums_3=p2p2(3,16);
    save pof2nums_3 p2nums_3
end

if exist('pof2nums_4.mat','file')
    load pof2nums_4    %  loads   'p2nums_4' array
else
    p2nums_4=p2p2(4,16);
    save pof2nums_4 p2nums_4
end

c=mean([minn,maxx]);

% if negative, make positive
c1=c*sign(c);
minn=minn*sign(c);
maxx=maxx*sign(c);

if minn>maxx, t=minn; minn=maxx; maxx=t; end

% minn,maxx
% 2222222
p=ceil(log2(maxx));
p2nums_2=p2nums_2/2^(15-p) ;
a=p2nums_2 (find( (p2nums_2>=minn) & (p2nums_2<=maxx) ) ) ;

if isempty(a)
	p2nums_3=p2nums_3/2^(15-p) ;
	a=p2nums_3 (find( (p2nums_3>=minn) & (p2nums_3<=maxx) ) ) ;
	if isempty(a)
		p2nums_4=p2nums_4/2^(15-p) ;
		a=p2nums_4 (find( (p2nums_4>=minn) & (p2nums_4<=maxx)) );
	end
end

[temp,ind]=min(abs(a-c1));

% restore sign
y=a(ind)*sign(c);