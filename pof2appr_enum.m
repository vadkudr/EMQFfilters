% function y=pof2appr_enum(a);

% pof2appr_enum.m
%
% Calculate power-of-two approximation for number

% if exist('gfv'), close(gfv); clear(gfv); end

% filter prototype

n=16;   % significant bits
ONE=2^14;

na1=3;
% ia1=1;

a1ind=(1:na1);
a1bit=zeros(1,n); a1bit(a1ind(1))=1;     % initialize
a1bit(a1ind(2:end))=-1;

p2=fliplr(2.^[0:n-1]) ;
p2nums=[];

	% loop for 'a1'
	while(1)

		p2nums=[p2nums p2*a1bit'];

		% calculate next p2number 'a1'
		if (all(a1ind==[n-na1+1:n])) & (all(a1bit(a1ind)==1)), break; end
		[a1bit,a1ind]=next_p2(a1bit,a1ind);
		if a1ind(1)==n-na1+1, break; end;
	end