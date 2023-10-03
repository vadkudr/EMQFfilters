function [b,y]=next_p2(a,x);

% next_p2.m
%   compute next p2 number
% Input:
% a - power-of-two representation
% x - indexes of non-zero elements
% Output:
% b - power-of-two representation
% y - indexes of non-zero elements

%% for debug
% a=[1 0 0 1 0 0 -1];
% x=[1 4 7];

n=length(a);
nx=length(x);

number=3.^[nx-2:-1:0]*(a(x(2:end))+1)';
if number~=3^(nx-1)-1,
	number=number+1;

	% 3-adic representation
	c=[];
	for i=1:nx-1,
		temp=rem(number,3);
		c=[temp c];
		number=(number-temp)/3;
	end
        
    c=c-1;
    a(x(2:end))=c;
else
	% new position
	if ~(all(x==[n-nx+1:n])) ,
		i=nx;
		while i>0,
			if x(i)~=n-nx+i,
				x(i)=x(i)+1;
				x(i+1:end)=[x(i)+1:x(i)+nx-i];
				break;
			else
				i=i-1;
			end
		end
		a=zeros(1,n);
		a(x)=[1 -ones(1,nx-1)];
	end
end

b=a;
y=x;