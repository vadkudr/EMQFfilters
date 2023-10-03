function b=p2p2(p,n);

% % debug
% p=4;        % number of non-zero bits
% n=16;        % number of significant bits


r=[];

% cycle for all  value combinations  (in ternary system)
for i=0:3^p-1,
    a=i;
	u=[];
	for j=1:p,
		u=[u rem(a,3)];
		a=floor(a/3);
    end
    u=u-1;
	a=zeros(n*ones(1,p));
    for j=1:p,
		b=u(j)*2.^(0:n-1)';                                 % values
		b=repmat(b,[1 n*ones(1,p-1)]);            % extend to other dimensions
		ind=1:p; ind(1)=j; ind(j)=1;
		b=permute(b,ind);                                         % exchange dimensions
		a=a+b;
    end
	a=reshape(a,1,prod(size(a)));                     % make one-dim array
	r=[r a];
end

%  postprocessing
r=r(find(r>=0));
r=sort(r);                                                                 %  sort it in scending order
t1=find(diff(r)~=0);
b=r([1 t1+1]);