function k  = bndrychng(min,max,C,deps); % boundaries for power change

C=(floor(log2(C)));

if (min==-Inf) | (imag(min)~=0),
    if (max==-Inf) | (imag(max)~=0),
        k=-Inf;
    else
        k=[-Inf C-deps:max];
    end
else
    if (max==-Inf) | (imag (max) ~=0),
        k=-Inf;
    else
        k=min:max;
    end
end
