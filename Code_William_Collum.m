%Input island number you want the complexity for into Geometric_Complexity,
%it outputs a value based on the kolmogorov complexity function
function y = Geometric_Complexity(IslandOutline) 

y=kolmogorov(IslandOutline);
end





function complexity=kolmogorov(s)
n=length(s);
c=1;
l=1;
i=0;
k=1;
k_max=1;
stop=0;
while stop==0
	if s(i+k)~=s(l+k)
        if k>k_max
            k_max=k;
        end
        i=i+1;
        
        if i==l
            c=c+1;
            l=l+k_max;
            if l+1>n
                stop=1;
            else
                i=0;
                k=1;
                k_max=1;
            end
        else
            k=1;
        end
	else
        k=k+1;
        if l+k>n
            c=c+1;
            stop=1;
        end
	end
end
b=n/log2(n);
complexity=c/b;
end