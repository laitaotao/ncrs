function v = guidedrandsample2(z,k,w)
v = zeros(1,k);
i = 0;
c = 0;
while(i<k)
    i = i+1;
    x = rand();
    a = find(z >x,1);
    if(find(v == a))
        i = i-1;
    else
        v(i)=a;
    end
    c = c + 1;
    if(c>= 5*k)
        v = guidedrandsample3(w,k);%,true,w,'false');
        break;
    end
end
end

%  z = cumsum(w(:)'/sum(w));


function v = guidedrandsample3(w,k)
z = cumsum(w+eps)/sum(w+eps);
v = zeros(1,k);
i = 0;

while(i<k)
    i = i+1;
    x = rand();
    a = find(z >x,1);
    if(find(v == a))
        i = i-1;
    else
        v(i)=a;        
        w = w+eps;
        w(v(1:i)) = 0;
        if(i<k)
            z = cumsum(w)/sum(w);
        end
    end
end
end