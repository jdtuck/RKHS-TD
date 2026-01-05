function f = gen_fourier(x, r, theta)
v = zeros(length(x), 10);
for j = 1:length(x)
    for i = 1:10
        if (mod(i,2) == 0)
            v(j,i) = theta(r,i) * sqrt(2) * cos((2-1)*pi*x(j));
        else
            v(j,i) = theta(r,i) * sqrt(2) * sin((2-1)*pi*x(j));
        end
    end
end

f = sum(v,2);