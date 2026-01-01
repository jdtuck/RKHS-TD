addpath("~/Documents/Research/tensor_toolbox/")
%% Generate Data
f = @(x) normpdf(linspace(0,1,99), sin(2*pi*x(1)^2)/4 - sqrt(x(1)*x(1))/10 + .5, 0.05) * x(2);

n = 100;
nt = 99;
p = 3;
x_train = rand(p, n, p);
e = randn(1, n*99);
A = zeros(p,n,nt);
for i = 1:p
    for j = 1:n
        A(i,j,:) = f(x_train(i,j,:));
    end
end

tt = linspace(0, 1, nt);

X = tensor(A);