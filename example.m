addpath("~/Documents/Research/tensor_toolbox/")
%% Generate Data
f = @(x) x(1)*cos(2*pi*linspace(0,1,99))+x(2)*sin(2*pi*linspace(0,1,99))+x(3)*cos(2*pi*2*linspace(0,1,99))+x(4)*sin(2*pi*2*linspace(0,1,99));

n = 100;
nt = 99;
p = 4;
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