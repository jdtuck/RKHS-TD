addpath("../tensor_toolbox/")
%% Generate Data

nt = 100;
tt = linspace(0, 1, nt);
p = 51;
n = 60;
R = 5;

a = zeros(n,R);
b = zeros(p,R);
for r = 1:R
    a(:,r) = rand(1,n);
    b(:,r) = rand(1,p);
end

theta = zeros(R, 10);
for i = 1:10
    theta(:,i) = -1/i + (1/i-(-1/i)).*rand(R,1);
end

A = zeros(n,p,nt);
for i = 1:n
    data_tmp = zeros(p, nt);
    for j = 1:p
        lambda = zeros(nt,1);
        for r = 1:R
            lambda = lambda + 10*sqrt(r) * a(i,r) * b(j,r) * gen_fourier(tt, r, theta);
        end
        A(i,j,:) = lambda;
    end
end

X = tensor(A);

[P, fititer] = cp_hifi(X, tt, 2);

plot(P.u{3})
