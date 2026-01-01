% inputs, X is a tensor
R = 3;  % rank
gamma = 50;
lam = .001;

fitchangetol = 1e-4;
maxiters = 50;
dimorder = 1:N;
printitn = true;

N = ndims(X);
normX = norm(X);

% create Kernel Matrix
K = createRBfKernelMatrix(X(:,:,3),gamma);

Uinit = cell(N,1);
for n = 1:N
    Uinit{n} = rand(size(X,n),R);
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;

% Store the last MTTKRP result to accelerate fitness computation.
U_mttkrp = zeros(size(X, dimorder(end)), R);

% All ones vector for inner product computation.
e = ones(1,size(X, dimorder(end)));

% Create temp storage for computing norm(P) efficiently.
tmpvecs = zeros(R^2,N+1);

UtU = zeros(R,R,N);
for n = 1:N
    if ~isempty(U{n})
        UtU(:,:,n) = U{n}'*U{n};
    end
end

for iter = 1:maxiters

    fitold = fit;

    for n = 1:3  % currently assume first two are discrete, third is functional
        if n < 3
            % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
            Unew = mttkrp(X,U,n);

            % Compute the matrix of coefficients for linear system
            Y = prod(UtU(:,:,[1:n-1 n+1:N]),3);
            Unew = Unew / Y;
            if issparse(Unew)
                Unew = full(Unew);   % for the case R=1
            end

            U{n} = Unew;
            UtU(:,:,n) = U{n}'*U{n};
        else
            % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
            Unew = mttkrp(X,U,n);
            y = vec(Unew);

            % Compute the matrix of coefficients for linear system
            Y = prod(UtU(:,:,[1:n-1 n+1:N]),3);

            M = kron(Y,K) + lam * eye(R*size(X(:,:,3),1));

            w = M / y;

            Unew = reshape(w,size(X(:,:,3),1), R);

            U_mttkrp = Unew;

            U{n} = Unew;
        end

    end

    if iter == 1
        lambda = sqrt(sum(Unew.^2,1))'; %2-norm
    else
        lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
    end

    P = ktensor(lambda,U);

    % This is equivalent to innerprod(X,P).
    iprod = e * (P.U{dimorder(end)} .* U_mttkrp) * lambda;

    % This is equivalent to norm(P)^2
    tmpvecs(:,1:N) = reshape(UtU,[],N);
    tmpvecs(:,N+1) = reshape(lambda*lambda',[],1);
    normPsqr = sum( prod(tmpvecs,2) ) ;

    if normX == 0
        fit = normPsqr - 2 * iprod;
    else
        normresidual = sqrt( normX^2 + normPsqr - 2 * iprod );
        fit = 1 - (normresidual / normX); %fraction explained by model
    end
    fitchange = abs(fitold - fit);

    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        flag = 0;
    else
        flag = 1;
    end

    if (mod(iter,printitn)==0) || ((printitn>0) && (flag==0))
        fprintf(' Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange);
    end

    % Check for convergence
    if (flag == 0)
        break;
    end

end
