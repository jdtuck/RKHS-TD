function K = createRBfKernelMatrix(X, gamma)
    % createRBFKernelMatrix computes the RBF kernel matrix for data X
    % X: N x D matrix, where N is the number of samples and D is the dimension
    % gamma: RBF kernel parameter (sigma relates to gamma, gamma = 1/(2*sigma^2))

    % Calculate the squared Euclidean distances between all pairs of points
    % This is a highly optimized vectorized method:
    % (a-b)^2 = a^2 + b^2 - 2ab
    % sum(X.^2, 2) calculates the squared L2 norm of each row (a^2 part)
    nsq = sum(X.^2, 2);
    
    % K now holds the squared distances matrix
    % bsxfun handles broadcasting for the matrix operations
    K = bsxfun(@minus, nsq, 2 * (X * X.')); % Computes a^2 - 2ab
    K = bsxfun(@plus, nsq.', K);           % Adds b^2: a^2 - 2ab + b^2 = (a-b)^2
    
    % Apply the exponential function to get the RBF kernel matrix
    K = exp(-gamma * K);
end