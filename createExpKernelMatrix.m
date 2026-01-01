function K = createExpKernelMatrix(X, c)
    % createExpKernelMatrix computes the RBF kernel matrix for data X
    % X: N x D matrix, where N is the number of samples and D is the dimension
    % c: kernel parameter

    % Compute pairwise Euclidean distances
    % pdist2(X, X) calculates the distance between every pair of rows in X.
    r = pdist2(X, X); 
    
    % Compute the kernel matrix using the exponential formula
    K = exp(-r / c); 
end