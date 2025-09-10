function [A, b] = build_sparse_A(Z0, iU)
    [n, p] = size(Z0);
    q = size(iU,2) - 1;
    dim = (p-1)*(q+1);
    A = sparse(n*p, p*dim);
    b = [];
    interM = [];
    for j=1:(q+1)
        interM=[interM Z0.*iU(:,j)];
    end
    for j=1:p
        Xj = interM;
        Xj(:,((0:q)*p+j)) = [];
        A(((j-1)*n+1):(j*n),((j-1)*dim+1):(j*dim)) = Xj;
        y1 = Z0(:,j);
        b = [b; y1 - mean(y1)];
    end
end