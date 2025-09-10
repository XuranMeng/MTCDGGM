function [y, beta] = mulGMMReg(Z0, iU,pools)
    [n, p] = size(Z0);
    q = size(iU,2) - 1;
    dim = (p-1)*(q+1);

    %% cross-validation
    lam1_max = NaN; 
    alpha = 0.75; 
    nl = 20; % lambda grid
    cut = repelem(1:5, n/5);
    Cerror = zeros(5, nl);


    lambda1 = NaN(1, nl);
    lambda2 = NaN(1, nl);

   
    parpool('local', pools);
    eigsopt.issym = 1;
    eigsopt.tol = 1e-3;
    for t = 1:5
        disp(['CV fold: ', num2str(t)]);


        idx_train = (cut ~= t);
        Z0_train = Z0(idx_train, :);
        iU_train = iU(idx_train, :);
        [Atrain, btrain] = build_sparse_A(Z0_train, iU_train);


        idx_test = (cut == t);
        Z0_test = Z0(idx_test, :);
        iU_test = iU(idx_test, :);
        [Atest, btest] = build_sparse_A(Z0_test, iU_test);
        
        disp(['Still Ok now']);
 
          
        lam1_max = norm(Atrain' * btrain, Inf); 
        lambda1 = lam1_max * exp(linspace(log(0.8), log(0.3), nl));
        lambda2 = lambda1 * (1 - alpha) / alpha;
     
        %% operator
        Amap = @(x) mexMatvec(Atrain, x, 0);
        ATmap = @(y) mexMatvec(Atrain, y, 1);
        AATmap = @(x) Amap(ATmap(x));

        
        disp(['Still Ok now']);
        Lip = normest(Atrain,0.1)^2;
        disp(['Still Ok now, have Lip']);
        opts.stoptol = 1e-4;
        opts.printyes = 0;
        opts.Lip = Lip;
        Ainput.A = Atrain;
        Ainput.Amap = Amap;
        Ainput.ATmap = ATmap;

        %% group structure
        int = repmat(ceil((1:dim)/(p-1)), 1, p);
        G = [];
        for i = 1:(q+1)
            G = [G; find(int == i)'];
        end
        ind = zeros(3, q+1);
        grpsize = round(size(Atrain,2)/(q+1)); 
        for i = 1:(q+1)
            if i == 1
                ind(1,1) = 1; ind(2,1) = grpsize; ind(3,1) = 0;
            else
                ind(1,i) = ind(2,i-1) + 1;
                ind(2,i) = i*grpsize;
                ind(3,i) = 1;
            end
        end

        %% solver
        disp(['Starting par']);
        parfor i = 1:nl
            disp([' Lambda grid: ', num2str(i)]);
            c = [lambda1(i); lambda2(i)];
            [obj, yval, z, x, info, runhist] = SGLasso_SSNAL(Ainput, btrain, size(Atrain,2), c, G, ind, opts);
            res = Atest * x - btest;
            Cerror(t,i) = Cerror(t,i) + sum(res.^2);
        end
    end

    delete(gcp('nocreate'));
    
    meanC = mean(Cerror, 1);
    [~, minid] = min(meanC);

    %% 
    [A, b] = build_sparse_A(Z0, iU);

    lam1_max = norm(A' * b, Inf); 
    lambda1 = lam1_max * exp(linspace(log(0.8), log(0.3), nl));
    lambda2 = lambda1 * (1 - alpha) / alpha;
    c = [lambda1(minid); lambda2(minid)];

    Amap = @(x) mexMatvec(A, x, 0);
    ATmap = @(y) mexMatvec(A, y, 1);
    AATmap = @(x) Amap(ATmap(x));

    Lip = normest(A,0.1)^2;

    opts.Lip = Lip;
    opts.stoptol = 1e-4;
    Ainput.A = A;
    Ainput.Amap = Amap;
    Ainput.ATmap = ATmap;

    int = repmat(ceil((1:dim)/(p-1)), 1, p);
    G = [];
    for i = 1:(q+1)
        G = [G; find(int == i)'];
    end
    ind = zeros(3, q+1);
    grpsize = round(size(A,2)/(q+1)); 
    for i = 1:(q+1)
        if i == 1
            ind(1,1) = 1; ind(2,1) = grpsize; ind(3,1) = 0;
        else
            ind(1,i) = ind(2,i-1) + 1;
            ind(2,i) = i*grpsize;
            ind(3,i) = 1;
        end
    end

    [obj, y, z, x, info, runhist] = SGLasso_SSNAL(Ainput, b, size(A,2), c, G, ind, opts);
    beta = x;
end