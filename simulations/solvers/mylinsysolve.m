    function q = mylinsysolve(L,r) 
       if strcmp(L.matfct_options,'chol') %Cholesky 分解
          q(L.perm,1) = mextriang(L.R, mextriang(L.R,r(L.perm),2) ,1);
       elseif strcmp(L.matfct_options,'spcholmatlab') %稀疏Cholesky 分解
          q(L.perm,1) = mexbwsolve(L.Rt,mexfwsolve(L.R,r(L.perm,1)));
       end
    end

    