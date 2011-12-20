function [x] = newton_solve(fun, x0)

    N = size(x0,1);
    n_x = size(x0,2);
        
    tol = 1e-8;
    
    err = 1;
    maxit = 1;
    it = 0;
    
    while err > tol && it < maxit
        
        [res, jac] = fun(x0);
        
        jac = permute(jac,[2,3,1]);
        minv = zeros(size(jac));
        for i = 1:N
            mat = jac(:,:,i);
            minv(:,:,i) = inv( mat );
        end
        minv = permute(minv,[3,1,2]);
        
        x = x0 - arraymult(  minv, res);
        
        err = abs(max(max(res)));
        it = it + 1;
        x0 = x;
    end
        
end