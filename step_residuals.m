function [ R, R_x ] = step_residuals( s, x, e, w, parms, f, g, model, coeff, cdef )
%RESIDUALS Summary of this function goes here
%   Detailed explanation goes here

    
%    [coeff,B]=funfitxy(cdef, s, x);
    
    %xnext=funeval(coeff, cdef, snext);
    
    nobs = size(s,1);
    n_s = size(s,2);
    n_x = size(x,2);
    
    n_e = size(e,1);
    
        
    ind   = (1:nobs);
    ind   = ind(ones(1,n_e),:);
    ss    = s(ind,:);   %%%% each observation in ss and xx
    xx    = x(ind,:);   %%%%  is repeated K times
    ee    = e(repmat(1:n_e,1,nobs),:);
    
    
    
    
    [SS, SS_ss, SS_xx] = g(ss, xx, ee, model);
    
    derivs = [zeros(1,n_x); eye(n_x)];
    
    M = funeval(coeff, cdef, SS, derivs);
    
    XX = M(:,:,1);
    XX_SS = M(:,:,2:end);
    
    
    [FF, FF_ss, FF_xx, FF_SS, FF_XX] = f(ss, xx, ee, SS, XX, model);
    
    RR = FF;
    RR_xx = FF_xx + arraymult( FF_SS, SS_xx) + arraymult(FF_XX, arraymult( XX_SS, SS_xx ));
    
    tmp = reshape( ss, n_e, nobs * n_s);
    
    
    R = reshape( w' * reshape( RR, n_e, nobs * n_x ),  nobs, n_x);
    R_x = reshape( w' * reshape( RR_xx, n_e, nobs * n_x * n_x ),  nobs, n_x, n_x);
    
    %for i = 1:n_e
    %   R = R 
    %end
    


end

