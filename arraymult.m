function c = arraymult(a, b, me, you, him, our)

    if ndims(b) == 3
        n = size(a,1);
        p = size(a,2);
        q = size(a,3);
        r = size(b,3);
        c = zeros( n, p, r );
        for i = 1:p
            for j = 1:r
                for k = 1:q
                    c(:,i,j) = c(:,i,j) + a(:,i,k) .* b(:,k,j);
                end
            end
        end
    end
    
    if ndims(b) == 2
        n = size(a,1);
        p = size(a,2);
        q = size(a,3);
        %r = size(b,2);
        c = zeros( n, p );
        for i = 1:p
            for j = 1:q
                c(:, i) = c(:, i) + a(:, i, j) .* b(:,j);
            end
        end
    end
end
