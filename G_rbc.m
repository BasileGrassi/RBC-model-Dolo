function [G, G_s, G_x, G_e] = G_rbc(s,x,e,model);
% 
%
% none=[];
% model=rbc_matlab('model',none,none,none,none,none,none,none); 

%rbc_matlab(flag,s,x,z,e,snext,xnext,p,out);
if nargout == 1
    G = rbc_matlab('g',s,x,[],e,[],[],model.params);
else
    [G, G_s, G_x] = rbc_matlab('g',s,x,[],e,[],[],model.params);
end


end