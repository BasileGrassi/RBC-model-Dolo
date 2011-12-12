function out=G_rbc(s,x,e,model);
% 
%
% none=[];
% model=rbc_matlab('model',none,none,none,none,none,none,none); 

%rbc_matlab(flag,s,x,z,e,snext,xnext,p,out);
out=rbc_matlab('g',s,x,[],e,[],[],model.params);

end