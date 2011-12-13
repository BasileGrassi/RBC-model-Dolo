function out=F_rbc(s,x,e,snext,xnext,model)




z=rbc_matlab('h',s,x,[],e,snext,xnext,model.params);

%f 'rond' h
out=rbc_matlab('f',s,x,z,e,[],[],model.params);


end