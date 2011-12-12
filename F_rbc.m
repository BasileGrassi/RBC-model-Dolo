function out=F_rbc(s,x,e,snext,xnext,model);

%none=[];
%model=rbc_matlab('model',none,none,none,none,none,none,none); 

%rbc_matlab(flag,s,x,z,e,snext,xnext,p,out);
h=rbc_matlab('h',s,x,[],e,snext,xnext,model.params);

%f 'rond' h
f=rbc_matlab('f',s,x,h,e,[],[],model.params);

out=f;
end