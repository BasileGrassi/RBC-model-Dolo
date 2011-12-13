function [F,F_s,F_x,F_S,F_X] =F_rbc(s,x,e,snext,xnext,model)

    if nargout == 1
        z=rbc_matlab('h',s,x,[],e,snext,xnext,model.params);
        F=rbc_matlab('f',s,x,z,e,[],[],model.params);
    else
        [z, z_s, z_x, z_S, z_X]=rbc_matlab('h',s,x,[],e,snext,xnext,model.params);
        [f,f_s,f_x,f_z] = rbc_matlab('f',s,x,z,[],[],[],model.params);
        F = f;
        F_s = f_s + arraymult(f_z, z_s);
        F_x = f_x + arraymult(f_z, z_x);
        F_S = arraymult(f_z, z_S);
        F_X = arraymult(f_z, z_X);
    end

end