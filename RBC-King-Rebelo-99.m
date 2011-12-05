function [out1,out2,out3,out4,out5] = mf_RBC-King-Rebelo-99(flag,s,x,z,e,snext,xnext,p,out);

output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'Jz',0,'hmult',0);

if nargin == 9
  output                     = catstruct(output,out);
  voidcell                   = cell(1,5);
  [out1,out2,out3,out4,out5] = voidcell{:};
else
  if nargout >= 2, output.Js = 1; end
  if nargout >= 3, output.Jx = 1; end
  if nargout >= 4
    if strcmpi(flag, 'f')
      output.Jz = 1;
    else
      output.Jsn = 1;
    end
  end
  if nargout >= 5, output.Jxn = 1; end
end


switch flag

  case 'b';
    n = size(s,1);
    out1 = zeros(n,2);
    out1(:,1) = 0;
    out1(:,2) = 0;
    out2 = zeros(n,2);
    out2(:,1) = inf;
    out2(:,2) = inf;

  case 'f';
    n = size(s,1);

    % f
    if output.F
      out1 = zeros(n,2);
      out1(:,1) = z(:,1) - 1;
      out1(:,2) = -p(4).*x(:,1).^p(2).*x(:,2).^p(3) + s(:,1).*(s(:,2)./x(:,2)).^p(6).*(-p(6) + 1);
    end

    % df/ds
    if output.Js
      out2 = zeros(n,2,2);
      out2(:,1,1) = 0; % d eq_1 w.r.t. z
      out2(:,1,2) = 0; % d eq_1 w.r.t. k
      out2(:,2,1) = (s(:,2)./x(:,2)).^p(6).*(-p(6) + 1); % d eq_2 w.r.t. z
      out2(:,2,2) = p(6).*s(:,1).*(s(:,2)./x(:,2)).^p(6).*(-p(6) + 1)./s(:,2); % d eq_2 w.r.t. k
    end

    % df/dx
    if output.Jx
      out3 = zeros(n,2,2);
      out3(:,1,1) = 0; % d eq_1 w.r.t. c
      out3(:,1,2) = 0; % d eq_1 w.r.t. n
      out3(:,2,1) = -p(4).*p(2).*x(:,1).^p(2).*x(:,2).^p(3)./x(:,1); % d eq_2 w.r.t. c
      out3(:,2,2) = -p(6).*s(:,1).*(s(:,2)./x(:,2)).^p(6).*(-p(6) + 1)./x(:,2) - p(4).*p(3).*x(:,1).^p(2).*x(:,2).^p(3)./x(:,2); % d eq_2 w.r.t. n
    end

    % df/dz
    if output.Jz
      out4 = zeros(n,2,1);
      out4(:,1,1) = 1; % d eq_1 w.r.t. Eh
      out4(:,2,1) = 0; % d eq_2 w.r.t. Eh
    end
        
  case 'g';
    n = size(s,1);

    % g
    if output.F
      out1 = zeros(n,2);
      out1(:,1) = p(7).*s(:,1) + p(8).*(-p(7) + 1) + e(:,1);
      out1(:,2) = -x(:,1) + s(:,2).*(p(6).*snext(:,1).*(xnext(:,2)./snext(:,2)).^(-p(6) + 1) - p(5) + 1);      
    end

    if output.Js
      out2 = zeros(n,2,2);
      out2(:,1,1) = p(7); % d eq_1 w.r.t. z
      out2(:,1,2) = 0; % d eq_1 w.r.t. k
      out2(:,2,1) = 0; % d eq_2 w.r.t. z
      out2(:,2,2) = p(6).*snext(:,1).*(xnext(:,2)./snext(:,2)).^(-p(6) + 1) - p(5) + 1; % d eq_2 w.r.t. k
    end

    if output.Jx
      out3 = zeros(n,2,2);
      out3(:,1,1) = 0; % d eq_1 w.r.t. c
      out3(:,1,2) = 0; % d eq_1 w.r.t. n
      out3(:,2,1) = -1; % d eq_2 w.r.t. c
      out3(:,2,2) = 0; % d eq_2 w.r.t. n
    end
        
  case 'h';
    n = size(snext,1);

    %h
    if output.F
      out1 = zeros(n,1);
      out1(:,1) = p(1).*(x(:,1)./xnext(:,1)).^p(2).*(p(6).*snext(:,1).*(xnext(:,2)./snext(:,2)).^(-p(6) + 1) - p(5) + 1);
    end

    if output.Js
      out2 = zeros(n,1,2);
      out2(:,1,1) = 0; % d eq_1 w.r.t. z
      out2(:,1,2) = 0; % d eq_1 w.r.t. k
    end

    if output.Jx
      out3 = zeros(n,1,2);
      out3(:,1,1) = p(1).*p(2).*(x(:,1)./xnext(:,1)).^p(2).*(p(6).*snext(:,1).*(xnext(:,2)./snext(:,2)).^(-p(6) + 1) - p(5) + 1)./x(:,1); % d eq_1 w.r.t. c
      out3(:,1,2) = 0; % d eq_1 w.r.t. n
    end

    if output.Jsn
      out4 = zeros(n,1,2);
      out4(:,1,1) = p(6).*p(1).*(x(:,1)./xnext(:,1)).^p(2).*(xnext(:,2)./snext(:,2)).^(-p(6) + 1); % d eq_1 w.r.t. z(1)
      out4(:,1,2) = -p(6).*p(1).*snext(:,1).*(x(:,1)./xnext(:,1)).^p(2).*(xnext(:,2)./snext(:,2)).^(-p(6) + 1).*(-p(6) + 1)./snext(:,2); % d eq_1 w.r.t. k(1)
    end

    if output.Jxn
      out5 = zeros(n,1,2);
      out5(:,1,1) = -p(1).*p(2).*(x(:,1)./xnext(:,1)).^p(2).*(p(6).*snext(:,1).*(xnext(:,2)./snext(:,2)).^(-p(6) + 1) - p(5) + 1)./xnext(:,1); % d eq_1 w.r.t. c(1)
      out5(:,1,2) = p(6).*p(1).*snext(:,1).*(x(:,1)./xnext(:,1)).^p(2).*(xnext(:,2)./snext(:,2)).^(-p(6) + 1).*(-p(6) + 1)./xnext(:,2); % d eq_1 w.r.t. n(1)
    end
        
  case 'e';
    warning('Euler equation errors not implemented in Dolo');

  case 'params';
    out1 = [0.99,1,1,4,0.025,0.33,0.95,1];

  case 'solution';


end
