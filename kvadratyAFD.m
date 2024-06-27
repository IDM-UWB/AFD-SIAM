function Vnew=kvadratyAFD(V)
  %
% proklada vlacek kvadratickymi funkcemi, u delsich dimenzi dvemi, jinak
% jednou
%
  N=length(V);
  Vnew=V;
  for n=1:N
    [R1,N1,R2]=size(V{n});
    if N1<10   %%%       kratke dimenze prokladam jednou parabolou
               %%%       asi to jsou ty kovariancni matice
      M1=[ones(N1,1),(0:(N1-1))',((0:(N1-1)).^2)'];
      M=inv(M1'*M1);
      for r1=1:R1
        for r2=1:R2
          Vnew{n}(r1,:,r2)=proloz(V{n}(r1,:,r2),M);
        end
      end
    else     %%%   dlouhe dimenze prokladam dvemi parabolami
             %%%   to je pro ty stredni hodnoty a pravdepodobnosti.
      N2=N1; N1=floor((N2+1)/2);
      M1=[ones(N1,1),(0:(N1-1))',((0:(N1-1)).^2)'];
      M=inv(M1'*M1);
      for r1=1:R1
        for r2=1:R2
          Vnew{n}(r1,1:N1,r2)=proloz(V{n}(r1,1:N1,r2),M);
          Vnew{n}(r1,N2-N1+1:N2,r2)=proloz(V{n}(r1,N2-N1+1:N2,r2),M);
        end
      end
    end
  end
end
%
function xr=proloz(x,M)
  %
  N1=length(x);
  x0=sum(x);
  x1=sum(x(:)'.*(0:N1-1));
  x2=sum(x(:)'.*(0:N1-1).^2);
  ab=M*[x0; x1; x2];
  xr=ab(3)*(0:N1-1).^2+ab(2)*(0:N1-1)+ab(1);
end
