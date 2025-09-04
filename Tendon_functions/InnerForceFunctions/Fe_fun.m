function Fe = Fe_fun(Fe0,uk,qk,qk0,phik,Phik,Fibers,Dvec,ElemDofs,PosDofs,Gint,Nint,detF0) 
L = Dvec(end);   % element length
Fe = Fe0;
for ii=1:Nint    % integration all over the element's volume
    Fe = Fe + FedV(Fe0,uk,qk,qk0,phik,Phik,Fibers,Dvec,ElemDofs,PosDofs,Gint(ii,1),Gint(ii,2),Gint(ii,3))*Gint(ii,4);
end
Fe=Fe*1/2*L*detF0;
