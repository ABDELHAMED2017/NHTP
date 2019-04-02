function [out1,out2] = compressed_sensing(x,flag,T1,T2,data)


    Tx = find(x); 
    Ax = data.A(:,Tx)*x(Tx);
 
    switch flag        
    case 'ObjGrad'
        Axb  = Ax-data.b;
        out1 = sum(Axb.*Axb)/2;                %objective function 
        if  nargout>1 
        out2 = data.At*Axb;                    %gradien  
        end
    case 'Hess'
        out1 = data.At(T1,:)*data.A(:,T1);     %submatrix  containing T1 rows and T1 columns of Hessian
        if nargout>1 
        out2 = data.At(T1,:)*data.A(:,T2);     %submatrix  containing T1 rows and T2 columns of Hessian
        end
  
    end

end


