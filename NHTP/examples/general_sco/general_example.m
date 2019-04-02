function [out1,out2] = general_example(z,flag,T1,T2,data)
    switch flag
    case 'ObjGrad'             
        out1 = data(z,'obj');  % objective function   
        if nargout>1
        out2 = data(z,'grad');  % gradient
        end       
    case 'Hess'
        H    = data(z,'hess'); 
        out1 = H(T1,T1);       % submatrix containing T1 rows and T1 columns of Hessian
        if nargout>1 
        out2 = H(T1,T2);       % submatrix containing T1 rows and T2 columns of Hessian
        end
    end
     
end



