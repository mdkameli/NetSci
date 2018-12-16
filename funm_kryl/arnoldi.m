function [ w,H,h,breakdown ] = arnoldi( A,m,H,s,param )
%ARNOLDI extend a given Arnoldi decomposition (V_big,H) of dimension s
%        up to dimension m



global V_big
H(m+1,m) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;

for k = s:m,
    
    w = V_big(:,k);
    if isnumeric(A),
        w = A*w;
    else
        w = A(w);
    end;
    
    sj = max([1,k-trunc+1]):k;
    if(k==s), sj = 1; end
    for r = 0:reo,
        for j = sj:k,
            ip = param.inner_product(w,V_big(:,j));
            H(j,k) = H(j,k) + ip;
            w = w - V_big(:,j)*ip;
        end;
    end;
    
    H(k+1,k) = sqrt(param.inner_product(w,w));
    
    if H(k+1,k) < k*eps,
        breakdown = 1;
        break;
    end
    
    w = (1/H(k+1,k))*w;
    if k < m
        V_big(:,k+1) = w;
    end;
end;

h = H(m+1,m);
H = H(1:m,1:m);
