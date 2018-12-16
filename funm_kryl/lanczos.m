function [ v1,H,h,breakdown ] = lanczos( A,m,H,s,param )
%LANCZOS extend a given Lanczos decomposition (V_big,H) of dimension s
%        up to dimension m

breakdown = 0;   % no checking
global V_big
H(m+1,m) = 0;

% orthog first mvp against ALL previous

v0 = V_big(:,s);
if isnumeric(A),
    w = A*v0;
else
    w = A(v0);
end;
    
for j = 1:s,
    ip = param.inner_product(w,V_big(:,j));
    H(j,s) = H(j,s) + ip;
    w = w - V_big(:,j)*ip;
end;
    
H(s+1,s) = sqrt(param.inner_product(w,w));
v1 = (1/H(s+1,s))*w;

% from here on: Lanczos with v0 and v1    
for k = s+1:m,
    H(k-1,k) = H(k,k-1);
    V_big(:,k) = v1;
    if isnumeric(A),
        w = A*v1;
    else
        w = A(v1);
    end;
    w = w - H(k-1,k)*v0;
    H(k,k) = param.inner_product(w,v1);
    w = w - H(k,k)*v1;
    H(k+1,k) = sqrt(param.inner_product(w,w));
    v0 = v1;
    v1 = w*(1/H(k+1,k));
end;
    
h = H(m+1,m);
H = H(1:m,1:m);
