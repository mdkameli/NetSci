function [ out ] = eval_parfrac( r,A,b )
%EVAL_PARFRAC Evaluate partial fraction returned by parfrac
%  EVAL_PARFRAC(r,A)    returns r(A)
%  EVAL_PARFRAC(r,A,b)  returns r(A)*b

if nargin == 2   % matrix case
    if issparse(A),
        I = speye(size(A));
    else
        I = eye(size(A));
    end;
    out = r.absterm*I;
    for j = 1:length(r.single_poles)
        out = out + r.single_coeff(j)*(r.single_poles(j)*I - A)^-1;
    end;
    if isreal(A),
        for j = 1:length(r.conj_poles)
            out = out + 2*real(r.conj_coeff(j)*(r.conj_poles(j)*I - A)^-1);
        end;
    else
        for j = 1:length(r.conj_poles)
            out = out + r.conj_coeff(j)*(r.conj_poles(j)*I - A)^-1;
            out = out + conj(r.conj_coeff(j))*(conj(r.conj_poles(j))*I - A)^-1;
        end;
    end;
end;

if nargin == 3   % matrix-vector case
    if issparse(A),
        I = speye(size(A));
    else
        I = eye(size(A));
    end;
    out = r.absterm*b;
    for j = 1:length(r.single_poles)
        out = out + r.single_coeff(j)*((r.single_poles(j)*I - A)\b);
    end;
    if isreal(A),
        for j = 1:length(r.conj_poles)
            out = out + 2*real(r.conj_coeff(j)*((r.conj_poles(j)*I - A)\b));
        end;
    else
        for j = 1:length(r.conj_poles)
            out = out + r.conj_coeff(j)*((r.conj_poles(j)*I - A)\b);
            out = out + conj(r.conj_coeff(j))*((conj(r.conj_poles(j))*I - A)\b);
        end;
    end;
end;
    
    