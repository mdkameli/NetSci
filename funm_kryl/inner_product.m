function [ w ] = inner_product( a,b )
%INNER_PRODUCT used by arnoldi or lanczos
%   the Krylov basis is orthonormalized with this

w = b'*a;