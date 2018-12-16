function [ ell,U,T,D ] = thick_sym( H,U,T,D,cycle,h_f )
%THICK_SYM returns reordered Schur form of Hermitian H
%   [U,T] = schur(H) and D = ordeig(H)
%   return reordered Schur form and eigenvalues such that
%   the ell "wanted" ones occur first
%
%   If the restart-length is large (>100), ordschur may have
%   stability problems.
%   If A is Hermitian, the Schur form is diagonal and
%   reordering can be done stably by permutation of columns.

ell = 10;

[ignore,jj] = sort(real(D),1,'ascend');

T = diag(T);
T = diag(T(jj));
U = U(:,jj);
D = D(jj);

