% a = sym('a','real');
% b = sym('b','real');
% c = sym('c','real');
% lambda = sym('lambda','real');
% syms z
% A=[1-a*z,-c*z;0,1-b*z]
% H=inv(A)
% S=H*H'
% fY_X=log(det(S(1,1))/det(H(1,1)*H(1,1)'))
% fYX_lambda=subs(fY_X, z, cos(lambda)-1i*sin(lambda))
% fYX_lambda=simplify(fYX_lambda,'Steps', 100)
%% frequency
a = sym('a','real');
b = sym('b','real');
c = sym('c','real');
lambda = sym('lambda','real');
syms z
d = sym('d','real');
A=[1-a*z,-c*z;0,1-b*z-d*z^2]
H=inv(A)
S=H*H'
fY_X=log(det(S(1,1))/det(H(1,1)*H(1,1)'))
fYX_lambda=subs(fY_X, z, cos(lambda)-1i*sin(lambda))
fYX_lambda=simplify(fYX_lambda,'Steps', 100)
fYX_lambda=simplify(fYX_lambda,'Steps', 1000)