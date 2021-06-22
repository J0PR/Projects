function Z= rightDivision(X,Y)
sx = size(X);
dx = ndims(X);
Xt = reshape(permute(X, [1 3:dx 2]), [prod(sx)/sx(2) sx(2)]);
Z = Xt/Y;
Z = permute(reshape(Z, sx([1 3:dx 2])), [1 dx 2:dx-1]);