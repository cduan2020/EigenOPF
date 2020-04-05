function X= sparsediag(vec)
n=length(vec);
X = spdiags(vec,0,n,n);
end