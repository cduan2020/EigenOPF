function maximumValues=maxk(Array,n)


Arraycopy = Array;
for j = 1:n
   [a, Index(j)] = max(Arraycopy);
   Arraycopy(Index(j)) = -inf;
end
maximumValues = Array(Index);

end