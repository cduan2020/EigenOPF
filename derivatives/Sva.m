function dS_dva=Sva(K,Y,V,lam)

% S=diag(K*V)*Y*V

volt=abs(V);
theta=angle(V);


dS_dva=sparsediag(K*V)*conj(Y)*sparsediag(lam)*conj(sparsediag(1j*exp(1j*theta)))+sparsediag(conj(Y)*conj(sparsediag(exp(1j*theta)))*lam)*K*sparsediag(1j*V)...
    +sparsediag(conj(Y*V))*K*sparsediag(lam)*sparsediag(1j*exp(1j*theta))+sparsediag(K*sparsediag(exp(1j*theta))*lam)*conj(Y)*conj(sparsediag(1j*V));

end