function dS_dav=Sav(K,Y,V,lam)

% S=diag(K*V)*Y*V

volt=abs(V);
theta=angle(V);


dS_dav=sparsediag(K*V)*conj(Y)*sparsediag(lam)*conj(sparsediag(1j*exp(1j*theta)))+sparsediag(conj(Y)*conj(diag(1j*V))*lam)*K*sparsediag(exp(1j*theta))...
    +sparsediag(conj(Y*V))*K*sparsediag(lam)*sparsediag(1j*exp(1j*theta))+sparsediag(K*diag(1j*V)*lam)*conj(Y)*conj(diag(exp(1j*theta)));

end