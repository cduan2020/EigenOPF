function dS_dvv=Svv(K,Y,V,lam)

% S=diag(K*V)*Y*V

volt=abs(V);
theta=angle(V);


dS_dvv=sparsediag(conj(Y)*conj(sparsediag(exp(1j*theta)))*lam)*K*sparsediag(exp(1j*theta))...
    +sparsediag(K*sparsediag(exp(1j*theta))*lam)*conj(Y)*conj(sparsediag(exp(1j*theta)));

end