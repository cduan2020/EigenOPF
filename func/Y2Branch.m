function [branch, Gs, Bs]=Y2Branch(Ybus)

n=size(Ybus,1);
branch=[];

for i=1:n
    for j=i+1:n
        if Ybus(i,j)~=0
            fbus=i;
            tbus=j;
            r=real(1/(-Ybus(i,j)));
            x=imag(1/(-Ybus(i,j)));
            b=0;
            rateA=Inf;
            rateB=Inf;
            rateC=Inf;
            ratio=0;
            angle=0;
            status=1;
            angmin=-360;
            angmax=360;
            newbrch=[fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax];
            branch=[branch ; newbrch];
        end
    end
end

Gs=real(sum(Ybus,2)); % add Gs to bus(:,5)
Bs=imag(sum(Ybus,2));  % add Gs to bus(:,6)

end
        