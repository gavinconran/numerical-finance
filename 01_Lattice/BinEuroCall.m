% BinEuroCall.m
## Not working

function price = BinEuroCall(S0, K, r, sigma, T, n)

delta_T = T / n;
u = exp(sigma*sqrt(delta_T));
d = 1/u;
p = (exp(r*delta_T) - d) / (u - d);

lattice = zeros(n+1, n+1);
printf("n = %d", n);

for i=0:n
    if i > 1
        delta_T_i = T / i;
        u_i = exp(sigma*sqrt(delta_T_i));
        #d_i = 1 / u_i;
        lattice(i+1,n+1) = max(0, u*d*S0 -K);
    else
        lattice(i+1,n+1) = K;
    endif;    

% Go backwards
for k=n-1:-1:0
    for j=0:k   
        lattice(j+1,k+1) = p*lattice(j+2,k+2) + (1-p)*lattice(i+1,k+2);
    end;
end;    

price = exp(-r*T)*lattice(1,1);

