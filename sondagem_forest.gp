P = polcyclo(16);
z = Mod(x, P);
r2 = z^2 - z^6;
I4 = z^4;
cj(a) = Mod(subst(lift(a), x, x^15), P);
dag(M) = matrix(2,2,i,j, cj(M[j,i]));
H = [1,1;1,-1]/r2;
S = [1,0;0,I4];
Sd = dag(S);
Uz(a) = [1,0;0,z^a];
Ux(a) = H*Uz(a)*H;
Uy(a) = S*Ux(a)*Sd;
Urot(p,a) = if(p==1, Ux(a), if(p==2, Uy(a), Uz(a)));
sig = [[0,1;1,0], [0,-I4;I4,0], [1,0;0,-1]];
bloch(U) = matrix(3,3,i,j, trace(sig[i]*U*sig[j]*dag(U))/2);
dexp(e) = if(e==0, -10^6, -valuation(norm(e),2)/2);
Nmax(B) = vecmax(apply(dexp, concat(Vec(B))));
print("sanidade: r2^2 = ", r2^2, " ; H*H = ", H*H, " ; N(Uz(2)) = ", Nmax(bloch(Uz(2))), " ; N(Uz(1)) = ", Nmax(bloch(Uz(1))), " ; N(Ux(3)) = ", Nmax(bloch(Ux(3))));
palavra(mmax) = { my(m = 1+random(mmax), U = matid(2), p = 0, q, a, nT = 0, nR = 0); for(i=1, m, q = 1+random(3); while(q==p, q = 1+random(3)); p = q; a = 1+random(3); if(a==2, nT++, nR++); U = U*Urot(p,a)); [Nmax(bloch(U)), nT, nR, m]; }
setrand(20261114);
{ bad = 0; for(t=1, 3000, w = palavra(12); if(w[1] != 2*w[2]+3*w[3], bad++; if(bad<=5, print("DIVERGE ", w)))); }
print("palavras = 3000 | divergencias de N = 2 nT + 3 nR: ", bad);
setrand(1);
{ bad2 = 0; for(t=1, 300, w = palavra(8); if(w[1] != w[2]+3*w[3], bad2++)); }
print("controle (formula errada nT + 3 nR): falhas = ", bad2, "/300");
quit
