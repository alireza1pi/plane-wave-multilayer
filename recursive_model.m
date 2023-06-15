function El = recursive_model(Ell, dl, kzl, kzll, mul, mull)

a_l = exp(-1i*kzl*dl);
b_l = exp(+1i*kzl*dl);

a_ll = exp(-1i*kzll*dl);
b_ll = exp(+1i*kzll*dl);

pl = kzl/mul;
pll = kzll/mull;

E0 = [a_ll + Ell*b_ll; pll*(a_ll - Ell*b_ll)];

AA_left = [a_l; pl*a_l];
AA_right = [b_l; -pl*b_l];

El_minus = det([AA_left, E0]);
El_plus =det( [E0, AA_right]);

 El = El_minus / El_plus;

 El =     1/ ((b_l/ a_l)     *    (b_ll*Ell*(pl-pll) + a_ll*(pl+pll))      /     (b_ll*Ell*(pl+pll) + a_ll*(pl-pll)))  ;
 
 AA = [a_l, b_l; pl*a_l, -pl*b_l] \ [a_ll+b_ll*Ell; pll*a_ll-pll*b_ll*Ell];
 
 El = AA(2)/AA(1);
end