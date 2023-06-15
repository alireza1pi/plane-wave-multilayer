function [Ell_plus, Ell_minus] = forward_calculation(El_plus, El_minus, dl, kzl, kzll, mul, mull)

a_l = exp(-1i*kzl*dl);
b_l = exp(+1i*kzl*dl);

a_ll = exp(-1i*kzll*dl);
b_ll = exp(+1i*kzll*dl);

pl = kzl/mul;
pll = kzll/mull;

AA_plus = [El_plus*a_l + El_minus*b_l, b_ll; pl*a_l*El_plus - pl*b_l*El_minus, -pll*b_ll];
AA_minus = [a_ll, El_plus*a_l + El_minus*b_l; pll*a_ll, pl*a_l*El_plus -  pl*b_l*El_minus];


DD = [a_ll, b_ll; pll*a_ll, -pll*b_ll];

 Ell_plus = det(AA_plus) / det(DD);
 Ell_minus = det(AA_minus) / det(DD);

% Ell_plus = +(El_plus*a_l*(pll+pl)  + El_minus*b_l*(pll - pl)) /  (2*pll*a_ll);

% Ell_minus = -(El_plus*a_l*(pll-pl)  + El_minus*b_l*(pll + pl)) /  (2*pll*b_ll);

AA = [a_ll b_ll; pll*a_ll -pll*b_ll] \ [El_plus*a_l + El_minus*b_l; pl*(El_plus*a_l - El_minus*b_l)];

Ell_plus = AA(1);
Ell_minus = AA(2);


end