using MTH229

cm = 1e-2;
m = 1;

ρ = 1

f(x) = x^2 + z^2
@syms x
@syms y
@syms z
Iyy = ρ * integrate(f(x), (x, 100cm, 1m + 100cm), (y, 0, 2m), (z, 0, 3m))
Iyy_old = ρ * integrate(f(x), (x, 0, 1m), (y, 0, 2m), (z, 0, 3m))
Iyy_old + ρ * 1m * 2 * 3 * (100cm)^2 # WRONG BECAUSE STEINER THEOREM ONLY WRT CoM