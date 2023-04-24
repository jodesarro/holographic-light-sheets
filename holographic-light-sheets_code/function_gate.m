function porta = function_gate (Z1,Z2)
syms z
sympref('HeavisideAtOrigin', 0);
if(Z2>Z1)
    porta(z) = heaviside(z-Z1)+heaviside(-z+Z2)-1;
else
    porta(z) = 0.*z;
end