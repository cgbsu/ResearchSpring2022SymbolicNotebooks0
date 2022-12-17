// function psi
//     input Real x;
//     parameter Real C_0;
//     parameter Real C_1;
//     parameter Real k;
//     output Real psi_x;
// algorithm
//     psi_x := ( C_0 * exp( k * x ) ) + ( C_1 * exp( -k * x ) );
// end psi;

model HarmonicConstant
    input Real V;
    parameter Real m;
    parameter Real E;
    parameter Real hbar;
    output Real k;
equation
    k = 2 * m * (E - V) / (hbar^2);
end HarmonicConstant;

model InfiniteWell_1
    input Real x(start = 0.0);
    parameter HarmonicConstant k;
    parameter Real V_0;
    parameter Real L;
    output Real psi_x;
equation
    der(der(x)) = -k.k * x;
    psi_x = der(der(x)) + k.k * x;
    // der(psi_x ^ 2) = der(1);
initial equation
    if 0 < x or L < x then
        k.V = V_0;
        // psi_x.psi_x = 0;
    else
        k.V = 0;
    end if;
end InfiniteWell_1;
