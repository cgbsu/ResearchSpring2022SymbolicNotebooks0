// function psi
//     input Real x;
//     parameter Real C_0;
//     parameter Real C_1;
//     parameter Real k;
//     output Real psi_x;
// algorithm
//     psi_x := ( C_0 * exp( k * x ) ) + ( C_1 * exp( -k * x ) );
// end psi;

model InfiniteWell
    input Real x(start = 0.0);
    parameter Real m;
    parameter Real E;
    parameter Real V_0;
    parameter Real hbar;
    parameter Real L;
    output Real k;
    output Real C_0;
    output Real C_1;
    output Real psi_x;
protected
    Real V_x;
equation
    k * k = (2 * (V_x - E) * m) / (hbar * hbar);
    der(der(psi_x)) + k * psi_x = 0;
    der(psi_x * psi_x) = der(1);
    psi_x = ( C_0 * exp( k * x ) ) + ( C_1 * exp( -k * x ) );
    (C_0 * exp(0)) + (C_1 * exp(0)) - (C_0 * exp(k * L)) + (C_1 * exp(-k * L)) = 0;
initial equation
    if 0 < x or x < L then
        V_x = V_0;
    else
        V_x = 0;
    end if;
end InfiniteWell;
