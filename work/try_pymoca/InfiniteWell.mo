function psi
    input Real x;
    parameter Real C_0;
    parameter Real C_1;
    parameter Real k;
    output Real psi_x;
algorithm
    psi_x := ( C_0 * exp( k * x ) ) + ( C_1 * exp( -k * x ) );
end psi;

model InfiniteWell
    input Real from;
    input Real to;
    input Real x;
    parameter Real m;
    parameter Real E;
    parameter Real V_0;
    parameter Real hbar;
    parameter Real L;
    Real V_x;
    output Real prob_amp;
    output Real k;
equation
    k * k = (2 * (V_x - E) * m) / (hbar * hbar);
    der(der(psi(x))) + k * psi(x) = 0;
    intg(psi(x) * psi(x)) = 1;
//    ((to - from) * (to - from)) = 1;
//    der(prob_amp) = 1;
initial equation
    if 0 < x or x < L then
        V_x = V_0;
    else
        V_x = 0;
    end if;
    
end InfiniteWell;
