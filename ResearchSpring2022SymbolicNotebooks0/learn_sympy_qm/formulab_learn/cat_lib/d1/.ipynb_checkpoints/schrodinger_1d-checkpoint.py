import sympy as sp
import sympy.physics.units.quantities as sq
#https://sympygamma.com/
PSI_FUNCTION = sp.Function( "psi" )
POTENTIAL_FUNCTION = sp.Function( "V" )
TOTAL_ENERGY_SYMBOL = sq.Symbol( 'E', nonzero = True, positive = True )
MASS_SYMBOL = sq.Quantity( 'm', positive = True, nonzero = True )
REDUCED_PLANCK_CONSTANT_SYMBOL = sq.Quantity( "hbar" )
POSITION_SYMBOL = sp.Symbol( 'x', positive = True )

def time_independent_schroedinger_equation( 
            psi = PSI_FUNCTION, 
            potential = POTENTIAL_FUNCTION, 
            total_energy = TOTAL_ENERGY_SYMBOL, 
            mass = MASS_SYMBOL, 
            reduced_planck_constant = REDUCED_PLANCK_CONSTANT_SYMBOL, 
            position = POSITION_SYMBOL 
        ): 
    return sp.Eq( 
            ( ( -( reduced_planck_constant ** 2 ) / ( 2 * mass ) ) \
                    * sp.Derivative( psi( position ), ( position, 2 ) ) )
                    + ( potential( position ) * psi( position ) ), 
            total_energy * psi( position ) 
        )
