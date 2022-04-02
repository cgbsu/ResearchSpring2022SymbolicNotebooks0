import sympy as sp
import sympy.physics.units.quantities as sq
from sympy.physics.quantum.constants import hbar

PSI_FUNCTION = sp.Function( "psi" )
POTENTIAL_FUNCTION = sp.Function( "V" )
TOTAL_ENERGY_SYMBOL = sq.Symbol( 'E', nonzero = True, positive = True )
MASS_SYMBOL = sq.Quantity( 'm', positive = True, nonzero = True )
POSITION_SYMBOL = sp.Symbol( 'x', positive = True )

def time_independent_schroedinger_equation( 
            psi = PSI_FUNCTION, 
            potential = POTENTIAL_FUNCTION, 
            total_energy = TOTAL_ENERGY_SYMBOL, 
            mass = MASS_SYMBOL, 
            reduced_planck_constant = hbar, #REDUCED_PLANCK_CONSTANT_SYMBOL, 
            position = POSITION_SYMBOL 
        ): 
    return sp.Eq( 
            ( ( -( reduced_planck_constant ** 2 ) / ( 2 * mass ) ) \
                    * sp.Derivative( psi( position ), ( position, 2 ) ) )
                    + ( potential( position ) * psi( position ) ), 
            total_energy * psi( position ) 
        )

class ZeroPotential( sp.Function ): 
    @classmethod
    def eval( cls, position ):
        return 0

#https://docs.sympy.org/latest/modules/assumptions/index.html#querying
#https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cxxnodes

class Potential( sp.Function ): 
    DEFAULT_POTENTIAL = sp.Symbol( "V_0" )
    @classmethod
    def eval( cls, position, potential = DEFAULT_POTENTIAL ): 
        return potential

class TunnelPotential( sp.Function ): 
    DEFAULT_WELL_LENGTH = sp.Symbol( 'L' )
    DEFAULT_POTENTIAL = Potential.DEFAULT_POTENTIAL
    DEFAULT_START = 0
    @classmethod
    def eval( cls, position, length = DEFAULT_WELL_LENGTH, 
             start = DEFAULT_START, potential = DEFAULT_POTENTIAL ): 
        if position < start or position > sp.simplify( length + start ): 
            return ZeroPotential.eval( position )
        return Potential.eval( position, potential )

class StairWell( sp.Function ): 
    UNIFORM_LENGTH_SYMBOL = sp.Symbol( 'L' )
    UNIFORM_POTENTIAL_SYMBOL = sp.Symbol( 'V' )
    UNIFORM_STAIR_LENGTHS = ( 
            UNIFORM_LENGTH_SYMBOL, 
            UNIFORM_LENGTH_SYMBOL, 
            UNIFORM_LENGTH_SYMBOL 
        )
    UNIFORM_POTENTIALS = ( 
            UNIFORM_POTENTIAL_SYMBOL, 
            2 * UNIFORM_POTENTIAL_SYMBOL / 3, 
            UNIFORM_POTENTIAL_SYMBOL / 3
        )
    NON_UNIFORM_STAIR_LENGTHS = sp.symbols( "L_0 L_1 L_2" )
    NON_UNIFORM_POTENTIALS = sp.symbols( "V_0 V_1 V_2" )
    #NON_UNIFORM_ASSUMPTIONS = sp.Q.gt( NON_UNIFORM_POTENTIALS[ 0 ], NON_UNIFORM_POTENTIALS[ 1 ] ) \
    #        & sp.Q.gt( NON_UNIFORM_POTENTIALS[ 1 ], NON_UNIFORM_POTENTIALS[ 2 ] ) \
    #        & sp.Q.gt( NON_UNIFORM_STAIR_LENGTHS[ 0 ] + NON_UNIFORM_STAIR_LENGTHS[ 1 ], NON_UNIFORM_STAIR_LENGTHS[ 0 ] ) \
    #        & sp.Q.gt( NON_UNIFORM_STAIR_LENGTHS[ 2 ] + NON_UNIFORM_STAIR_LENGTHS[ 1 ], NON_UNIFORM_STAIR_LENGTHS[ 2 ] )  
    #DEFAULT_ASSUMPTIONS = NON_UNIFORM_ASSUMPTIONS \
    #        & sp.Q.positive( UNIFORM_LENGTH_SYMBOL ) \
    #        & sp.Q.positive( UNIFORM_POTENTIAL_SYMBOL )
    DEFAULT_START = 0
    
    @classmethod
    def eval( 
                cls, 
                position, 
                start = DEFAULT_START, 
                potentials = NON_UNIFORM_POTENTIALS, 
                lengths = NON_UNIFORM_STAIR_LENGTHS#, 
                #assumptions = DEFAULT_ASSUMPTIONS 
            ): 
        position = position + start
        if position < lengths[ 0 ]: 
            return potentials[ 0 ]
        elif position < ( lengths[ 0 ] + lengths[ 1 ] ): 
            return potentials[ 1 ]
        elif position < ( lengths[ 0 ] + lengths[ 1 ] + lengths[ 2 ] ): 
            return potentials[ 2 ]
