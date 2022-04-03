import sympy as sp
import sympy.physics.units.quantities as sq
from sympy.physics.quantum.constants import hbar
from custom_libraries.stepper import *

PSI_FUNCTION = sp.Function( "psi" )
POTENTIAL_FUNCTION = sp.Function( "V" )
TOTAL_ENERGY_SYMBOL = sq.Symbol( 'E', nonzero = True, positive = True )
MASS_SYMBOL = sq.Quantity( 'm', positive = True, nonzero = True )
POSITION_SYMBOL = sp.Symbol( 'x', positive = True )

def time_independent_schroedinger_equation_1d( 
            psi = PSI_FUNCTION, 
            potential = POTENTIAL_FUNCTION, 
            total_energy = TOTAL_ENERGY_SYMBOL, 
            mass = MASS_SYMBOL, 
            reduced_planck_constant = hbar, 
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
    
    def default_non_uniform_length_potential_table(): 
        return { 
                    StairWell.NON_UNIFORM_STAIR_LENGTHS[ ii ] : StairWell.NON_UNIFORM_POTENTIALS[ ii ] 
                    for ii in range( len( StairWell.NON_UNIFORM_POTENTIALS ) ) 
            }
    
    def default_uniform_length_potential_table(): 
        return { 
                    StairWell.UNIFORM_STAIR_LENGTHS[ ii ] : StairWell.UNIFORM_POTENTIALS[ ii ] 
                    for ii in range( len( StairWell.UNIFORM_POTENTIALS ) ) 
            }

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

def zip_to_table( keys, values ): 
    assert len( keys ) == len( values )
    keys = list( keys )
    values = list( values )
    return { keys[ ii ] : values[ ii ] for ii in range( len( keys ) ) }

def table_to_pairs( table ): 
    keys = list( table.keys() )
    values = list( table.values() )
    return [ ( keys[ ii ], values[ ii ] ) for ii in range( len( keys ) ) ]

make_psi_parameter_constrained = lambda psi_function, region_potential_table, position : \
        [ psi_function( position < region ) for region in region_potential_table.keys() ]

make_numbered_psi_function = lambda psi_function, region_potential_table : [ \
        sp.Function( '\\' + str( psi_function ) + "_{" + str( ii ) + '}' ) \
        for ii in range( len( region_potential_table ) ) \
    ]

make_psi_numbered = lambda psi_function, region_potential_table, position : [ \
        psi( position ) for psi in make_numbered_psi_function( psi_function, region_potential_table ) \
    ]

def make_psi_numbered_parameter_constrained( psi_function, region_potential_table, position ): 
    psis = make_numbered_psi_function( psi_function, region_potential_table )
    regions = region_potential_table.keys()
    return [ psis[ ii ]( position < regions[ ii ] ) for ii in range( len( region_potential_table ) ) ]

class TimeIndependentSchrodingerConstantPotentials1D: 

    DEFAULT_CHECK_POINT_NAME_BASE = "TimeIndependentSchrodingerConstantPotentials1DCheckPoint"
    DEFAULT_CONSTANT_NAME_BASE = 'k'
    
    BEFORE_SOLVE_HARMONIC = "BeforeSolveHarmonic"
    SOLVED_HARMONIC = "SolveHarmonicSolved"
    
    def __init__( 
                self, 
                region_potential_table, 
                repeating = False, 
                inverse_repeating = False, 
                psi_function = PSI_FUNCTION, 
                potential_function = POTENTIAL_FUNCTION, 
                total_energy = TOTAL_ENERGY_SYMBOL, 
                mass = MASS_SYMBOL, 
                reduced_planck_constant = hbar, 
                position = POSITION_SYMBOL, 
                make_psis = make_psi_parameter_constrained, 
                check_point_name_base = DEFAULT_CHECK_POINT_NAME_BASE, 
                constant_name_base = DEFAULT_CONSTANT_NAME_BASE, 
                check_point_count= 0
            ): 
        self.region_potential_table = region_potential_table
        self.repeating = repeating 
        self.inverse_repeating = inverse_repeating 
        self.psi_function = psi_function 
        self.potential_function = potential_function 
        self.total_energy = total_energy 
        self.mass = mass 
        self.reduced_planck_constant = reduced_planck_constant 
        self.position = position
        self.make_psis = make_psis
        self.check_point_name_base = check_point_name_base
        self.constant_name_base = constant_name_base
        self.check_point_count = check_point_count
        self.psis = make_psis( self.psi_function, self.region_potential_table, self.position )
        region_psi_table = self.region_psi_table()
        
        self.vanilla_schrodinger_equation_1d = \
                Stepper( time_independent_schroedinger_equation_1d( 
                        self.psi_function, 
                        self.potential_function, 
                        self.total_energy, 
                        self.mass, 
                        self.reduced_planck_constant, 
                        self.position
                    ) )

        self.equations = [
                Stepper( self.vanilla_schrodinger_equation_1d.last_step() \
                        .subs( self.potential_function( self.position ), potential ) \
                        .replace( self.psi_function( self.position ), region_psi_table[ region ] ) ) \
                for region, potential in self.region_potentials()
            ]
    
    def region_potentials( self ): 
        return table_to_pairs( self.region_potential_table )
    
    def region_psi_table( self ): 
        return zip_to_table( self.region_potential_table.keys(), self.psis )
    
    def second_derivative( self, equation_index : int ): 
        return sp.Derivative( self.psis[ equation_index ], ( self.position, 2 ) )

    def _new_check_point_number( self ): 
        self.check_point_count += 1
        return str( self.check_point_count - 1 )
        
    def _harmonic_constant_name( self, equation_index, name_base = None ): 
        name_base = not_none_value( name_base, self.constant_name_base )
        return name_base + '_' + str( equation_index )
    
    def make_harmonic_constant( self, equation_index : int, name_base = None ): 
        equation = self.equations[ equation_index ]
        psi = self.psis[ equation_index ]
        before_check_point = self.check_point_name_base \
                + TimeIndependentSchrodingerConstantPotentials1D.BEFORE_SOLVE_HARMONIC \
                + self._new_check_point_number()
        equation.check_point( before_check_point )
        psi_solution = sp.solve( equation.last_step(), psi )
        psi_solution = psi_solution[ 0 ] if type( psi_solution ) is list else psi_solution
        #sp.pprint( psi_solution )
        equation.add_step( sp.Eq( psi, psi_solution ) )
        # For Stepper's / and /= are the same
        second_deriviative = self.second_derivative( equation_index )
        #equation.add_step( sp.Eq( 
        #        psi / second_deriviative, 
        #        psi_solution / self.second_derivative( equation_index ) 
        #    ) )
        #sp.pprint( equation.last_step() )
        equation / second_deriviative
        equation ** ( -1 )
        constant_name = self._harmonic_constant_name( equation_index, name_base )
        equation.right_to_constant( constant_name )
        equation.check_point( 
                self.check_point_name_base \
                    + TimeIndependentSchrodingerConstantPotentials1D.SOLVED_HARMONIC \
                    + self._new_check_point_number()
            )
        equation.restore_from_check_point( before_check_point, True )
        #print( dir( equation.symbols() ) )
        return equation.constant_symbols().symbol_by_string_name( constant_name )
