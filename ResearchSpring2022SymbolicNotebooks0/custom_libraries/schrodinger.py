import sympy as sp
import sympy.physics.units.quantities as sq
from sympy.physics.quantum.constants import hbar
from custom_libraries.stepper import *
from custom_libraries.boundries import *
import copy

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

DEFAULT_NORMALIZATION_VALUE = 1
INDEFINITE_NORMALIZATION_INTEGRAL = None
DEFAULT_CONJUGATE_NOT_SQUARED_ABSOLUTE_VALUE = True

def simple_wave_function_normalization( 
            from_, to_or_indefinite, 
            psi = PSI_FUNCTION, 
            position = POSITION_SYMBOL, 
            normalization_value = DEFAULT_NORMALIZATION_VALUE, 
            conjugate_not_squared_absolute_value 
                    = DEFAULT_CONJUGATE_NOT_SQUARED_ABSOLUTE_VALUE 
        ): 
    integral_function = psi.func( position ) * sp.conjugate( psi.func( position ) ) \
            if conjugate_not_squared_absolute_value \
            else sp.Abs( psi.func( position ) ) ** 2
    if to_or_indefinite: 
        return sp.Eq( 
                sp.Integral( 
                        integral_function, 
                        ( position, from_, to_or_indefinite ) 
                    ),
                normalization_value 
            )
    else:
        return sp.Eq( 
                sp.Integral( integral_function ),
                normalization_value 
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
    UNIFORM_LENGTH_SYMBOL = sp.Symbol( 'L', real = True, finite = True, nonzero = True )
    UNIFORM_POTENTIAL_SYMBOL = sp.Symbol( 'V', real = True, finite = True, nonzero = True )
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
    NON_UNIFORM_STAIR_LENGTHS = sp.symbols( "L_0 L_1 L_2", real = True, finite = True, nonzero = True )
    NON_UNIFORM_POTENTIALS = sp.symbols( "V_0 V_1 V_2", real = True, finite = True, nonzero = True )
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

non_list_to_list = lambda canidate : canidate if type( canidate ) is list else [ canidate ]
is_type_len_gt_0 = lambda type_, canidate : len( canidate ) > 0 if type( canidate ) is type_ else False


to_functions = lambda functions_with_parameters : tuple( function for function in functions_with_parameters )
default_boundry_constant_name_base = lambda _ : "B"

standard_harmonic_assumptions = lambda equation : { 'finite' : True, 'nonzero' : True }

region_warning = lambda region : f"""Error: {region} is a symbol that does not have 
                    the property of being in the set of real numbers ("Real"), this can lead to normalizations 
                    getting quite complicated. You may make `{region}` Real by setting the real property to 
                    'True' in Symbol's constructor, e.g {region} = sympy.Symbol( "L_{0}", real = True ). If 
                    you know what your doing do not want to assume that your {region} is real, please set 
                    `ensure_lengths_are_real` to False"""

def warn_about_region( region_key, region ): 
    assert region_key.assumptions0[ 'real' ] == True if type( region_key ) is sp.Symbol else True, region_warning( region )

class TimeIndependentSchrodingerConstantPotentials1D( Symbols ): 

    DEFAULT_CHECK_POINT_NAME_BASE = "TimeIndependentSchrodingerConstantPotentials1DCheckPoint"
    DEFAULT_CONSTANT_NAME_BASE = 'k'
    
    CHECK_POINT_BEFORE_SOLVE_HARMONIC_CONSTANT = "BeforeSolveHarmonicConstant"
    CHECK_POINT_SOLVED_HARMONIC_CONSTANT = "SolveHarmonicConstantSolved"
    CHECK_POINT_SUBSTITUTE_HARMONIC_CONSTANT = "SubstitutingHarmonicConstant"
    CHECK_POINT_HARMONIC_SOLUTION_TO_CANONOCAL_FORM = "HarmonicSolutionToCanonicalForm"
    CHECK_POINT_SOLVERS_ODE_DSOLVE = "SolveODEsWithSolversOdeDsolve"
    CHECK_POINT_BOUNDRY_TO_CONSTANT_SUBSTITUTION = "BoundryToConstantSubstitution"
    CHECK_POINT_SUBSTUTE_WAVE_FUNCTIONS_INTO_NORMALIATIONS = "SubstituteWaveFunctionsIntoNormalizations"
    CHECK_POINT_INTEGRATE_NORMALIZATION = "IntegrateNormalization"

    BOUNDRY_CONTINUITY_CONDITIONS = "ContinuityConditions"
    BOUNDRY_REPEATING_POTENTIALS_CONDITION = "RepeatingPotentialsCondition"
    BOUNDRY_ALL = "LastUpdatedAllBoundryConditions"
    BOUNDRY_CONSTANT_SUBSTITUTIONS = "ConstantSubstitutions"
    BOUNDRY_ZERO_CONDITIONS = "BoundryZeroConditions"
    BOUNDRY_BOUNDRY_SIMPLIFICATION_TABLE = "BoundrySimplificationTable"
    BOUNDRY_BOUNDRY_CONSTANT_TABLE = "BoundryConstantTable"
    
    DEFAULT_NORMALIZATION_CONSTANT_BASE_NAME = 'N'
    DEFAULT_TOTAL_NORMALIZATION_VALUE = 1 
    
    COMMIT_CHECK_POINT_PREFIX_BEFORE = "Before"
    COMMIT_CHECK_POINT_PREFIX_POST = "Post"
    
    
    def __init__( 
                self, 
                region_potential_table, 
                region_end = None, 
                region_start = None, 
                initial_boundries = None, 
                repeating = False, 
                inverse_repeating = False, 
                psi_function = PSI_FUNCTION, 
                potential_function = POTENTIAL_FUNCTION, 
                total_energy = TOTAL_ENERGY_SYMBOL, 
                mass = MASS_SYMBOL, 
                reduced_planck_constant = hbar, 
                position = POSITION_SYMBOL, 
                make_psis = make_psi_numbered, 
                check_point_name_base = DEFAULT_CHECK_POINT_NAME_BASE, 
                constant_name_base = DEFAULT_CONSTANT_NAME_BASE, 
                check_point_count = 0, 
                psi_parameter_based = None, 
                create_normalization = simple_wave_function_normalization, 
                normalization_conjugate_not_squared_absolute_value = DEFAULT_CONJUGATE_NOT_SQUARED_ABSOLUTE_VALUE, 
                normaliation_constant_base_name = DEFAULT_NORMALIZATION_CONSTANT_BASE_NAME, 
                create_constant_table = equations_to_constant_table, 
                total_normalization_value = DEFAULT_TOTAL_NORMALIZATION_VALUE, 
                create_boundry_constant_name_base = default_boundry_constant_name_base, 
                as_distances = False, 
                harmonics_assumptions = standard_harmonic_assumptions, 
                ensure_lengths_are_real = True, 
            ): 
        super().__init__()
        self.region_potential_table = region_potential_table
        self.region_start = not_none_value( region_start, 0 )
        self.region_end = not_none_value( region_end, tuple( region_potential_table.keys() )[ -1 ] )
        if ensure_lengths_are_real: 
            warn_about_region( self.region_start, "region_start" )
            warn_about_region( self.region_end, "region_end" )
            for region in self.region_potential_table: 
                warn_about_region( self.region_start, f"region_potential_table_key({region})" )
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
        self.psi_parameter_based = not_none_value( psi_parameter_based, self.make_psis == make_psi_parameter_constrained )
        self.create_constant_table = create_constant_table
        self.psis = self.make_psis( self.psi_function, self.region_potential_table, self.position )
        self.boundries = Boundries( initial_boundries )
        self.create_normalization = create_normalization 
        self.normaliation_constant_base_name = normaliation_constant_base_name
        self.normalization_conjugate_not_squared_absolute_value = normalization_conjugate_not_squared_absolute_value
        self.total_normalization_value = total_normalization_value
        self.create_boundry_constant_name_base = create_boundry_constant_name_base
        self.boundry_constant_symbols = []
        self.distance_constant_symbols = []
        self.normalizations = []
        self.constant_solutions = {}
        self.harmonics_assumptions = harmonics_assumptions
        self.as_distances = as_distances
        self._create_schrodinger_equations()
        self._create_normalizations()
        self.harmonic_constants = self._make_harmonic_constants()
        self._impose_continuity_conditions()
        if repeating == True: 
            self.impose_repeating_potentials_condition()
    
    def _create_schrodinger_equations( self ): 
        region_psi_table = self.region_psi_table()
        self.vanilla_schrodinger_equation_1d = \
                time_independent_schroedinger_equation_1d( 
                        self.psi_function, 
                        self.potential_function, 
                        self.total_energy, 
                        self.mass, 
                        self.reduced_planck_constant, 
                        self.position
                    )
        self.equations = [
                Stepper( self.vanilla_schrodinger_equation_1d \
                        .subs( self.potential_function( self.position ), potential ) \
                        .replace( self.psi_function( self.position ), region_psi_table[ region ] ) ) \
                for region, potential in self.region_potentials()
            ]
        for psi in self.psis: 
            self.add_symbol( psi )
    
    def _create_normalizations( self ): 
        previous_from = self.region_start
        regions = self.regions()
        self.normalization_symbols = []
        #Yes, specifying the parameters makes this less customizable
        for ii in range( len( self.equations ) ): 
            self.normalization_symbols.append( 
                    sp.Symbol( self.normaliation_constant_base_name + "_" + str( ii ) )#+ '}' )
                )
            self.add_symbol( self.normalization_symbols[ ii ] )
            self.normalizations.append( self.create_normalization( 
                        0 if self.as_distances else previous_from, regions[ ii ], 
                        psi = self.psis[ ii ], 
                        position = self.position, 
                        normalization_value = self.normalization_symbols[ ii ], 
                        conjugate_not_squared_absolute_value = \
                                self.normalization_conjugate_not_squared_absolute_value
                    ).reversed
                )
            previous_from = regions[ ii ]
        self.total_normalization = Stepper( 
                sp.Eq( 
                        self.total_normalization_value, 
                        sum( [ symbol for symbol in self.normalization_symbols ] ) 
                    ), 
                constants = self.create_constant_table( self.normalizations )
            )
        self.normalizations = [ self.total_normalization.constants[ constant ] 
                               for constant in self.total_normalization.constants ]
            
    def regions( self ): 
        return tuple( self.region_potential_table.keys() )
    
    def potentials( self ): 
        return ( self.region_potential_table[ region ] for region in self.region_potential_table )
    
    def region_potentials( self ): 
        return table_to_pairs( self.region_potential_table )
    
    def region_psi_table( self ): 
        return zip_to_table( self.region_potential_table.keys(), self.psis )
    
    def second_derivative( self, equation_index : int ): 
        return sp.Derivative( self.psis[ equation_index ], ( self.position, 2 ) )

    def _new_check_point_number( self ): 
        self.check_point_count += 1
        return str( self.check_point_count - 1 )
    
    def _equation_new_check_point( self, equation, check_point, check_point_number = None, pass_check_point_name = False ): 
        check_point_number = not_none_value( check_point_number, self._new_check_point_number() )
        return equation.check_point( 
                self.check_point_name_base \
                        + check_point \
                        + str( check_point_number ), 
                pass_check_point_name = pass_check_point_name 
            )
    
    def _equations_new_check_point( self, check_point, check_point_number = None ): 
        check_point_number = not_none_value( check_point_number, self._new_check_point_number() )
        check_point_names = [ [], [] ]
        new_check_point = lambda equation, index : check_point_names[ index ].append( 
                        self._equation_new_check_point( 
                                equation, 
                                check_point, 
                                check_point_number, 
                                pass_check_point_name = True 
                            )
                    )
        for ii in range( len( self.equations ) ): 
            new_check_point( self.equations[ ii ], 0 )
            new_check_point( self.normalizations[ ii ], 1 )
        return check_point_names
    
    def _harmonic_constant_name( self, equation_index, name_base = None ): 
        name_base = not_none_value( name_base, self.constant_name_base )
        return name_base + '_' + str( equation_index )
    
    def _make_harmonic_constant( self, equation_index : int, name_base = None, restore_before_checkpoint = False ): 
        equation = self.equations[ equation_index ]
        psi = self.psis[ equation_index ]
        self._equations_new_check_point( 
                TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_BEFORE_SOLVE_HARMONIC_CONSTANT 
            )
        psi_solution = sp.solve( equation.last_step(), psi )
        psi_solution = psi_solution[ 0 ] if type( psi_solution ) is list else psi_solution
        equation.add_step( sp.Eq( psi, psi_solution ) )
        second_deriviative = self.second_derivative( equation_index )
        # For Stepper's / and ** are the same as /= **=
        equation /= second_deriviative
        equation **= ( -1 )
        equation.root( 2 )
        constant_name = self._harmonic_constant_name( equation_index, name_base )
        equation.right_to_constant( constant_name, assumptions = self.harmonics_assumptions( equation_index ) )
        self._equations_new_check_point( 
                TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_SOLVED_HARMONIC_CONSTANT 
            )
        equation.substitute_constant( equation.constants_as_symbols().symbol_by_string_name( constant_name ) )
        self._equations_new_check_point( 
                TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_SUBSTITUTE_HARMONIC_CONSTANT 
            )
        equation **= 2
        equation *= psi
        equation -= equation.right()
        self._equations_new_check_point( 
                TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_HARMONIC_SOLUTION_TO_CANONOCAL_FORM 
            )
        if restore_before_checkpoint: 
            equation.restore_from_check_point( before_check_point, True )
        return equation.constant_symbols().symbol_by_string_name( constant_name )
    
    def _make_harmonic_constants( self, name_base = None ): 
        return [ self._make_harmonic_constant( ii, name_base ) for ii in range( len( self.region_potential_table ) ) ]
    
    def _impose_continuity_conditions( self, parameter_based = None ): 
        """Run automatically in __init__ constructor, standard requirment of Quantum mechanichs where the 
        wave function must be continuous, realized here by the symbolic representations of what is 
        (hopefully) functionally equivelent too: 
        let there be a wave functions psi indexed by i, with corresponding distances D indexed by i
        let \psi_{i}( L ) = \psi_{i + 1}( L ) where L = \sum_0^i{D_i}
        This will be imposed for every wave function in the set that has a wave function following it
        this is non-circular, see `impose_repeating_potentials_condition` for making 
        a repeating condition
        
        return: name of "boundry set" (str) added and boundry set (dict)
        rtype: tuple( str, dict )
        """
        assert not not_none_value( self.psi_parameter_based, parameter_based ), """
        Parameter based continuity condition not yet supported! Try using numbered instead (default)
        TODO: Requires using infentesmials
        """
        regions = tuple( self.region_psi_table().keys() )
        assert len( self.psis ) == len( regions )
        length = len( self.psis )
        return self.boundries.add_boundries( 
                TimeIndependentSchrodingerConstantPotentials1D.BOUNDRY_CONTINUITY_CONDITIONS, {
                        self.psis[ ii ].func( regions[ ii ] ) : self.psis[ ii + 1 ].func( 0 if self.as_distances else regions[ ii ] )
                                for ii in range( length ) if ( ii + 1 ) < length
            } )

    def impose_zero_conditions_to_ends( self ): 
        return tuple( [ 
                self.impose_zero_condition_to_boundry( sp.sympify( 0 ), self.psis[ 0 ] ), 
                self.impose_zero_condition_to_boundry( self.regions()[ -1 ], self.psis[ -1 ] ) 
            ] )
    
    def impose_constant_to_ends( self, constant, second_constant = None ): 
        return tuple( [ 
                self.impose_constant_to_boundry( sp.sympify( 0 ), self.psis[ 0 ], constant ), 
                self.impose_constant_to_boundry( self.regions()[ -1 ], 
                        self.psis[ -1 ], sp.sympify( not_none_value( second_constant, constant ) )
                   ) 
            ] )
    
    def impose_zero_condition_to_boundry( self, region, wave_function ): 
        return self.impose_constant_to_boundry( region, wave_function, 0 )

    def impose_constant_to_boundry( self, region, wave_function, constant ): 
        return self.boundries.add_boundries( 
                TimeIndependentSchrodingerConstantPotentials1D.BOUNDRY_ZERO_CONDITIONS, {
                        wave_function.func( region ) : constant 
            }, automatically_append = True )
    
    def update_harmonic_constants( self, name_base = None ): 
        self.harmonic_constants = self._make_harmonic_constants( name_base )
        return self.harmonic_constants
    
    def psis_to_functions( self ): 
        return to_functions( self.psis )
    
    def normalization_steppers( self ): 
        return [ self.total_normalization.constants[ constant ] 
                for constant in self.total_normalization.constants ]
    
    def normalization_constants( self ): 
        return self.total_normalization.constants
    
    def impose_repeating_potentials_condition( self, start = None, end = None ): 
        """Impose the assumption that the first wave function 
        of the first (zeroth/0th) region at `start` will 
        be equal to the last ("nth"/"n'th") wave function 
        at `end`, `start` defaults to `self.region_start`, 
        and `end` to `self.region_end`
        
        Arguments: 
        
        start: The parameter to the first/0th/zeroth wave function, 
                physical location of the beggining of the first potential 
                region, defaults to `self.region_start`, the zeroth 
                wave function of this parameter will be set to the nth/last 
                wave function of `end`
        end: The parameter to the last/nth wave function and the physical 
                location of the end of the last potential region, defaults 
                to `self.region_end`. the last/nth wave function of this 
                parameter will be set to the first/0th/zeroth wave function 
                of `end`
        return: name of "boundry set" (str) added and boundry set (dict)
        rtype: tuple( str, dict )
        """
        start = not_none_value( start, self.region_start )
        end = not_none_value( end, self.region_end )
        psis = self.psis_to_functions() 
        return self.boundries.add_boundries( 
                TimeIndependentSchrodingerConstantPotentials1D.BOUNDRY_REPEATING_POTENTIALS_CONDITION, { 
                        psis[ 0 ].func( start ) : psis[ -1 ].func( end ) 
            } )
    
    def solve_odes( self, input_into_normalizations = False ):
        solutions = []
        for ii in range( len( self.equations ) ): 
            current_function = self.psis[ ii ]
            boundries = self.boundries.boundries_with( current_function.func )
            position = self.equations[ ii ].symbols().symbol_by_string_name( self.position )
            solutions.append( sp.solvers.ode.dsolve( 
                    self.equations[ ii ].last_step(), 
                    current_function.func( position ), 
                    ics = boundries
                ) )
            self.equations[ ii ].add_step( solutions[ ii ] )
        self._equations_new_check_point( 
                TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_SOLVERS_ODE_DSOLVE 
            )
        if input_into_normalizations == True: 
            self.substitute_wave_functions_into_normalizations()
        return solutions
    
    def set_imposed_total_engineered_normalization_values( self, normalization_values : list ): 
        """This is the "total" engineered normalization because it MUST be applied to ALL 
        normalizations for ALL wave functions (this is for simplicity of implentation, 
        possibly more comprehensive options in the future).
        
        Arguments: 
        
        normalization_values: a list of values from 0 - self.total_normalization_value (inclusive), where 
                the sum of the list MUST be self.total_normalization_value (self.total_normalization_value representing 100%). 
                These are the probabilities for each region that the particle will be found in the region. This does not shape the whole 
                wave function, but it can be shaped by creating more regions."""
        assert sum( normalization_values ) == self.total_normalization_value, """
                Total probability of normalization values not equal to self.total_normalization_value!
                See help( TimeIndependentSchrodingerConstantPotentials1D.set_imposed_total_engineered_normalization_values ) for details"""
        assert len( normalization_values ) == len( self.psis ), """The number of wave functions does not match 
                the number of normalization values! The method has a 1-to-1 correlation between the two, see 
                See help( TimeIndependentSchrodingerConstantPotentials1D.set_imposed_total_engineered_normalization_values ) for details"""
        for ii in range( len( normalization_values ) ): 
            self.normalizations[ ii ].operate( lambda step : 
                    step.subs( { self.normalization_symbols[ ii ] : normalizations_values[ ii ] } ) 
                )
        return self.normalizations

    def boundries_in_expression_to_constants( 
                self, 
                equation = None, 
                automatically_append = True,  
                simplifcication_table_name = BOUNDRY_BOUNDRY_SIMPLIFICATION_TABLE, 
                constant_table_name = BOUNDRY_BOUNDRY_CONSTANT_TABLE, 
                before_prefix = COMMIT_CHECK_POINT_PREFIX_BEFORE, 
                after_prefix = COMMIT_CHECK_POINT_PREFIX_POST, 
                exclude_substituting_numbers = True 
            ):
        """Dont change the dafaults after `automatically_append`, they are just there to make the code cleaner really"""
        original_boundry_set_name, boundry_simplification_list = self.boundries.update_all_boundry_conditions()
        boundry_constant_base_name = self.create_boundry_constant_name_base( original_boundry_set_name )
        self.boundries.commit( simplifcication_table_name, before_prefix )
        self.boundries.add_boundries( 
                simplifcication_table_name, 
                boundry_simplification_list, 
                automatically_append = automatically_append
            )
        boundry_number = 0
        constant_substitution_table = getattr( self.boundries, constant_table_name ) \
                if hasattr( self.boundries, constant_table_name ) else {}
        for key in self.boundries.boundries[ simplifcication_table_name ]: 
            new_substitution_key = boundry_simplification_list[ key ]
            if not new_substitution_key in constant_substitution_table: 
                constant_substitution_table[ new_substitution_key ] = \
                        sp.Symbol( boundry_constant_base_name + "_{" + str( boundry_number ) + '}' )
                boundry_number += 1
        self.boundries.commit( constant_table_name, before_prefix )
        check_points = {}
        check_points[ before_prefix ] = self._equations_new_check_point( 
                before_prefix + TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_BOUNDRY_TO_CONSTANT_SUBSTITUTION
            )
        self.boundries.add_boundries( 
                constant_table_name, 
                constant_substitution_table, 
                automatically_append = automatically_append 
            )
        self.boundries.commit( constant_table_name, after_prefix )
        self.boundry_constant_symbols = []
        for ii in range( len( self.equations ) ): 
            self.equations[ ii ] \
                    .operate( lambda step : step.subs( boundry_simplification_list ), chain = True ).last_step()
            self.normalizations[ ii ] \
                    .operate( lambda step : step.subs( boundry_simplification_list ), chain = True )
        for key in constant_substitution_table: 
            constant = constant_substitution_table[ key ]
            self.boundry_constant_symbols.append( constant )
            for ii in range( len( self.equations ) ): 
                is_number = type( key ) is int \
                        or type( key ) is float \
                        or ( key.is_Number if hasattr( key, 'is_Number' ) else False )                
                if ( is_number == False ) if exclude_substituting_numbers else True: 
                    try: 
                        self.equations[ ii ].replace_with_constant( key, constant )
                    except Exception: 
                        pass
                    try: 
                        self.normalizations[ ii ].replace_with_constant( key, constant )
                    except Exception: 
                        pass
        check_points[ after_prefix ] = self._equations_new_check_point( 
                after_prefix + TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_BOUNDRY_TO_CONSTANT_SUBSTITUTION
            )
        return self.equations, self.normalizations, check_points
    
    def make_substitution_solution( 
                self, 
                position = None, 
                other_key = None, 
                before_prefix = COMMIT_CHECK_POINT_PREFIX_BEFORE, 
                after_prefix = COMMIT_CHECK_POINT_PREFIX_POST 
            ): 
        substitute_regions = position == None
        to_replace = not_none_value( other_key, self.position )
        substitutions = [ [], [] ]
        def substitute( equation ): 
                equation.operate( lambda step : step.subs( { to_replace : position } ) )
        def to_solution( equation ): 
                print( type( equation.last_step() ) )
                if type( equation.last_step() ) is sp.Eq: 
                    return equation.append_solutions_to_sets( 
                            solve_for = str( { to_replace : position } ), 
                            solutions = [ equation.last_step() ], 
                            automatically_make_new_solution_sets = True
                        )
                return []
        check_points = {}
        check_points[ before_prefix ] = self._equations_new_check_point( 
                before_prefix + "Sub" + str( to_replace ) + ':' + str( position )
            )
        for ii in range( len( self.equations ) ): 
            position = self.regions()[ ii ] if substitute_regions else position
            substitute( self.equations[ ii ] )
            substitute( self.normalizations[ ii ] )
        self.boundries_in_expression_to_constants()
        check_points[ after_prefix ] = self._equations_new_check_point( 
                after_prefix + "Sub" + str( to_replace ) + ':' + str( position )
            )
        for ii in range( len( self.equations ) ): 
            substitutions[ 0 ] += to_solution( self.equations[ ii ] )
            self.equations[ ii ].restore_from_check_point( check_points[ before_prefix ][ 0 ][ ii ] )
            substitutions[ 1 ] += to_solution( self.normalizations[ ii ] )
            self.normalizations[ ii ].restore_from_check_point( check_points[ before_prefix ][ 1 ][ ii ] )
        return substitutions
    
    def substitute_wave_functions_into_normalizations( 
                self, 
                before_prefix = COMMIT_CHECK_POINT_PREFIX_BEFORE, 
                after_prefix = COMMIT_CHECK_POINT_PREFIX_POST
            ): 
        self._equations_new_check_point( 
                before_prefix + TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_SUBSTUTE_WAVE_FUNCTIONS_INTO_NORMALIATIONS 
            )
        subtitution_table = { equation.last_step().lhs : equation.last_step().rhs for equation in self.equations }
        for normalization in self.normalizations: 
            normalization.operate( lambda step : step.subs( subtitution_table ) )
        self._equations_new_check_point( 
                after_prefix + TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_SUBSTUTE_WAVE_FUNCTIONS_INTO_NORMALIATIONS 
            )
        return self.normalizations
    
    # Probably can solve for constants from solving the differential equation later as well. #
    def solve_boundry_constants_from_equation( 
                self, 
                equation, 
                assume_2_solutions_squared_absolute_value = False, 
                transform = key_value_to_stepper, 
                auto_integrate = True # If you dont do this, 
                                      # it can take a long time to 
                                      # solve or even crash this is tested when 
                                      # lengths are real/finite and harmonic constants 
                                      # is real. 
            ): 
        solved_for = {}
        self.solved_boundries = []
        def debug( attempt, sols ): 
            display( attempt )
            if sols != None: 
                print( "I have sols!" )
            else: 
                print( "NO SOLS!!" )
            if type( sols ) is dict: 
                for key in sols: 
                    display( key )
                    display( sols[ key ] )
            elif type( sols ) is list: 
                for sol in sols: 
                    
                    display( sol )
            else: 
                display( sols )
        
        self._equation_new_check_point( equation, TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_INTEGRATE_NORMALIZATION )
        if self.normalizations.index( equation ) >= 0 and auto_integrate == True: 
            equation.operate( lambda step : step.expand().doit() )
        for boundry in self.boundry_constant_symbols: 
            if equation.last_step().has( boundry ): 
                boundry_solution = []
                print( "Has ", str( boundry ) )
                solutions = sp.solve( equation.last_step(), boundry )
                print( "Solutions: ", solutions )
                debug( boundry, solutions )
                if solutions == None: 
                    print( "FAILED solve for ", str( boundry ) )
                    continue
                else: 
                    print( "Solved for ", str( boundry ) )
                if is_type_len_gt_0( list, solutions ): 
                    boundry_solution += self.add_constant_solutions( 
                            equation, 
                            boundry, 
                            solutions, 
                            assume_2_solutions_squared_absolute_value, 
                            transform 
                        )
                elif is_type_len_gt_0( dict, solutions ): 
                    boundry_solution += non_list_to_list( enter_lists_dict_of_list( 
                            self.constant_solutions, 
                            solutions, 
                            transform = transform
                        ) )
                    equation.append_solutions_to_sets( 
                                solve_for = boundry, 
                                solutions = solutions, 
                                automatically_make_new_solution_sets = True, 
                                transform = transform 
                            )
                else: 
                    boundry_solution += non_list_to_list( enter_dict_of_list( 
                            self.constant_solutions, 
                            boundry, 
                            solutions, 
                            transform 
                        ) )
                    equation.append_solutions_to_sets( 
                                solve_for = boundry, 
                                solutions = solutions, 
                                automatically_make_new_solution_sets = True, 
                                transform = transform 
                            )
                enter_dict_of_list( solved_for, boundry, boundry_solution )
            return solved_for
                    
        #self.boundries.boundries[ TimeIndependentSchrodingerConstantPotentials1D.BOUNDRY_BOUNDRY_CONSTANT_TABLE ]

    def add_constant_solutions( 
                    self, 
                    equation, 
                    constant_symbol, 
                    solutions, 
                    assume_2_solutions_squared_absolute_value = False, 
                    transform = key_value_to_stepper
                ): 
            entered = non_list_to_list( enter_dict_of_list( 
                    self.constant_solutions, 
                    constant_symbol, 
                    solutions, 
                    transform = transform 
                ) )
            equation.append_solutions_to_sets( 
                        solve_for = constant_symbol, 
                        solutions = solutions, 
                        automatically_make_new_solution_sets = True, 
                        transform = transform 
                    )
            if len( solutions ) == 2 and assume_2_solutions_squared_absolute_value == True: 
                last_step = lambda step : step.last_step() if type( step ) is Stepper else step
                boundry_squared = sp.Abs( boundry ) ** 2
                entered.append( boundry_squared )
                boundry_squared_solution = ( last_step( solution[ 0 ] ) * last_step( solution[ 1 ] ) ).simplify() 
                equation.append_solutions_to_sets( 
                            solve_for = boundry_squared, 
                            solutions = boundry_squared_solution, 
                            automatically_make_new_solution_sets = True, 
                            transform = transform 
                        )
                entered += non_list_to_list( enter_dict_of_list( 
                        self.constant_solutions, 
                        boundry_squared, 
                        boundry_squared_solution, 
                        transform = transform 
                    ) )
            return entered
