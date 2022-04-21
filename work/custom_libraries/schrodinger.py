import sympy as sp
import sympy.physics.units.quantities as sq
from sympy.physics.quantum.constants import hbar
from custom_libraries.stepper import *
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

potential_none_to_list = lambda none_canidate : [] if none_canidate == None else none_canidate

class Boundries: 
    """TODO: (Speculative)
        * : Add relations between objects in boundry, such as =, >, <, <=, >=, and, or, xor
        * : Make no "special" group, but make it be able to have groups of groups of groups for ex.
        * : More granular remove_boundries, possibly that can work with boundries_with
        * : Somehow make boundries associative/communative <- sortof accomplished using append_boundries 
                and update_all_boundry_conditions*
    """
    
    BOUNDRY_ALL = "LastUpdatedAllBoundryConditions"
    BOUNDRY_QUERY_RESPONCE = "LastBoundryWithQueryResponce"
    
    INITIAL_COMMIT = "InitialCommit"
    COMMIT = "Commit"
    DEFAULT_ADD_BOUNDRIES_DEFAULT_COMMIT_PREFIX = "BeforeAdd"
    DEFAULT_APPEND_BOUNDRIES_DEFAULT_COMMIT_PREFIX = "BeforeAPPEND"
    DEFAULT_POST_ADD_BOUNDRIES_DEFAULT_COMMIT_PREFIX = "PostAdd"
    DEFAULT_POST_APPEND_BOUNDRIES_DEFAULT_COMMIT_PREFIX = "PostAPPEND"

    DEFAULT_REMOVE_SET_BOUNDRIES_DEFAULT_COMMIT_PREFIX = "BeforeRemoveSET"
    DEFAULT_REMOVE_DEFAULT_COMMIT_PREFIX = "BeforeRemove"
    DEFAULT_POST_REMOVE_SET_BOUNDRIES_DEFAULT_COMMIT_PREFIX = "PostRemoveSET"
    DEFAULT_POST_REMOVE_DEFAULT_COMMIT_PREFIX = "PostRemove"
    
    DEFAULT_NO_COMMENT_COMMIT = "NoComment"

    
    def __init__( self, initial_boundries = None ): 
        self.boundries = not_none_value( initial_boundries, {} )
        self.log = [ ( copy.deepcopy( self.boundries ), Boundries.INITIAL_COMMIT ) ]
    
    def commit( self, comment = None, prefix = None ): 
        self.log.append( ( 
                copy.deepcopy( self.boundries ), 
                str( prefix if prefix else "" ) \
                        + str( comment if comment else boundries.DEFAULT_NO_COMMENT_COMMIT ) \
                        + Boundries.COMMIT 
            ) )
        return self.log[ -1 ]
    
    def commit_log( self ): 
        return [ commit[ 1 ] for commit in self.log ]
    
    def add_boundries( self, set_name, boundries, automatically_append = False, append_override = False, 
                before_commit_prefix = DEFAULT_ADD_BOUNDRIES_DEFAULT_COMMIT_PREFIX, 
                post_commit_prefix = DEFAULT_POST_ADD_BOUNDRIES_DEFAULT_COMMIT_PREFIX ): 
        set_name_in_boundries = set_name in self.boundries
        """Prevents overwriting boundries unless explicitly specified to do so.
        It is computed twice, but if the assertion is thrown, the client may see 
        the explicit problem."""
        assert True if automatically_append == True else not set_name_in_boundries, """
                Attempt to set boundries when set already exists, set `automatically_append` 
                to `True` or call `append_boundries` to append
                """
        if set_name_in_boundries == True and automatically_append == True: 
            self.append_boundries( set_name, boundries, append_override )
        elif set_name_in_boundries == False: 
            self.commit( str( set_name ), before_commit_prefix )
            self.boundries[ set_name ] = boundries
            setattr( self, set_name, self.boundries[ set_name ] )
            self.commit( str( set_name ), post_commit_prefix )
        else: 
            assert set_name_in_boundries, "Attempt to overwrite existing boundries: " + set_name
        return set_name, self.boundries[ set_name ]
    
    def append_boundries( self, set_name, to_append, override = False, 
                enable_logging = True, 
                before_commit_prefix = DEFAULT_APPEND_BOUNDRIES_DEFAULT_COMMIT_PREFIX, 
                post_commit_prefix = DEFAULT_POST_APPEND_BOUNDRIES_DEFAULT_COMMIT_PREFIX ): 
            if enable_logging == True: 
                self.commit( str( set_name ), before_commit_prefix )
            boundries = self.boundries[ set_name ]
            to_append_keys = tuple( to_append )
            for key in to_append_keys: 
                value = to_append[ key ]
                if key in boundries: 
                    if value == boundries[ key ]: 
                        continue
                    elif not ( value in boundries ): 
                        boundries[ value ] = key
                    elif override: 
                        boundries[ key ] = value
                    else: 
                        assert key in boundries and ( not ( value in boundries ) ) and override
                else: 
                        boundries[ key ] = value
            if enable_logging == True: 
                self.commit( str( set_name ), post_commit_prefix )
            return set_name, boundries

    def combind_boundry_conditions( self, first_set, second_set, new_name = None, allow_update = False, allow_reset = False, enable_logging = False ): 
        first_set_is_string = type( first_set ) is str
        second_set_is_string = type( second_set ) is str
        passed_boundry_set = not ( first_set_is_string or second_set_is_string )
        assert passed_boundry_set if new_name == None else True, """
        If passing a boundry set, not the name of the boundry set, please specify a name of the resulting boundry_set
        """
        name = new_name if passed_boundry_set else new_name if new_name != None else first_set + second_set
        name_in_boundries = name in self.boundries
        assert name_in_boundries or new_name or allow_update or allow_reset, """
        Name already in boundry set, either specify 
        `True` for `allow_update` or `allow_reset` or specify new name
        """
        if not name_in_boundries or ( allow_reset and name_in_boundries ): 
            self.boundries[ name ] = {}
        assert not ( name_in_boundries and not allow_update ), """
        About to update boundry, when no update is allowed, to allow updates, please specify `allow_update` as `True`
        """
        self.append_boundries( name, self.boundries[ first_set ] if first_set_is_string else first_set, 
                enable_logging = enable_logging ), "first set", self.boundries[ first_set ]
        return self.append_boundries( name, self.boundries[ second_set ] if second_set_is_string else second_set, 
                enable_logging = enable_logging )
    
    def update_all_boundry_conditions( self, boundry_all_name = BOUNDRY_ALL, return_old = False ): 
        old = {}
        if boundry_all_name in self.boundries: 
            old = self.boundries.pop( boundry_all_name )
        boundry_keys = tuple( self.boundries.keys() )
        self.boundries[ boundry_all_name ] = {}
        for key in boundry_keys: 
            if key != boundry_all_name: 
                self.combind_boundry_conditions( boundry_all_name, key, boundry_all_name, True, False, False )
        return boundry_all_name, self.boundries[ boundry_all_name ]
    
    def remove_boundry_set( self, set_name, 
                pre_commit_prefix = DEFAULT_REMOVE_SET_BOUNDRIES_DEFAULT_COMMIT_PREFIX, 
                post_commit_prefix = DEFAULT_POST_REMOVE_SET_BOUNDRIES_DEFAULT_COMMIT_PREFIX ): 
        self.commit( str( set_name ), pre_commit_prefix )
        return_data = set_name, self.boundries.pop( set_name )
        self.commit( str( set_name ), post_commit_prefix )
        return return_data

    
    def remove_boundries( 
                self, 
                set_name, 
                keys : list, 
                values : list, 
                key_or_value : list, 
                pre_commit_prefix = DEFAULT_REMOVE_DEFAULT_COMMIT_PREFIX, 
                post_commit_prefix = DEFAULT_POST_REMOVE_DEFAULT_COMMIT_PREFIX 
            ): 
        self.commit( str( set_name ), pre_commit_prefix )
        removed = []
        boundry_set = self.boundries[ set_name ]
        keys = tuple( potential_none_to_list( keys ) + potential_none_to_list( key_or_value ) )
        removed = [ boundry_set.pop( key ) for key in keys ]
        values = tuple( potential_none_to_list( values ) + potential_none_to_list( key_or_value ) )
        removed = removed + [ boundry_set[ key ] for key in boundry_set.keys() if boundry_set[ key ] in values ]
        self.commit( str( set_name ), post_commit_prefix )
        return set_name, removed
    
    def has_set( self, set_name ): 
        return set_name in self.boundries
    
    def display( self ): 
        for boundry_set in self.boundries: 
            display( boundry_set )
            for boundry in self.boundries[ boundry_set ]: 
                display( sp.Eq( boundry, self.boundries[ boundry_set ][ boundry ] ) )
    
    def boundries_with_in_set( self, set_name, query ): 
        boundry_set = self.boundries[ set_name ]
        result = {}
        for key in boundry_set: 
            value = boundry_set[ key ]
            if key.has( query ): 
                result[ key ] = value
            elif value.has( query ): 
                result[ value ] = key
        return result
    
    def boundries_with( self, query, query_name = BOUNDRY_QUERY_RESPONCE ): 
        self.update_all_boundry_conditions( query_name )
        return self.boundries_with_in_set( query_name, query )                

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
                    sp.Symbol( self.normaliation_constant_base_name + "_{" + str( ii ) + '}' )
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
    
    def _equation_new_check_point( self, equation, check_point, check_point_number = None): 
        check_point_number = not_none_value( check_point_number, self._new_check_point_number() )
        return equation.check_point( 
                self.check_point_name_base \
                        + check_point \
                        + check_point_number 
            )
    
    def _equations_new_check_point( self, check_point, check_point_number = None ): 
        check_point_number = not_none_value( check_point_number, self._new_check_point_number() )
        for ii in range( len( self.equations ) ): 
            self._equation_new_check_point( self.equations[ ii ], check_point, check_point_number )
            self._equation_new_check_point( self.normalizations[ ii ], check_point, check_point_number )
    
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
                automatically_append = True,  
                simplifcication_table_name = BOUNDRY_BOUNDRY_SIMPLIFICATION_TABLE, 
                constant_table_name = BOUNDRY_BOUNDRY_CONSTANT_TABLE, 
                before_prefix = COMMIT_CHECK_POINT_PREFIX_BEFORE, 
                after_prefix = COMMIT_CHECK_POINT_PREFIX_POST
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
        constant_substitution_table = {}
        for key in self.boundries.boundries[ simplifcication_table_name ]: 
            constant_substitution_table[ boundry_simplification_list[ key ] ] = \
                    sp.Symbol( boundry_constant_base_name + "_{" + str( boundry_number ) + '}' )
            boundry_number += 1
        self.boundries.commit( constant_table_name, before_prefix )
        self._equations_new_check_point( 
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
                    .operate( lambda step : step.subs( boundry_simplification_list ), chain = True )
            self.normalizations[ ii ] \
                    .operate( lambda step : step.subs( boundry_simplification_list ), chain = True )
        for key in constant_substitution_table: 
            constant = constant_substitution_table[ key ]
            self.boundry_constant_symbols.append( constant )
            for ii in range( len( self.equations ) ): 
                try: 
                    self.equations[ ii ].replace_with_constant( key, constant )
                except Exception: 
                    pass
                try: 
                    self.normalizations[ ii ].replace_with_constant( key, constant )
                except Exception: 
                    pass
        self._equations_new_check_point( 
                after_prefix + TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_BOUNDRY_TO_CONSTANT_SUBSTITUTION
            )
        return self.equations, self.normalizations
    
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
        solved_for = []
        self.solved_boundries = []
        self._equation_new_check_point( equation, TimeIndependentSchrodingerConstantPotentials1D.CHECK_POINT_INTEGRATE_NORMALIZATION )
        if self.normalizations.index( equation ) >= 0 and auto_integrate == True: 
            equation.operate( lambda step : step.expand().doit() )
        for boundry in self.boundry_constant_symbols: 
            if equation.last_step().has( boundry ): 
                solutions = sp.solve( equation.last_step(), boundry )
                if solutions == None: 
                    continue
                if is_type_len_gt_0( list, solutions ): 
                    solved_for += self.add_constant_solutions( 
                            equation, 
                            boundry, 
                            solutions, 
                            assume_2_solutions_squared_absolute_value, 
                            transform 
                        )
                elif is_type_len_gt_0( dict, solutions ): 
                    solved_for += non_list_to_list( enter_lists_dict_of_list( 
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
                    solved_for += non_list_to_list( enter_dict_of_list( 
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
