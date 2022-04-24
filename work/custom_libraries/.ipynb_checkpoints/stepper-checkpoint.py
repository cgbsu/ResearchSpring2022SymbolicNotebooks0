import sympy as sp
import sympy.physics.units.quantities as sq
from sympy.physics.quantum.constants import hbar
import re
import copy

set_equal = lambda to_set, value : sp.Eq( to_set, value )
both_sides = lambda equation, operation : sp.Eq( operation( equation.lhs ), operation( equation.rhs ) )
equation_to_dict = lambda equation : { equation.lhs : equation.rhs }

not_none_value = lambda value, default : value if value != None else default

def display_steps( steppers, step = None ): 
    for stepper in steppers: 
        display( stepper.last_step( step ) )

def display_in( iterable ): 
    for ii in iterable: 
        display( ii )

transform_noop = lambda key, value : value

def set_key_to_equation( key, value ): 
    check_point = "set_key_to_equation : " + str( key ) + " : " + str( value )
    if check_point in value.check_points: 
        check_point += "#"
    value.check_point( "Before:" + check_point )
    value.add_step( sp.Eq( value.last_step(), key ) )
    value.check_point( "After:" + check_point )
    return value

def key_value_to_stepper( key, value, fix_stepper_equations = None, fix_sympy_equations = None, fix_type = sp.Eq ): 
    if type( value ) is Stepper: 
        if type( value.last_step() ) is fix_type: 
            return value
        elif fix_stepper_equations: 
            return fix_equations( key, value )
        return value
    elif type( value ) is fix_type: 
        if key == value.rhs or key == value.lhs: 
            return Stepper( value )
        elif fix_sympy_equations: 
            return Stepper( fix_sympy_equations( key, value ) )
    else: 
        return Stepper( value )
            

def append_list_in_dict_of_list( table : dict, key, value, transform = transform_noop ): 
        if type( value ) is list: 
            return table[ key ] + [ transform( key, current_value ) for current_value in value ]
        else: 
            return table[ key ].append( transform( key, value ) )

def enter_dict_of_list( table : dict, key, value, transform = transform_noop ): 
        if not key in table: 
            table[ key ] = [ transform( key, current_value ) for current_value in value ]
            return table[ key ]
        else: 
            return append_list_in_dict_of_list( table, key, value, transform )

def enter_lists_dict_of_list( table : dict, new_table : dict, enter = enter_dict_of_list, transform = transform_noop ): 
    return { key : enter( table, key, new_table[ key ], transform ) for key in new_table }

def braces_paired( 
            text, 
            open_brace_character = '(', 
            close_brace_character = ')' 
        ): 
    unmatched_braces = 0
    convert = { True : 1, False : 0 }
    for character in text: 
        unmatched_braces += convert[ character == open_brace_character ]
        unmatched_braces -= convert[ character == close_brace_character ]
    return unmatched_braces == 0

class FunctionSymbol: 
    def __init__( self, func, function_name, parameter_value_table = None ): 
        self.func = func
        self.function_name = function_name
        self.parameter_value_table = parameter_value_table \
                if type( parameter_value_table ) is dict \
                else { str( parameter_value_table ) : sp.Function( function_name )( *parameter_value_table ) } \
                        if parameter_value_table \
                        else {}
    
    def add_entry( self, parameter, value = None ): 
        key = parameter if type( parameter ) is tuple else ( parameter, )
        self.parameter_value_table[ key ] = \
                sp.Function( function_name )( *parameter_value_table ) \
                if value == None \
                else value
    
    def __call__( self, parameter ): 
        key = parameter if type( parameter ) is tuple else ( parameter, )
        return self.parameter_value_table[ str( key ) ]

class Symbols: 
    
    FORBIDDEN_IN_SYMBOL = r"[^a-zA-Z0-9]"
    NUMBER_REGEX = r"[0-9]"
    PERMISSABLE_PREFIXES = r"[a-zA-Z_]"
    IS_FUNCTION_REGEX = r"(.+\(+.+\)).+"
    FUNCTION_WITH_PARAMETER_REGEX = r".+\(.+\)"
    NUMBER_PREFIX = "_number_"
    MULTIPLY_REPLACE = "_multiply_"
    DIVIDE_REPLACE = "_divide_"
    ADD_REPLACE = "_add_"
    SUBTRACT_REPLACE = "_subtract_"
    RAISE_REPLACE = "_raise_"
    MATH_OPERATIONS_REPLACE = { 
            '^' : RAISE_REPLACE, 
            '-' : SUBTRACT_REPLACE, 
            '+' : ADD_REPLACE, 
            '*' : MULTIPLY_REPLACE, 
            '/' : DIVIDE_REPLACE
        }
        
    def _sanitize( 
                symbol, 
                particular_replace : dict = MATH_OPERATIONS_REPLACE, 
                prefix_fix_search = NUMBER_REGEX, 
                prefix_fix_replace = NUMBER_PREFIX, 
                universal_replace_search = FORBIDDEN_IN_SYMBOL, 
                universal_replace = '_', 
                admissable_prefix_sanity = PERMISSABLE_PREFIXES 
            ): 
        for to_replace, replace_with in particular_replace.items(): 
            symbol = symbol.replace( to_replace, replace_with )
        if re.match( prefix_fix_search, symbol[ 0 ] ):
            symbol = prefix_fix_replace + symbol
        symbol = re.sub( universal_replace_search, universal_replace, symbol )
        assert re.match( admissable_prefix_sanity, symbol[ 0 ] ) and len( symbol ) > 0
        return symbol
    
    def _push_function( 
                symbol, 
                symbol_value, 
                is_function = FUNCTION_WITH_PARAMETER_REGEX
            ):
        original_symbol = symbol
        symbol = str( symbol )
        function_match = re.match( is_function, symbol )
        function_match = ( symbol[ function_match.start() : function_match.end() ] if function_match else False )
        if function_match == symbol: 
            symbol = str( original_symbol.func )
            if symbol_value == original_symbol: 
                symbol_value = FunctionSymbol( original_symbol.func, symbol, original_symbol.args )
            elif type( symbol_value ) is dict: 
                symbol_value = FunctionSymbol( original_symbol.func, symbol, symbol_value )
            else: 
                symbol_value = FunctionSymbol( original_symbol.func, symbol, { original_symbol.args : symbol_value } )
            return symbol, symbol_value
        return None, None
    
    def __init__( self, *symbols, table = None ): 
        self.symbols = []
        if symbols or len( symbols ) > 0: 
            self.add_symbols( symbols )
        if table: 
            self.table_to_symbols( table )
    
    def add_symbol( 
                self, 
                symbol, 
                value = None, 
                sanitize = _sanitize, 
                push_function = _push_function
            ): 
        symbol_value = value if value != None else symbol
        original_symbol = symbol
        symbol = str( symbol )
        # More then one "function call" can match, so we see if only one does
        function_symbol, function_symbol_value = push_function( original_symbol, symbol_value )
        if function_symbol: 
            function_symbol = sanitize( function_symbol )
        symbol = sanitize( symbol )
        if function_symbol in dir( self ): 
            function_canidate = getattr( self, function_symbol )
            if type( function_canidate ) is FunctionSymbol: 
                function_canidate.add_entry( 
                        original_symbol.args, 
                        function_symbol_value( original_symbol.args ) 
                    )
        else: 
            setattr( self, str( function_symbol ), function_symbol_value )
        setattr( self, str( symbol ), symbol_value )
        self.symbols.append( symbol_value )
        return self
    
    def add_symbols( self, symbols ):
        for symbol in symbols: 
            self.add_symbol( symbol )
        return self
    
    def table_to_symbols( self, symbol_table : dict ): 
        for symbol in symbol_table.keys(): 
            self.add_symbol( symbol, symbol_table[ symbol ] )
        return self
    
    def symbol_by_string_name( 
                self, 
                symbol, 
                sanitize = _sanitize 
            ): 
        return getattr( self, sanitize( str( symbol ) ) )

equations_to_constant_table = lambda equations : { equation.lhs : Stepper( equation.reversed ) for equation in equations }

class Stepper: 

    LEFT = "LEFT", 
    RIGHT = "RIGHT"
    DEFAULT_APPLY = both_sides
    CONSTANT_SUBSTITUTION = lambda step, constant : step.subs( equation_to_dict( constant ) )    
    DEFAULT_GET_ELEMENT = lambda step, element : { 
            Stepper.LEFT : lambda step : step.lhs, 
            Stepper.RIGHT : lambda step : step.rhs, 
            0 : lambda step : step.lhs, 
            1 : lambda step : step.rhs, 
        }[ element ]( step )
    as_expressions = lambda table : { 
            symbol : table[ symbol ].last_step().rhs 
            for symbol in table.keys() 
        }
    as_equations = lambda table : { 
            symbol : table[ symbol ].last_step() 
            for symbol in table.keys() 
        }
    as_steppers = lambda table : { 
            symbol : table[ symbol ] 
            for symbol in table.keys() 
        }
    
    # Will be subscripted
    DEFAULT_CONSTANT_NAME_BASE = 'K'
    # Will be subscripted
    DEFAULT_CHECKPOINT_NAME_BASE = "checkpoint"
    
    def __init__( self, 
                first_step, new_steps = None, 
                apply = DEFAULT_APPLY, 
                constant_substitution = CONSTANT_SUBSTITUTION, 
                get_element = DEFAULT_GET_ELEMENT, 
                default_constant_name_base = DEFAULT_CONSTANT_NAME_BASE, 
                default_checkpoint_name_base = DEFAULT_CHECKPOINT_NAME_BASE, 
                default_assumptions = [], 
                constants = None, 
                closed = True, 
                solution_sets = None 
            ): 
        self.steps = new_steps if new_steps else []
        self.steps.append( first_step )
        self.chaining = False
        self.apply = both_sides
        self.constant_substitution = constant_substitution
        self.get_element = get_element
        self.check_points = {}
        self.check_point_steps = {}
        self.constants = not_none_value( constants, {} )
        self.default_constant_name_base = default_constant_name_base
        self.default_checkpoint_name_base = default_checkpoint_name_base
        self.assumptions = tuple( default_assumptions )
        self.closed = closed
        self.solution_sets = not_none_value( solution_sets, {} )
    
    def step_number( self, step = None ): 
        return ( len( self.steps ) - 1 ) if not step else step
    
    def last_step( self, step = None ):
        return self.steps[ self.step_number( step ) ]
    
    def _return_chain( self, step, chain ): 
        return self if ( chain or self.chaining ) else step        
    
    def add_step( self, new_step, chain = False ):
        self.steps.append( new_step )
        return self._return_chain( self.steps[ -1 ], chain )

    def operate( self, operation, step = None, chain = False ):
        self.steps.append( operation( self.last_step( step ) ) )
        return self._return_chain( self.steps[ -1 ], chain )

    def _apply_operation( self, apply ): 
        return self.apply if not apply else apply
    
    def _constant_substitution_operation( self, constant_substiution ): 
        return self.constant_substitution if not constant_substiution else constant_substiution

    def _get_element_operation( self, get_element ): 
        return self.get_element if not get_element else get_element
    
    def manipulate( self, operation, chain = False, apply = None ): 
        return self._return_chain( self.add_step( 
                self._apply_operation( apply )( self.last_step(), operation ) ), 
                chain 
            )
    
    def delete_step( self, step ): 
        old_step = self.steps[ step ]
        del self.steps[ step ]
        return old_step
    
    def undo( self, step = None, chain = False ): 
        return self._return_chain( 
                self.delete_step( self.step_number( step ) ), 
                chain 
            )
    
    def new_solution_set( self, solve_for, from_step = None, enable_overwrite = False ): 
        assert enable_overwrite if solve_for in self.solution_sets else True, f"""Error: Attempt to 
        overwrite existing solution set, `{str( solve_for )}`, if this is intentional, please 
        set the keyword argument `enable_overwrite` to `True` when calling `Stepper.add_to_solution_set`"""
        self.solution_sets[ solve_for ] = [ self.clone( from_step ) ]
        return self.solution_sets[ solve_for ]
    
    def append_solutions_to_sets( self, 
                solve_for = None, 
                solutions = None, 
                from_step = None, 
                automatically_make_new_solution_sets = False, 
                transform = None
            ): 
        enter = enter_dict_of_list if automatically_make_new_solution_sets else append_list_in_dict_of_list
        solutions = not_none_value( solutions, self.clone( from_step ) )
        assert solve_for or solutions, """Error: In order to run `Stepper.append_solutions_to_sets` you must specify 
                the `solve_for` and/or `solutions` parameter. If `solve_for` is not specified and `solutions` is not a `dict` 
                I dont know where to put solutions, and I cant append a default value to a set. If `solutions` and `solve_for` are 
                not specified I dont know what to do!"""
        if type( solutions ) == dict: 
            return enter_lists_dict_of_list( self.solution_sets, solutions, enter )
        assert solve_for, """Error: If you specify `solutions` (not as a dict { constant : solution }) 
                and not what `solve_for` in `Stepper.append_solutions_to_sets` I dont know where to put the solution!"""
        return enter( self.solution_sets, solve_for, solutions )
    
    # Not mutable
    def to_solution( self, solve_for, from_step = None, other_equations = None ): 
        last_step = self.last_step( from_step )
        return sp.solve( not_none_value( other_equations, [] ) + [ last_step ], solve_for )
    
    def clone( self, from_step = None ):
        return self.branch( lambda blank : blank, from_step )
    
    def branch( self, operation = None, from_step = None, apply = None ): 
        from_step = self.step_number( from_step )
        return Stepper( 
                self._apply_operation( apply )( self.last_step( from_step ), operation ), 
                self.steps[ : self.step_number( from_step ) : ], 
                self.apply, 
                self.constant_substitution, 
                self.get_element, 
                self.default_constant_name_base, 
                self.default_checkpoint_name_base, 
                copy.deepcopy( self.assumptions ), 
                copy.deepcopy( self.constants ), 
                copy.deepcopy( self.closed ), 
                copy.deepcopy( self.solution_sets )
            )
    
    def substitute_constant( self, constant, chain = False, constant_substitute = None ):
            in_constants = constant in self.constants
            if in_constants:
                return self.operate( lambda step : 
                        self._constant_substitution_operation( constant_substitute )( 
                                step, 
                                self.constants[ constant ].last_step().reversed 
                            ), 
                        chain = chain 
                    )
            else: 
                assert in_constants if self.closed else True 
                self.operate( lambda step : 
                        self._constant_substitution_operation( constant_substitute )( step, constant.reversed ), 
                        chain = chain 
                    ) 
    
    def replace_with_constant( self, to_replace, constant_name, last_step = None, chain = False ): 
        last_step = self.last_step( last_step )
        replace_result = last_step.replace( to_replace, constant_name, map = True )
        constant_name = sp.Symbol( str( constant_name ) )
        assert not constant_name in self.constants.keys(), """Error, attempt add a constant with a name that already 
        is entered""" + str( constant_name )
        assert len( replace_result[ 1 ] ) > 0, """Error, attempt to perform replacement failed because the desired 
        expression to replace with a constant does not exist!"""
        self.constants[ constant_name ] = Stepper( sp.Eq( to_replace, constant_name ) )
        self.add_step( replace_result[ 0 ] )
        return self._return_chain( self.constants[ constant_name ], chain )
    
    def constants_as_symbols( self, last_step = None ): 
        last_step = self.last_step( last_step )
        return Symbols( *tuple( self.constants.keys() ) )
    
    def constant_symbols( 
                self, 
                last_step = None,
                format_as = as_equations
            ): 
        return Symbols(
                table = format_as( self.constants )
            )
    
    def symbols( self, last_step = None ): 
        last_step = self.last_step( last_step )
        symbol_list = list( last_step.atoms() ) + list( last_step.atoms( sp.Function ) )
        return Symbols( *tuple( symbol_list ) )
    
    def rename( self, old_symbol, new_symbol_name, chain = False ): 
        last_step = self.last_step()
        assert old_symbol in last_step.atoms() or old_symbol in last_step.atoms( sp.Function )
        return self.add_step( last_step.subs( 
            old_symbol if type( old_symbol ) == dict \
                    else { old_symbol : new_symbol_name } ), chain 
                )
    
    def chain( self, operation ):
        self.chaining = True
        operation( self )
        self.chaining = False
        return self
    
    def begin_chain( self ): 
        self.chaining = True
        return self
    
    def end_chain( self ): 
        self.chaining = False
        return self
    
    def retrieve_element( self, element, from_step = None, get_element = None ): 
        return self._get_element_operation( get_element )( self.last_step( from_step ), element )
    
    # Creature comforts
    def left( self, from_step = None, get_element = None ): 
        return self.retrieve_element( Stepper.LEFT, from_step, get_element )

    def right( self, from_step = None, get_element = None ): 
        return self.retrieve_element( Stepper.RIGHT, from_step, get_element )
    
    def _default_name( self, name, name_base, number ): 
        return name_base + "_" + str( number ) if not name else name

    def _default_constant_name( self, name ): 
        return self._default_name( name, self.default_constant_name_base, len( self.constants ) )

    def _default_check_point_name( self, name ): 
        return self._default_name( name, self.default_checkpoint_name_base, len( self.check_points ) )
    
    def _default_name( self, name, name_base, number ): 
        return name_base + "_" + str( number ) if not name else name

    def add_assumption( self, assumption, chain = False ): 
        self.assumptions = tuple( list( self.assumptions ) + [ assumption ] )
        return self._return_chain( self.assumptions, chain )
    
    def element_to_constant( self, element, constant_name = None, from_step = None, chain = False, get_element = None, assumptions = None ): 
        assert not constant_name in self.constants.keys(), """Error, attempt add a constant with a name that already 
        is entered""" + str( constant_name )
        assumptions = not_none_value( assumptions, {} )
        constant = sp.Symbol( self._default_constant_name( constant_name ), **assumptions )
        self.constants[ constant ] = Stepper( sp.Eq( constant, self.retrieve_element( element, from_step, get_element ) ) )
        return self._return_chain( self.constants[ constant ], chain )

    def left_to_constant( self, constant_name = None, from_step = None, chain = False, get_element = None, assumptions = None ): 
        return self.element_to_constant( Stepper.LEFT, constant_name, from_step, chain, get_element, assumptions )

    def right_to_constant( self, constant_name = None, from_step = None, chain = False, get_element = None, assumptions = None ): 
        return self.element_to_constant( Stepper.RIGHT, constant_name, from_step, chain, get_element, assumptions )
    
    def check_point( self, name_or_marker = None, from_step = None, chain = False, pass_check_point_name = False ): 
        checkpoint_name = self._default_check_point_name( name_or_marker )
        self.check_points[ checkpoint_name ] = self.last_step( from_step )
        self.check_point_steps[ checkpoint_name ] = self.steps.index( self.last_step( from_step ) )
        return self._return_chain( checkpoint_name if pass_check_point_name else self.check_points[ checkpoint_name ], chain )
    
    def restore_from_check_point( self, name_or_marker, add_not_undo = True, chain = False ): 
        assert name_or_marker in self.check_points
        step = None
        if add_not_undo: 
            step = self.add_step( self.check_points[ name_or_marker ] )
        else: 
            check_point_index = self.check_point_steps[ name_or_marker ]
            if check_point_index >= 0: 
                for ii in range( self.step_number() - check_point_index ): 
                    self.undo()
                for key in self.check_point_steps.keys(): 
                    if self.check_point_steps[ key ] > check_point_index: 
                        self.check_point_steps[ key ] = -1
        return self._return_chain( self.check_points[ name_or_marker ], chain )
    
    def constant_names( self ): 
        return list( self.constants.keys() )
    
    def check_point_markers( self ): 
        return list( self.check_points.keys() )
    
    # If you are wondering where I got the funny names from: 
    # https://en.wikipedia.org/wiki/Multiplication
    
    def multipy( self, multiplicand, chain = False, apply = None ): 
        return self.manipulate( lambda side : side * multiplicand, chain, apply )
    
    def divide( self, denomonator, chain = False, apply = None ): 
        return self.manipulate( lambda side : side / denomonator, chain, apply )
    
    def add( self, addend, chain = False, apply = None ): 
        return self.manipulate( lambda side : side + addend, chain, apply )
    
    def subtract( self, subtrahend, chain = False, apply = None ): 
        return self.manipulate( lambda side : side - subtrahend, chain, apply )
    
    def raise_to_power( self, power, chain = False, apply = None ): 
        return self.manipulate( lambda side : side ** power, chain, apply )
    
    def root( self, degree, chain = False, apply = None ): 
        return self.manipulate( lambda side : sp.root( side, degree ), chain, apply )

    def __imul__( self, multiplicand ): 
        return self.manipulate( lambda side : side * multiplicand, True )
    
    def __itruediv__( self, denomonator ): 
        return self.manipulate( lambda side : side / denomonator, True )
    
    def __iadd__( self, addend ): 
        return self.manipulate( lambda side : side + addend, True )
    
    def __isub__( self, subtrahend ): 
        return self.manipulate( lambda side : side - subtrahend, True )
    
    def __ipow__( self, power ): 
        return self.manipulate( lambda side : side ** power, True )
