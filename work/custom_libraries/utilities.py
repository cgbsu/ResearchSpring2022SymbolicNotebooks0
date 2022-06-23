import sympy as sp

set_equal = lambda to_set, value : sp.Eq( to_set, value )
both_sides = lambda equation, operation : sp.Eq( operation( equation.lhs ), operation( equation.rhs ) )
both_sides_no_evaluation = lambda equation, operation : sp.Eq( operation( equation.lhs ), operation( equation.rhs ), evaluate = False )
equation_to_dict = lambda equation : { equation.lhs : equation.rhs }

not_none_value = lambda value, default : value if value != None else default

default_format_constant_group_constants = lambda constant_name_base, constants : str( constant_name_base ) + "_{" + str( len( constants ) ) + '}'

flip_dictonary = lambda dictonary : { dictonary[ key ] : key for key in dictonary }

def make_constant( 
                expression, 
                constant_name_base, 
                constants, 
                assumptions, 
                format_constant 
            ): 
        constant = sp.Symbol( format_constant( constant_name_base, constants ), **assumptions )
        if expression in constants: 
            constant = constants[ expression ]
        constants[ expression ] = constant
        return constant, True

def bubble_constants( 
            expression, 
            free_variables : list, 
            constant_name_base : str, 
            constants : dict, 
            assumptions, 
            format_constant, 
            tabs = 0, 
            debug = False 
        ): 
    is_atom = type( expression ) is sp.Symbol
    is_communative = expression.is_commutative
    if not getattr( expression, 'has' ): 
        return make_constant( expression, constant_name_base, constants, assumptions, format_constant )
    if expression.has( *tuple( free_variables ) ) and not is_atom: 
        substitution_table = {}
        new_expression = expression
        new_arguments = []
        replacable = []
        for argument in expression.args: 
            replacement_argument, replacable_atom = bubble_constants( 
                    argument, 
                    free_variables, 
                    constant_name_base, 
                    constants,  
                    assumptions, 
                    format_constant, 
                    tabs, 
                    debug 
                )
            if replacable_atom and is_communative and not replacement_argument.has( *tuple( free_variables ) ): 
                replacable.append( replacement_argument )
            else: 
                new_arguments.append( replacement_argument )
        if len( replacable ) > 0: 
            new_constant = make_constant( 
                    type( expression )( *tuple( replacable ) ), 
                    constant_name_base, 
                    constants, 
                    assumptions, 
                    format_constant 
                )[ 0 ]
            new_arguments.append( new_constant )
        return type( expression )( *tuple( new_arguments ) ), True
    elif expression.is_Number: 
        return expression, False
    elif is_atom: 
        return expression, True
    else: 
        return make_constant( expression, constant_name_base, constants, assumptions, format_constant )

def group_constants( 
            expression, 
            free_variables : list, 
            constant_name_base : str = 'S', 
            constants : dict = None, 
            assumptions = { 'real' : True, 'finite' : True }, 
            format_constant = default_format_constant_group_constants 
        ): 
    constants = not_none_value( constants, {} )
    return bubble_constants( 
                    expression, 
                    free_variables, 
                    constant_name_base, 
                    constants, 
                    assumptions, 
                    format_constant, 
                )[ 0 ], constants

def ungroup_constants( expression, constants ): 
    flipped = flip_dictonary( constants )
    while expression.has( *tuple( flipped.keys() ) ) == True:
        expression = expression.subs( flipped )
    return expression
