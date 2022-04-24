import sympy as sp

def make_constant( expression, constant_name_base, constants ): 
        constant = sp.Symbol( str( constant_name_base ) + "_{" + str( len( constants ) ) + '}' )
        if expression in constants: 
            constant = constants[ expression ]
        constants[ expression ] = constant
        return constant, True

def bubble_constants( expression, free_variables : list, constant_name_base : str, constants : dict, tabs = 0, debug = False ): 
    is_atom = type( expression ) is sp.Symbol
    is_communative = expression.is_commutative
    if not getattr( expression, 'has' ): 
        return make_constant( expression, constant_name_base, constants )
    if expression.has( *tuple( free_variables ) ) and not is_atom: 
        substitution_table = {}
        new_expression = expression
        new_arguments = []
        replacable = []
        for argument in expression.args: 
            replacement_argument, replacable_atom = bubble_constants( argument, free_variables, constant_name_base, constants, tabs, debug )
            if replacable_atom and is_communative and not replacement_argument.has( *tuple( free_variables ) ): 
                replacable.append( replacement_argument )
            else: 
                new_arguments.append( replacement_argument )
        if len( replacable ) > 0: 
            new_constant = make_constant( 
                    type( expression )( *tuple( replacable ) ), 
                    constant_name_base, 
                    constants 
                )[ 0 ]
            new_arguments.append( new_constant )
        return type( expression )( *tuple( new_arguments ) ), True
    elif expression.is_Number: 
        return expression, False
    elif is_atom: 
        return expression, True
    else: 
        return make_constant( expression, constant_name_base, constants )

def group_constants( expression, free_variables : list, constant_name_base : str = 'S' ): 
    constants = {}
    return bubble_constants( 
                    expression, 
                    free_variables, 
                    constant_name_base, 
                    constants 
                )[ 0 ], constants
