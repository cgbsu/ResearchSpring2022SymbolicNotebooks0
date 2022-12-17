from custom_libraries.stepper import *

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

