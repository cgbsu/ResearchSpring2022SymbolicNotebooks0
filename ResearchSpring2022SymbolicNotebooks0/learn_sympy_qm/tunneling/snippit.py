def iter_tree( expression, attribute ): 
    #print( "Iter" )
    #sp.pprint( expression )
    #object_directory = dir( expression )
    if attribute in dir( expression ): 
        if "keys" in dir( expression ):
            data = getattr( expression, attribute )
            #print( data.keys() )
            if len( data.keys() ) > 0:
                print( data )
        for i in dir( expression ): 
            iter_tree( getattr( expression, i ), attribute )
            
            

with sp.assuming( k_0_squared, sp.Q.real( total_energy ), sp.Q.positive( total_energy ), 
              sp.Q.real( mass ), sp.Q.positive( mass ), sp.Q.nonzero( mass ), 
              sp.Q.real( total_energy ), sp.Q.positive( total_energy ), sp.Q.nonzero( total_energy ), 
              total_energy > Potential.DEFAULT_POTENTIAL, 
              sp.Q.real( Potential.DEFAULT_POTENTIAL ), sp.Q.positive( Potential.DEFAULT_POTENTIAL ), 
              sp.Q.nonzero( Potential.DEFAULT_POTENTIAL ), 
              sp.Q.positive( hbar ), sp.Q.real( hbar ), sp.Q.nonzero( hbar ) ): 
    print( sp.ask( sp.Q.positive( k_0 ) ) )