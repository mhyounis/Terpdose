
#  ==================================
#    General utilities for Terpdose  
#  ==================================

class VersionError (RuntimeError):
    pass

def list_attributes (obj):
    
    # Lists a given object's attributes in a pretty way, ignoring '_'.
    # Returns an array of lines for use with the __repr__ function of 
    # a given class.
    
    import inspect
    
    lines = []
    lines.append('   Data')
    for key, val in vars(obj).items():
        if inspect.isfunction(val):
            continue
        else:
            if key.startswith('_'):
                continue
            if hasattr(val, 'shape'):
                shape = getattr(val, 'shape', None)
                dtype = getattr(val, 'dtype', None)
                lines.append(f'      {key} : shape={shape}, dtype={dtype}')
            elif isinstance(val, (int, float, complex, str, bool, bytes)):
                lines.append(f'      {key} : {val}')
            else:
                lines.append(f'      {key} : {type(val).__name__}')
    
    # If I ever reintroduce the cache thsi will be great to have.
    # # Include the cache as part of Data
    # if hasattr(obj, '_cache'):
    #     for key, val in obj._cache:
    #         if hasattr(val, 'shape'):
    #             shape = getattr(val, 'shape', None)
    #             dtype = getattr(val, 'dtype', None)
    #             lines.append(f'      {key} (cached) : shape={shape}, dtype={dtype}')
    #         elif isinstance(val, (int, float, complex, str, bool, bytes)):
    #             lines.append(f'      {key} (cached) : {val}')
    #         else:
    #             lines.append(f'      {key} (cached) : {type(val).__name__}')
    
    lines.append('   Methods')
    for key, val in inspect.getmembers(obj):
        if key.startswith('_'):
            continue
        if inspect.isroutine(val):
            lines.append(f'      {key}')
    
    return lines

def isinstancelist (obj, d) -> bool:
    
    ''' This function checks whether or not a list is of a particular data type. '''
    
    return isinstance (obj, list) and all(isinstance(x, d) for x in obj)

def lbsorted (keys):
    
    ''' This function takes a set of HDF5 keys of the form somestr_1, somestr_10, ..., and sorts    
        them by the trailing integer, so you'd have somestr_1, somestr_2, etc.                      
        somestr may itself contain _, so this function specifically defines the index as coming     
        after the final _.                                                                       '''
    
    def suffix_int (k):
        prefix, _, suffix = k.rpartition('_')
        return int(suffix)
    
    return sorted(keys, key=suffix_int)

def generic_repr (obj):
    
    # Lists all available data. Access by printing the object.
    
    lines = [f'<{type(obj).__name__}> Object:'] + list_attributes (obj)
    
    return '\n'.join(lines)