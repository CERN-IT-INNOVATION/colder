

def check_physical_parameter_function(ff, dd):
    func, der = None, None
    if callable(ff):
        func = ff
        if callable(dd):
            der = dd
        elif dd is None:
            raise Exception('must provide a callable function for derivative')
        else:
            der = lambda any: dd
    else:
        func = lambda any: ff
        der = lambda any: 0
        
    return func, der