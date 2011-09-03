#$ID:   $
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['BayesianBlocks'])
    env.Tool('astroLib')
    env.Tool('st_facilitiesLib')
    env.Tool('LikelihoodLib')
    env.Tool('irfInterfaceLib')
    #probably also need irfLoaderLib

def exists(env):
    return 1
