from operatortospwf import *
from string         import Template

densities =[]
dimensions=[]

def initdensities():
    global densities, dimensions
    
    densities      = ['rho', 'tau']
    dimensions     = [  0  ,   2  ]

def processdensities(fname, src, target):
    #===========================================================================
    # Writing the code to compute the densities
    
    #===========================================================================
    # Some templates
    decformat  = '    real(KIND=dp), allocatable :: DENSITY(:,:,:,EXTRADIM)'
    
    iniformat  =             '    if(.not.allocated(DENSITY)) then \n'
    iniformat  = iniformat + '       allocate(DENSITY(nx,ny,nz,EXTRADIM)) \n'
    iniformat  = iniformat + '    endif \n'
    iniformat  = iniformat + '    DENSITY = 0.0d0 '
    
    compformat = '          DENSITY(i,1,1,EXTRADIM) = DENSITY(i,1,1,EXTRADIM)  '
    #===========================================================================

    
    declaration   =''
    compute       =''
    initialization= ''
    
    loopindices = ['mu', 'nu', 'ka']
    
    for den in densities:
    
        dim = dimensions[densities.index(den)]
    
        #-----------------------------------------------------------------------
        #Fixing correct allocation
        extradim=''
        for j in range(dim):
            extradim = extradim + ':,'
        extradim = extradim + ':'
        
        adddec = decformat.replace('DENSITY', den).replace('EXTRADIM', extradim)
        declaration    = declaration    + adddec + '\n'
        
        #-----------------------------------------------------------------------
        #Fixing correct initialisation
        extradim=''
        for j in range(dim):
            extradim = extradim + '3,'
        extradim = extradim + '2'
        
        addini  = iniformat.replace('DENSITY', den).replace('EXTRADIM', extradim)
        initialization = initialization + addini + '\n'
        
        #-----------------------------------------------------------------------
        # Fixing correct calculation

        if (den =='rho'):
            left  = 'HFBasis'
            right = 'HFBasis'
        elif ( den  == 'tau' ):
            left  = 'HFderiv'
            right = 'HFderiv'
        
        if(dim == 0) :
            #Scalar density
            extradim = 'it'
            compute        = compute  + compformat.replace('DENSITY', den).replace('EXTRADIM', extradim)
            compute        = compute  + GenDensityExpression( left, 'j', right, 'j', Identity, Identity)
            compute        = compute  + '\n'
        elif(dim == 1):
            # vector density
            extradim = 'mu,it'
            
            for mu in [1,2,3] :                
                compute        = compute  + compformat.replace('DENSITY', den).replace('EXTRADIM', extradim)
                compute        = compute  + GenDensityExpression( left, 'mu,j', right, 'mu,j', Identity, Identity)
                compute        = compute  + '\n'
                
        elif(dim == 2):
            # rank 2 tensor density
            for mu in [1,2,3] :
                for nu in [1,2,3] :
                    extradim = '%d,%d,it'%(mu,nu)
                    compute        = compute  + compformat.replace('DENSITY', den).replace('EXTRADIM', extradim)
                    compute        = compute  + GenDensityExpression( left, '%d,j'%mu, right, '%d,j'%nu, Identity, Identity)
                    compute        = compute  + '\n'
        else :
            print 'Defined density with too high dimension'
            exit()                      
        
    dic={}
    dic['DECLARATION']    = declaration
    dic['INITIALIZATION'] = initialization
    dic['EXPRESSION']     = compute    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  
                
  
