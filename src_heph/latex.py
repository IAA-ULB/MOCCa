#-------------------------------------------------------------------------------
# Little module that writes LateX code about the functional that was input 
# into Tantalus.
#
#-------------------------------------------------------------------------------
from heph_densities  import OrderOfDen, ParseOperators
from heph_densities  import Densities_needed
from heph_fields     import Fields_needed, ParseOperatorsField
from heph_fields     import Pairing_Fields_needed 
from heph_functional import Functional_terms , ParseDensities
from heph_functional import coupling_constants_0, coupling_constants_1, rreplace
from string          import Template
import os

greek = [r'\alpha', r'\beta', r'\gamma', r'\delta', r'\eta']

tex_location      = 'tex/'
template_location = 'tex_templates/'

def Build(funcfile, description):
    #---------------------------------------------------------------------------
    # Write .tex files containing info on the functional read from file.
    #---------------------------------------------------------------------------
    
    os.system('cp ' + template_location + 'main.tex ' + tex_location)        
    
    # comment characters do not play nice 
    # and it should be interpreted as a raw string
    raw = r'%s'%description.replace('#', '')
    with open(tex_location + 'info.tex', 'w') as f:
        f.write('Obtained from %s \n'%funcfile)
        f.write(r'\begin{verbatim}'  + '\n')        
        f.write(raw )     
        f.write(r'\end{verbatim}')  

    FillDensities()
    FillTerms()
    FillFields()
    FillSphamil()

    exit()

def FillDensities():
        #-----------------------------------------------------------------------
        # Write the densities.tex file that can be \input by main.tex.
        #-----------------------------------------------------------------------
        density_template = Template(r'\item $$ $DEN $$')

        exp = ''
        dic = {}
        latex = []
        for den in Densities_needed:
            #  Don't care for cross-couplings
            latex.append(LateXDensity(den.replace('x', '')))

        # Remove duplicates
        latex = list(set(latex))

        for l in latex:   
            dic['DEN'] = l
            exp = exp + density_template.substitute(dic) + '\n'

        with open(tex_location +'densities.tex', 'w') as f:
            f.write(exp)

def LateXDensity(density, ext_coup=[],nosub = 0): 
        #-----------------------------------------------------------------------
        # Returns a LateX expression for the density, represented in a  string.

        operatordic = {}
        operatordic['I'] = r'1'
        operatordic['N'] = r'\nabla'
        operatordic['S'] = r'\sigma'
        operatordic['C'] = ''
        operatordic['T'] = ''

        (der,lap,left,right,coup,cross) = ParseOperators(density)

        #-----------------------------------------------------------------------
        superscript = ''
        subscript   = ''
        
        prefix = ''
        for i in range(der):
            prefix = prefix + r'\nabla'
        for i in range(lap):
            prefix = prefix + r'\Delta'

        latex = prefix
        #  The 'G' is for when this routine is used to generate fields
        if('C' in density or 'G' in density):
            latex = latex + ' C'
        else:
            latex = latex + ' D'       
        
        i = 0 
        order = OrderOfDen(density)
        # We don't need C's or I's anymore
        leftreplaced  = left.replace('C', '') #.replace('I','')
        rightreplaced = right.replace('C', '')#.replace('I','')
        for i, op in enumerate(leftreplaced+ rightreplaced):
            st = operatordic[op]
            for c in coup:
                if(i == c[0]):
                    st =  '(' + st
                elif(i == c[1]):
                    st = st + ')'
            if(i == len(leftreplaced)-1):
                st = st + ','

            
            superscript = superscript + st
        #-----------------------------------------------------------------------
        #If just the density needs to be printed, then include free subscripts.
        if(nosub !=1 ):
            if(len(ext_coup) == 0):
                order = OrderOfDen(density)            
                for i in range(order):  
                    subscript = subscript + greek[i] 

        latex = latex + '^{%s}'%superscript + '_{%s}'%subscript 

        return(latex)

def FillTerms():
    #---------------------------------------------------------------------------
    # Write
    # a)  The functional terms to the functional.tex file.
    # b)  Coupling constant expressions to the constants.tex file
    
    term_template = Template(r'\item $$ $CPL  $TERM $$'  + '\n')
    cte_0_template  = Template(r' ${CPL}_0 &= ${EXP_0}' + r' \nonumber \\' + '\n')
    cte_1_template  = Template(r' ${CPL}_1 &= ${EXP_1}' + r'\nonumber \\' + '\n')

    exp = ''
    cte_0 = ''
    cte_1 = ''

    dic = {}
    for i, term in enumerate(Functional_terms):
        (den,coupl) = ParseDensities(term)
        
        cpl = 'C^{' 
        aux = ''        
        for d in den:
            aux = aux + LateXDensity(d)  + ' '
            cpl = cpl + LateXDensity(d,nosub=1)

        cpl = cpl + '}'

        dic['CPL']  = cpl
        dic['TERM'] = aux
        dic['EXP_0']= coupling_constants_0[i].replace('*', '')
        dic['EXP_1']= coupling_constants_1[i].replace('*', '')
 

        exp = exp + term_template.substitute(dic)
        cte_0 = cte_0 + cte_0_template.substitute(dic) 
        cte_1 = cte_1 + cte_1_template.substitute(dic)

    cte_0 = rreplace(cte_0, r'\\', '', 1)
    cte_1 = rreplace(cte_1, r'\\', '', 1)

    # a) Functional terms
    with open(tex_location +'functional.tex', 'w') as f:
            f.write(exp)
    # b) Coupling constants
    with open(tex_location +'constants_0.tex', 'w') as f:
            f.write(cte_0)
    with open(tex_location +'constants_1.tex', 'w') as f:
            f.write(cte_1)

def FillFields():
    #---------------------------------------------------------------------------
    #
    #---------------------------------------------------------------------------
        
    fields_template = Template(r'$FIELD &= \nonumber \\' + '\n')    

    exp = ''
    dic = {}

    for field in Fields_needed:
        (left, right, coupling, cross) = ParseOperatorsField(field)
        
        if(len(cross) > 0):
            continue        
        
        dic['FIELD'] = LateXField(field)
    
        exp = exp + fields_template.substitute(dic)    
    
    exp = rreplace(exp, r'\\', '', 1)
    with open(tex_location +'potentials.tex', 'w') as f:
        f.write(exp)

def LateXField(field):
    #---------------------------------------------------------------------------
    # Generate a LateX expression for a field.
    #
    #---------------------------------------------------------------------------

    den = LateXDensity(field)
    latex = den.replace('C', 'G').replace('D', 'F')

    return latex

def FillSphamil():
    #
    #
    #
    with open(tex_location +'sphamil.tex', 'w') as f:
        f.write('x^2 = x')
