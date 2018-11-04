#-------------------------------------------------------------------------------
# Little module that writes LateX code about the functional that was input 
# into Tantalus.
#
#-------------------------------------------------------------------------------
from heph_densities import OrderOfDen, ParseOperators
from heph_densities import Densities_needed
from heph_fields import Fields_needed 
from heph_fields import Pairing_Fields_needed 
from heph_functional import Functional_terms 
from string        import Template

greek = [r'\alpha', r'\beta', r'\gamma', r'\delta', r'\eta']

operatordic = {}
operatordic['I'] = r'1'
operatordic['N'] = r'\nabla'
operatordic['S'] = r'\sigma'
operatordic['C'] = ''
operatordic['T'] = ''

template_location = 'tex_templates/'

def FillDensities():

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

        with open(template_location +'densities.tex', 'w') as f:
            f.write(exp)

def LateXDensity(density, ext_coup=[]): 
        #-----------------------------------------------------------------------
        # Returns a LateX expression for the density, represented in a  string.
        (der,lap,left,right,coup,cross) = ParseOperators(density)

        latex = ''
        if('C' in density):
            latex = 'C'
        else:
            latex = 'D'       
        #-----------------------------------------------------------------------
        superscript = ''
        subscript   = ''
        i = 0 
        # We don't need C's or I's anymore
        for i, op in enumerate((left + right).replace('C', '').replace('I','')):
            st = operatordic[op]
            for c in coup:
                if(i == c[0]):
                    st =  '(' + st
                elif(i == c[1]):
                    st = st + ')'
            if(i == len(left)-1):
                st = st + ','

            
            superscript = superscript + st
        #-----------------------------------------------------------------------
        #If just the density needs to be printed, then include free subscripts.
        if(len(ext_coup) == 0):
            order = OrderOfDen(density)            
            for i in range(order):  
                subscript = subscript + greek[i] 
            
    

        latex = latex + '^{%s}'%superscript + '_{%s}'%subscript 

        print latex
        return(latex)
