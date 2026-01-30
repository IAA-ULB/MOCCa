"""
  This module provides functionality to generate custom code that MOCCaPy can rely on.
"""
from importlib import import_module
from pathlib import Path
from string import Template

# TODO: how to properly pass this around?
spaces_per_indent = 4

# def generate_fortran(heph_input, output_dir):
#     """
#       Generate fortran source code for several purposes.

#       Args :

#         heph_input : dictionary
#           a dictionary containing all relevant input for a Hephaestos run

#           It currently requires the following keys
#          'func_file' : the file specifying the EDF being generated

#         output_dir : string
#           directory name where the various .f90 source code
#           files will be written to.
#     """

#     return

class FunctionalGenerator:
    """Create a EDF class from a `.func` file or a `.func.json` file.

    Remark:
        This class wraps the methods provided by Wouter.
    """
    def __init__(self, func_path:Path):
        """
        Args:
            func_path : Path to a `.func` file or a `.func.json` file
        """
        self.func_path = func_path
        self.edf_specification = read_functional_from_file(func_path)

    def generate_module(self, location:Path=Path('.')):
        """Generate a Python module with the EDF class from the specification read in the constructor.
        Args:
            location : Path to folder where Python module is written.
        """

        class_name = edf_class_name(self.edf_specification['func_name'])
        self.edf_specification['class_name' ] = class_name
        self.edf_specification['module_name'] = class_name.lower()
        self.edf_specification['func_path'  ] = self.func_path

        generate_EDF_class(edf_specification=self.edf_specification, location=location)

    def load(self):
        """Load the python module.

        Returns:
            a reference to the python module loaded.
        Raises:
            ModuleNotFoundError if the python module could not be loaded. It may not exist, or
            its location may not be on the search path.
        """
        mod = import_module(self.edf_specification['module_name'])
        return mod
#--------------------------------------------------------------------------------------------------
# functions provided by Wouter.

def edf_class_name(func_name:str):
    """
    Args:
        func_name : string name of the .func file (stem, without the extension)
    Returns:
        string: dashes and blanks replaced with underscores. First character is capitalized
    Raises:
        ValueError: if func_name does not start with a alphabetical character.
    """
    if func_name[0].isalpha():
        class_name = func_name[0].upper() + func_name[1:]
    else:
        raise ValueError(
            f"The name of a `.func` file must start with an alphabetical character.\n"
            f"Otherwise it cannot be converted to a class name. '{func_name}' does not comply."
        )
    return class_name.replace('-', '_').replace(' ', '_')


def generate_EDF_class(edf_specification, location:Path=Path('.')):
    """
    Generate Python source code for several purposes.

    Args :
        edf_specification : a dictionary containing all relevant input for a Hephaestos run
            It currently requires the following keys
            'func_file' : the file specifying the EDF being generated

        location: Path to the folder where Python module is written.

    Remark:
        the original output_file argument is removed because the module name and the class
        name to which the new EDF class file will be written is now
        automatically generated from the .func filename
    """

    # Relevant keys currently
    # $DESCRIPTION : a string to place in the comments
    # $CLASS_NAME  : name of the EDF class, derived froM the name of the `.func` file, agreeing
    # $FUNC_NAME   : the name of the `.func` file, agreeing

    # with standard Python naming conventions. i.e. 'BXL-N2LO-full' becomes  BXL_N2LO_full.

    if not 'module_name' in edf_specification:
        class_name = edf_class_name(edf_specification['func_name'])
        edf_specification['class_name' ] = class_name
        edf_specification['module_name'] = class_name.lower()

    module_path = location / (edf_specification['module_name'] + '.py')

    generated_code = {}

    generated_code['FUNC_NAME'  ] = edf_specification['func_name']
    generated_code['CLASS_NAME' ] = edf_specification['class_name']
    generated_code['DESCRIPTION'] = '\n'.join(edf_specification['description'])

    coupling_constants_calculation = ''

    # Hephaestos first writes some code to shorten the expressions that come after
    for param in edf_specification['parameters']:
        coupling_constants_calculation +=  2*spaces_per_indent*' ' \
                                        +  '%-10s  = self.param.%s\n'%(param,param)

    coupling_constants_calculation += '\n'
    for cc_fortran, term, iso, den_dep in zip(edf_specification['coupling_constants'],
                                              edf_specification['functional_terms']  ,
                                              edf_specification['isospin_indices']   ,
                                              edf_specification['density_dependence'],
                                              ):
        term_string = "'" + term + "'"
        cc_python = parse_cc_expression(cc_fortran, edf_specification)
        den_dep_str = "'" + den_dep + "'"

        expression = 2*spaces_per_indent*' ' \
                   + "self.coupling_constants[(%-40s,%-15s,%-10s)] = %s"%(term_string, tuple(iso), den_dep_str, cc_python)
        coupling_constants_calculation +=  expression + '\n'

    # Remove the final newline
    coupling_constants_calculation   = coupling_constants_calculation[:-1]
    generated_code['CC_CALCULATION'] = coupling_constants_calculation

    # As it is assumed that the link between a .func file and the Python module
    # with the corresponding EDF class is unique we do NOT check that the Python module
    # exists. If it does, it is simply overwritten.
    with open(module_path, mode='w') as f:
        text = substitute('mocca/hephaestos/templates/EDF_class.py', generated_code)
        f.write(text)


    # An implementation of a specific EDF implementation should have
    # - a way to read a parameterisation file
    # - a way to calculate coupling constants; results stored inside a param object
    # - some way to retain the information on what densities are relevant
    # - some way to retain the information on the corresponding potentials
    # - a way to calculate the action of sphamil on single-particle wavefunctions
    # - a way to calculate the energy of a given DensityVector

def parse_cc_expression(fortran_exp, edf_specification):
    """
        Args:
            fortran_exp :
                fortran-valid expression for the calculation of the
                coupling constant
            edf_specification: dictionary
                complete specification of the entire EDF, required to
                identify parameters
        Returns:
            python_exp:
                a valid Python expression for the calculation of the
                coupling constant
    """

    python_exp = fortran_exp
    # replace all calls to the functions CC() and CT() with
    #  CC ->  skyrme_cc
    #  CT ->  skyrme_ct
    python_exp = python_exp.replace('Cc(', 'skyrme_cc(')
    python_exp = python_exp.replace('Ct(', 'skyrme_ct(')

    # Remove all 'd0'
    python_exp = python_exp.replace('d0', '')

    return python_exp

# def build(source_dir):
#     """
#       Construct Python interfaces for all fortran source code encountered in
#       a given source directory.

#       Args:
#         source_dir : string
#                     directory containing all relevant .f90 source code files.
#     """

#     return

def read_functional_from_file(fname):
    """
     Read the details of the EDF terms, their structure and coupling constants,
     from the file named fname.

     Args:
        fname : string or Path
            the path to the file containing the functional specification

     Returns:
        edf_specification : dictionary
            this contains the following keys

            description : string
                a short description of the functional, useful for human
                identification purposes
            parameters  : list of strings
                a list of the (names of) EDF parameters
                e.g. [ 't0', 'x0', ... ]
            parameter_types : list of strings
                a list of the types of the parametrs in the parameters list
                e.g. [ 'real', 'integer', ... ]
            functional_terms : list of strings
                a list of functional terms
                e.g. [ 'E_D_I_I_D_I_I', ....]
            coupling_constants: list of strings
                a list of coupling constants
                e.g. [ 'Cc(t0,x0,+1,0,0)', ... ]
            isospin_indices : a list of lists of strings
                a list of the isospin indices of the terms
                each index can be '0', '1', 'n', or 'p'
                e.g. [ ['0', '0'], .... ]
            density_dependence: a list of strings
                a list of the exponent to apply to the first density in a given term
                e.g. [ '1', '1', ..., 'beta', ...]
            extra_calls : a list of strings
                a list of any extra function calls to be performed when evaluating
                a coupling constant; mostly useful to accomodate the use of
                microscopic pairing recipes
                e.g. ['', '', ..., 'vmicro()', ...]

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     An EDF specification file should be composed of

      Part 1: description of the functional (in words)
      Part 2: enumeration of the parameters
      Part 3: !TERMS => signalling the end of Part 2 and the start of Part 1
      Part 4: term specification

     All of these parts can be interspersed with lines starting with '!'; these
     are comment lines and do not influence the code generation in any way.

      Part 1: description of the type of functional, which will be included in
              the Hephaestos output.
              Lines need to start with '#'

      Part 2: enumeration of all parameters that should be read from a .param
              file. Entries should be separated by ';'.

      Part 3: '!TERMS' (no modification, EVER)

      Part 4: line-by-line specification of all the terms in the functional.
              These should have the form

              E_[D1]_[D2]_[D3]_[D4] ; C ; alpha ; iso_1 ; iso_2 ; iso_3 ; iso_4
                 (1)                 (2)  (3)      (4)

              (1)    enumeration of the densities in the term, including the way
                     they are coupled. Example:

                      E_D_I_Sm_D_I_Sm =  sum_{mu=x/y/z} s_mu(r) s_mu(r)

                     this version of the code allows for
                      (a) bilinear (two densities)
                      (b) trilinear
                      (c) quadrilinear terms

               (2)   coupling constant of the term. Can be given in terms of the
                     Cc function coded in the Fortran templates.

               (3)   density dependence of the FIRST density, D1.
                     If alpha != 1, the code will enforce iso_1 to be zero

               (4)   isospin indices of the densities; 0, 1, 'p' or 'n'.
                     normal densities should have isospin indices (0,1) and
                     pairing densities should have p/n indices ('p', 'n')

     Special exception: for pairing terms, the code can read ONE ADDITIONAL
     entry into the row, useable to make the code call an additional extra
     routine to calculate energies and fields.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """
    fpath = Path(fname).resolve()
    with open(fpath, 'r') as f:

        edf_specification = {
            'func_name'          : fpath.stem,
            'description'        : [],
            'parameters'         : [],
            'parameter_types'    : [],
            'functional_terms'   : [],
            'coupling_constants' : [],
            'isospin_indices'    : [],
            'density_dependence' : [],
            'extra_calls'        : [],
        }

        termsstart = False
        for line in f:
            try:
                if len(line.split()) == 0:
                    # Forget about empty lines
                    continue

                if line[0] == '#':
                    # '#'-signs indicate PART 1, the tdescription of the EDF
                    if line.startswith('# '):
                        line = line[2:]
                    line = line.strip()
                    if len(line) > 1:
                        edf_specification['description'].append(line)
                    continue

                if line[0:6] == '!TERMS' :
                    # Signal that the parameter specification, PART 2 is over.
                    termsstart = True
                    continue

                if line[0] == '!' :
                    # comment line, don't do anything with it
                    continue

                if termsstart == False :
                    # Parse the parameters of the functional in PART 2
                    param_descriptions = [w.split() for w in line.split(';')]
                    for type,name in param_descriptions:
                        if type not in ('RI'):
                            raise ValueError(
                                f"Hephaestos encountered invalid parameter type '{type}', expecting 'R'|'I'.\n"
                                f"  while processing line:\n"
                                f"  '{line}'\n"
                                f"  in file:\n"
                                f"  '{fpath}'."
                            )
                        edf_specification['parameters'].append(name)
                        edf_specification['parameter_types'].append(type)

                elif termsstart == True :
                    #  Start the actual terms of the functional in PART 4
                    split = [w.strip() for w in line.split(';')]
                    edf_specification['functional_terms'  ].append(split[0])
                    edf_specification['coupling_constants'].append(split[1])
                    edf_specification['density_dependence'].append(split[2])

                    # check how many densities are in this particular term...
                    densities,_ = parse_edf_term(split[0])
                    # ... such that we know how many isospin indices to expect!
                    iso_list = []
                    for k in range(3,3+len(densities)):
                        iso = split[k].strip()
                        if iso not in ('0', '1', 'p', 'n'):
                            raise IndexError
                        iso_list.append(iso)
                    edf_specification['isospin_indices'].append(iso_list)

                    if len(split)>3+len(densities):
                        extra = split[-1].strip()
                        edf_specification['extra_calls'].append(extra)
                    else:
                        edf_specification['extra_calls'].append('')

            except Exception:
                raise RuntimeError(
                    f"Hephaestos failed interpreting line:\n"
                    f"    {line}\n"
                    f"in `.func` file {fname}"
                )

    return edf_specification

def identify_param(param_string):
    """
    Identify a parameter from a string read from a functional file.
    More specifically, this function expects a string such as

            R x0

    which means that the code will depend on a REAL (R) parameter named x0.
    The alternative is

            I sigma

    which indicates an integer parameter sigma.

    Args:
        param_string: string

    Returns:
        param_type  : string
            type of the parameter (integer/real)
        param       : string
            name for the parameter
    """
    # First, strip sole "R" and "I"
    if 'R ' == param_string[0:2]:
        param_type = 'real'
    elif 'I ' == param_string[0:2]:
        param_type = 'integer'
    else:
        raise ValueError(
            f"Hephaestos faced an unrecognized parameter type '{param_type}', expecting 'R' or 'I'."
        )
        # print ('Offending entry: ', param_string)
        # sys.exit(1)

    param = param_string[2:].strip()
    return param_type, param

# replace with a.strip()
# def clean(a):
#     """
#         Quick'n'dirty string cleaning routine, strips spaces and newlines
#     """
#     return a.replace(' ', '').replace('\n', '')

def parse_edf_term(term):
    """
      We deconstruct a term in the functional.

      Args:
        term : string
            Specifies a term in the EDF

      Returns:
        densities: list of strings
            a list of the densities participating in a term
        coupling : list of tuples
            a list of the 'linked' indices in a given term

      Example:
            E_D_I_Sm_Derxm_C_I_Nxm
      leads to
          densities : D_I_S, Der_C_I_N
          coupling  : [(0,1,2)]

    """
    # TODO: what is the best way to propagate these 'hardcoded' things into
    #       the Hephaestos module?
    sumindices      = ['m', 'k', 'q', 'o', 'l']
    crossindices    = ['x', 'y', 'z']
    lapstring       = 'Lap'
    derstring       = 'Der'

    #  Split the input string along the underscores; extract the densities
    densities = []
    temp      = ''
    split     = term.replace('_DD', '').split('_')
    for i,_ in enumerate(split):
        if split[i][0:3] == derstring or split[i] == lapstring:
            temp  = temp + split[i] + '_'

        if split[i] == 'D' or split[i] == 'C' \
                           or split[i] == 'DP' or split[i] =='CP':
            temp = temp + split[i] + '_' + split[i+1] + '_' + split[i+2]
            densities.append(temp)
            temp = ''
    # Find all couplings by looping over all possible accepted summation letters
    coupling  = []
    foundx    = []
    for contracted_index in sumindices:
        c   = ()
        ind = 0
        foundx.append(0)
        for i in range(len(term)):
            if term[i] == contracted_index:
                c = c+ (ind,)
            if term[i] in sumindices:
                ind = ind + 1
        if len(c) > 0 :
            coupling.append(c)

    # Don't propagate couplings into the name that are not between
    # left and right operators
    for contracted_index in sumindices:
        for i in range(len(densities)):
            if derstring + contracted_index in densities[i] :
                # Remove the coupling if it involves derivatives
                densities[i] = densities[i].replace(contracted_index, '')
                for j in range(len(densities)):
                    densities[j] = densities[j].replace(contracted_index, '')
            for j in range(len(densities)):
                if i == j:
                    pass
                elif contracted_index in densities[i] and contracted_index in densities[j]:
                    # Remove the coupling if it is between more densities
                    densities[i] = densities[i].replace(contracted_index, '')
                    densities[j] = densities[j].replace(contracted_index, '')
                else:
                    pass

    for contracted_index in sumindices:
        for x in crossindices:
            for i in range(len(densities)):
                if derstring + x + contracted_index in densities[i] :
                    # Remove the coupling if it involves derivatives
                    densities[i] = densities[i].replace(contracted_index, '')
                    for j in range(len(densities)):
                        densities[j] = densities[j].replace(contracted_index, '')

    # Remove any vector coupling indices that might remain
    for x in crossindices:
        for i in range(len(densities)):
            densities[i] = densities[i].replace(x, '')

    return densities, coupling

def substitute(src_template, dictionary):
    """
    Substitute strings in a template file and write the result to a target file,
    taking care to place annotations that indicate where substitutions have been made.

    Input:
        src: Source file (template)
        dictionary: Dictionary with substitutions

    Output:
        None (writes to target file)
    """
    code = ''
    with open(src_template, 'r') as template:
        for line in template:
            newline = Template(line).substitute(dictionary)
            #if( newline != line):
            #    # Some string substitution happened
            #    generated.write(annotate(newline))
            #else:
            # No substitution happened, just copy the line
            code += newline

    return code

# def annotate(line):
#     TODO: this function is not suited to Python YET
#     """
#     Annotate a line that has undergone substitution with comments indicating
#     that a substitution has occurred. If the line contains multiple lines,
#     each line except the first and last will be on its own line with a comment
#     indicating the start and end of the substitution block.
#
#     Input:
#         line: The line to annotate
#     Output:
#         newline: the annotated line.
#     """
#     split = line.splitlines()
#     if(len(split)>1):
#         newline = split[0] + ' !Hephaestos<<< \n'
#         for k in range(1, len(split)-1):
#             newline = newline + split[k] + '\n'
#         newline = newline + split[-1] + ' !Hephaestos>>> \n'
#     else:
#         try:
#             if('#' == split[0][0]):
#                 newline = split[0] + ' /* Hephaestos substitution */ \n' # Preprocessor directives have a different comment syntax.
#             else:
#                 newline = split[0] + ' ! Hephaestos substitution \n'
#         except IndexError:
#             newline = line
#     return newline
