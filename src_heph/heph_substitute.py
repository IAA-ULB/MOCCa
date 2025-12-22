#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
#
#===============================================================================

def substitute(src, target, dic, make_notes=True):
  """
    Substitute strings in a template file and write the result to a target file, 
    taking care to place annotations that indicate where substitutions have been made.

    Input:
        src: Source file (template)
        target: Target file (generated)
        dic: Dictionary with substitutions

    Output:
        None (writes to target file)
    """

  from string import Template

  with open(src, 'r') as template:
    with open(target, 'w') as generated:
        for line in template:
            newline = Template(line).substitute(dic)
            if( newline != line):
                # Some string substitution happened
                if make_notes:
                    generated.write(annotate(newline))
                else:
                    generated.write(newline)
            else:
                # No substitution happened, just copy the line
                generated.write(line)

def annotate(line):
    """
    Annotate a line that has undergone substitution with comments indicating
    that a substitution has occurred. If the line contains multiple lines,
    each line except the first and last will be on its own line with a comment
    indicating the start and end of the substitution block.
    
    Input:
        line: The line to annotate
    Output: 
        newline: the annotated line.
    """ 
    split = line.splitlines()
    if(len(split)>1):
        newline = split[0] + ' !Hephaestos<<< \n'
        for k in range(1, len(split)-1):
            newline = newline + split[k] + '\n'
        newline = newline + split[-1] + ' !Hephaestos>>> \n'
    else:
        try:
            if('#' == split[0][0]):
                newline = split[0] + ' /* Hephaestos substitution */ \n' # Preprocessor directives have a different comment syntax.
            else:
                newline = split[0] + ' ! Hephaestos substitution \n'
        except IndexError:
            newline = line
    return newline