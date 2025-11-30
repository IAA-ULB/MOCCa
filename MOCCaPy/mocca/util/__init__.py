
def title_line(text='', width=80, char='*', start=4, above=False, below=False, cr=True):
    """
    Function for creating title lines.

    Args:
        text: the title text
        width: width of the title line
        char: character used for the line
        start: start point of the title text, a negative number indicates that the text be centered
        above: print a line above or not
        above: print a line below or not
        cr: add a carriage return at the end of the line

    Some examples::

        >>> print(title_line("A title", 30))
        *** A title ******************
        >>> print(title_line("A title",30,start=8 ,char='-',above=True,below=True))
        ------------------------------
        ------- A title --------------
        ------------------------------
        >>> print(title_line(width=30))
        ******************************
    """
    w = int(width/len(char))
    line0 = w*char
    if cr:
        line0 += '\n'
    if text:
        text = ' '+text+' '
    n =len(text)
    if start < 0:
        start = (width - len(text) - 2)//2
    line = line0[:start-1]+text+line0[n+start-1:]
    if above:
        line = line0+line
    if below:
        line = line+line0
    return line

def dbg_assert(d):
    """Call dbg_assert() for all objects in `d` which have this attribute. Typically called
    as:
    >>> dbg_assert(locals())

    Args:
        d: dict "variable_name : instance".
    """
    for instance in d.values():
        try:
            instance.dbg_assert()
        except AttributeError:
            pass # instance does not have `dbg_assert()` method.