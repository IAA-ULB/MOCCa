

def started_finished(func):
    """A decorator that prints a message with the function name
    before and after the function is called.
    """
    def wrapper():
        print()
        print( f" Started {func.__name__}".rjust(120, '-'))

        func()

        print(f" Finished {func.__name__}".rjust(120, '-'))
        print()

    return wrapper


# The functions below can be used when the function uses a fixture.
# The starting_finished decorator inhibits the resolution of fixtures.

def started(func_name):
    """
    Call from within a function as `started(inspect.stack()[0][3])`.
    Requires 'import inspect`
    """
    print()
    print(f" Started {func_name}".rjust(120, '-'))


def finished(func_name):
    """
    Call from within a function as `finished(inspect.stack()[0][3])`.
    Requires 'import inspect`
    """
    print(f" Finished {func_name}".rjust(120, '-'))
    print()
