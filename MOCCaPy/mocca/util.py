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
    print()
    print(f" Started {func_name}".rjust(120, '-'))


def finished(func_name):
    print(f" Finished {func_name}".rjust(120, '-'))
    print()
