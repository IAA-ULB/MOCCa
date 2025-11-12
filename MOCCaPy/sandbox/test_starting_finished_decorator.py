def starting_finished(func):
    def wrapper():
        print(f"\nStarting {func.__name__}")
        func()
        print(f"finished {func.__name__}\n")
    return wrapper

@starting_finished
def test_fun():
    print("Hello world")