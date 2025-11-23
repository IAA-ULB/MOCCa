import time
from contextlib import ContextDecorator
from dataclasses import dataclass, field
from typing import Any, Callable, ClassVar, Dict, Optional
from tabulate import tabulate

class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

@dataclass
class Timer(ContextDecorator):
    """Time your code using a class, context manager, or decorator"""

    timers = {}

    def __init__(self, name = ''):
        """Initialize the Timer object"""
        self.name = name
        self._start_time = None
        self.timers[name] = [.0, 0]

    def start(self) -> None:
        """Start the timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self) -> float:
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        # Calculate elapsed time
        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        # Report elapsed time
        self.timers[self.name][0] += elapsed_time
        self.timers[self.name][1] += 1

        return elapsed_time

    def __enter__(self) -> "Timer":
        """Start a new timer as a context manager"""
        self.start()
        return self

    def __exit__(self, *exc_info: Any) -> None:
        """Stop the context manager timer"""
        self.stop()

    @classmethod
    def report(cls):
        print("Timers")
        table = [ [key,*value] for key,value in cls.timers.items()]
        print(tabulate(table, tablefmt="fancy_grid", headers=['name','time [s]','count']))
