import time
import sys
from contextlib import ContextDecorator
from dataclasses import dataclass, field
from typing import Any, Callable, ClassVar, Dict, Optional
from tabulate import tabulate

class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""


class Timer(ContextDecorator):
    """Time your code using a Timer instance, a context manager, or a decorator.

    ```python
        from mocca.util.timer import Timer

        @Timer('myfun') # Create Timer 'myfun' for timing the function using as a decorator
        def myfun():
            # code
            tmr = Timer('snippet') #create Timer 'snippet' using a Timer instance
            tmr.start()
            # code snippet to be timed
            tmr.stop()

            with Timer('snippet2'): # Use Timer as q context manager
                # code snippet 2 to be timed

            # Timer 'snippet2' is now stopped.
    ```
    """

    timers = {}

    def __init__(self, name = ''):
        """Initialize the Timer object"""
        self.name = name
        self._start_time = None
        self.timers[name] = [ 0, .0, .0, sys.float_info.max, .0 ]
        #   [count, sum, ssq, min, max]

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
        rec = self.timers[self.name]
        rec[0] += 1
        rec[1] += elapsed_time
        rec[2] += elapsed_time**2
        rec[3] = min(elapsed_time, rec[3])
        rec[4] = max(elapsed_time, rec[4])

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
        print("Timers [s}")
        table = []
        for name,value in cls.timers.items():
            count = value[0]
            total = value[1]
            mean = total/count
            stddev = (value[2] - total*total/count) / (count-1) if (count > 1) else None
            mn = value[3]
            mx = value[4]
            table.append([name, count, total, mean, stddev, mn, mx])
        table.sort(key=lambda x: x[2], reverse=True)
        print(tabulate(table, tablefmt="fancy_grid", headers=['name', 'count', 'total', 'mean', 'stddev', 'min', 'max' ]))
