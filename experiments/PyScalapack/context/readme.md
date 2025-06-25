# creating context with a variable rather than constant?

with `mpirun -n 4 python script.py` this script seems to work
```python
with scalapack(context_order, 2, 2) as context:
    print(f"({context.rank.value},{context.size.value})")
```

but this does not:

```python
nprow = npcol = 2
with scalapack(context_order, nprow, npcol) as context:
    print(f"({context.rank.value},{context.size.value})")
```

Both work correctly! what did not work is this:
```python
with scalapack(context_order, 1, 1) as context0:
    if context0: # this test sets nprow and npcol on rank 0 only!
                 # as a consequence `with scalapack(context_order, nprow, nprow) as context:` 
                 # will fail on all ranks but rank 0! MISTAKES ARE ALWAYS SILLY!

        nranks = context0.size.value
        if nranks == 1:
            nprow = npcol = 1
        elif nranks == 4:
            nprow = npcol = 2
        else:
            raise ValueError("expecting a power of 2 for the number of processes")
        print(f"{context0.rank.value}/{context0.size.value} {nprow=}, {npcol=}")
print(f"{nprow=}, {npcol=}")

with scalapack(context_order, 2, 2) as context: # distributed system on all ranks
    print(f"({context.rank.value},{context.size.value})")
```