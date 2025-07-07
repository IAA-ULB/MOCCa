# Scaling the system size n

approach:

- fix the number of nodes
- choose n, the size of the nxn system matrix
- increase n by a factor of 1.4 (roughly sqrt(2)) until the matrix is too large to be accomodated 

on 16x16 n=59217 triggers an allocation error, 42298 passes
on 32x32 n=185356 triggers an allocation error,

| processor grid |     n | walltime [s] |
|:--------------:| -----:| ------------:|
| 16x16=256      |  2048 |     410 |
| 16x16=256      |  2867 |     549 |
| 16x16=256      |  4013 |     855 |
| 16x16=256      |  5618 |    2554 |
| 16x16=256      |  7865 |    5283 |
| 16x16=256      | 11011 |   42031 |
| 16x16=256      | 15415 |   82902 |
| 16x16=256      | 21581 |  201178 |
| 16x16=256      | 30213 |  508036 |
| 16x16=256      | 42298 | 1138075 |
| 32x32=1024     | 

![16x16=256](slurm.out/16x16%20processes.png)

On 32x32 n=65534 fails, but n=46340 passes.

This probably means that n=46340 is close to the limit. The critical issue is probably that on rank 0, the full matrix is distributed, which consumes most of the memory.
we must adapt our test script to not allocate the ful matrix on node 0.
That work is executed in experiments/PyScalapack/pdsyev