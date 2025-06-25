# Scaling tests

We solve a rather small problem sequentially and distributed with different blocking factors.
The matrix is represented fully on rank 0 and distributed over all ranks. 
we restrict to a single node 2x2, 4x4 and 8x8 processor grids

a Vaughan node has 256 GB RAM, for 64 cores, that is 4 GB per core. 
1 GB is `1073741824` bytes, so 1 core can accomodate `4 * 1 073 741 824 / 8 = 536 870 912` doubles.

this means 

| n     | np x np | #DP values per proces |
|-------|---------|-----------------------|
| 512   | 2x2 =4  |  