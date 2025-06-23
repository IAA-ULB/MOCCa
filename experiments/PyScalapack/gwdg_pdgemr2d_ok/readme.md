A `PyScalapack` attempt to same example as ../ScaLAPACK/gwdg_pdgemr2d_ok 

![7x5 matrix distributed in 2x2 blocks over 2x2 processes](../../../experiments/ScaLAPACK/image.png)

And this works too!

``` shell
> ./cl.sh 
context0: rank 0/4
(5, 7)
  C_CONTIGUOUS : False
  F_CONTIGUOUS : True
  OWNDATA : True
  WRITEABLE : True
  ALIGNED : True
  WRITEBACKIFCOPY : False

[[11. 12. 13. 14. 15. 16. 17.]
 [21. 22. 23. 24. 25. 26. 27.]
 [31. 32. 33. 34. 35. 36. 37.]
 [41. 42. 43. 44. 45. 46. 47.]
 [51. 52. 53. 54. 55. 56. 57.]]

context: rank 0/4 = (0,0): A_sub    |   context: rank 2/4 = (0,1): A_sub
[[0. 0. 0. 0.]                      |   [[0. 0. 0.]
 [0. 0. 0. 0.]                      |    [0. 0. 0.]
 [0. 0. 0. 0.]]                     |    [0. 0. 0.]]
context: rank 0/4 = (0,0): A_sub    |   context: rank 2/4 = (0,1): A_sub
[[11. 12. 15. 16.]                  |   [[13. 14. 17.]
 [21. 22. 25. 26.]                  |    [23. 24. 27.]
 [51. 52. 55. 56.]]                 |    [53. 54. 57.]]

context: rank 1/4 = (1,0): A_sub    |   context: rank 3/4 = (1,1): A_sub
[[0. 0. 0. 0.]                      |   [[0. 0. 0.]
 [0. 0. 0. 0.]]                     |    [0. 0. 0.]]
context: rank 1/4 = (1,0): A_sub    |   context: rank 3/4 = (1,1): A_sub
[[31. 32. 35. 36.]                  |   [[33. 34. 37.]
 [41. 42. 45. 46.]]                 |    [43. 44. 47.]]

context: rank 2/4 = (0,1): A_sub
[[0. 0. 0.]
 [0. 0. 0.]
 [0. 0. 0.]]
context: rank 2/4 = (0,1): A_sub
[[13. 14. 17.]
 [23. 24. 27.]
 [53. 54. 57.]]

context: rank 3/4 = (1,1): A_sub
[[0. 0. 0.]
 [0. 0. 0.]]
context: rank 3/4 = (1,1): A_sub
[[33. 34. 37.]
 [43. 44. 47.]]
```
The output is edited a bit for clarity.