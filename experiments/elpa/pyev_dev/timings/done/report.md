# Scaling tests

## Software

ELPA comes with a python wrapper `pyelpa`. For ScaLAPACK a python wrapper was found on pypi and forked inside tantalus (`pyscalapack`). The latter was wrapped a second time (`pyev`) to provide the same interface as `pyelpa`. This allows for the selection of a backend (ELPA or ScaLAPACK):

```python
    from pyev.helpers import use_pyelpa
    if use_pyelpa:
        # Use pyelpa directly
        from pyelpa import ProcessorLayout, DistributedMatrix, Elpa
    else:
        # Use the same interface as pyelpa, with Scalapack as a backend
        from pyev   import ProcessorLayout, DistributedMatrix, Elpa
    # Application code below is independent of the chosen backend.
    ...
```

## Scaling tests

Run on LUMI-C where compute nodes have 128 cores and approximately 2 GB of RAM per core.

### Single node tests

![](png/single%20node%20scaling%20(scalapack).png)

![](png/single%20node%20scaling%20(elpa).png)

![](png/single%20node%20scaling%20(scalapack%20vs%20elpa).png)

![](png/single%20node%20scaling%20(scalapack%20over%20elpa%20ratio).png)

From the figure above it is already clear that for systems with `na > 30000` ELPA is significantly faster than ScaLAPACK.

### Multi-node tests


![](png/multi-node%20scaling%20(scalapack).png)

![](png/multi-node%20scaling%20(elpa).png)

![](png/multi-node%20scaling%20(scalapack%20over%20elpa%20ratio).png)

The advantage of elpa ranges from 2.5 times faster than ScaLAPACK for small systems to between 12.5 and 20 times faster for larger systems (na >= 50000).

### Strong scaling

![](png/strong%20scaling%20(elpa).png)

Strong scaling is obviously missing, because the work of the algorithm increases as the third power of the number of rows.

> :memo: **Note:** With approximately 2 GB RAM per core we estimate the upper bound for the local matrix size as `math.sqrt(2*1024**3) ~= 46 000`.