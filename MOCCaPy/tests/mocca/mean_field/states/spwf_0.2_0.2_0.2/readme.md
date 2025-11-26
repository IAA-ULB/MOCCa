### contents of this directory

On a fully reduced mesh with M=30, d=0.8:

```python
    mesh = LagrangeMesh(M=30, d=.8, reduced=True)
```
with  single particle wave functions for $\mathrm{Ca}^{40}$ initialized as:
```python
	n_neutrons=20, n_protons=20,  
	nwn, nwp = 15, 15  
	osc_freq = (0.2, 0.2, 0.2)  
	hfpsi = states.HFPsi(  
	    n_neutrons=n_neutrons, n_protons=n_protons,  
	    n_proton_wf=nwp, n_neutron_wf=nwn,  
	    mesh=mesh,  
	    init='nilsson',osc_freq=osc_freq,  
		# orthogonalize=True, normalize=True
	)	
```
Following output is extracted from `MOCCa`
#### `nilsson_wfs.txt`
The single particle wave functions as delivered by "nilsson".  $15 \times 15 \times 15$ rows, $3 + 4 \times 30$ columns
```
hfpsi.mesh.xGrid[:] hfpsi.mesh.yGrid[:] hfpsi.mesh.zGrid[:] hfpsi.data[:,:]
```

#### `nillson_nabla_<x|y|x>_wfs.txt`
first order derivative of the single particle wave functions as delivered by "nilsson". Same format as `nilsson_wfs.txt`.

#### `nillson_laplacian_wfs.txt`
the Laplacian of the single particle wave functions as delivered by "nilsson". Same format as `nilsson_wfs.txt`.

#### `after-ortho_*.txt`
the same as `nilsson_*.txt` but after orthonormalization the spwfs as delivered by "nilsson". Same format as `nilsson_wfs.txt`.

#### `lagx.txt`, `lagy.txt`, `lagz.txt`
The Lagrange matrix for $\frac{d}{dx}$, $\frac{d}{dy}$, resp. $\frac{d}{dz}$,   on a reduced mesh:


