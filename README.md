# MOCCa

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

**Mean-field and linear response code for nuclear physics using Skyrme energy density functionals.**

---

## About MOCCa

MOCCa is a FORTRAN code that performs mean-field and linear response calculations using Skyrme-type energy density functionals. Single-particle wavefunctions are represented on a Lagrange mesh.

**Hephaestos** is a Python code generator and preprocessor for MOCCa. Based on your choices for:
- Types of functionals
- Self-consistent symmetries
- Many-body calculation type (static mean-field or linear response)

Hephaestos generates a customized MOCCa version tailored to your specific needs.

MOCCa builds on early work by P.-H. Heenen, Flocard, and others.

---

## Quickstart

### Prerequisites
- Python 3 with numpy and scipy
- Fortran compiler (gfortran, ifort)
- BLAS and LAPACK (or SCALAPACK for MPI)
- HDF5 (optional, if you want to use it)

### Compilation
```bash
git clone https://github.com/IAA-nuclear/MOCCa.git
cd MOCCa
cp make_include/make.inc.gnu-serial make.inc
make mf 
```

The last command will 

 1. run Hephaestos for the default configuration file which is suitable for NLO EDFs,
 2. compile the resulting MOCCa source code with gfortran, 
 3. place the executable `MOCCa.default.exe` in the `exec/` directory.

See the [compilation section](docs/user-guide/compilation.md) of the full documentation for advanced options.


### Running MOCCa

MOCCa takes its main user input from STDIN, e.g. you can type 

```bash
./exec/MOCCa.BXL.exe < MOCCa.in 
```

which will start a calculation. Note that essentially all relevant output will be printed to STDOUT, so you might want to save that.

A minimal example of `MOCCa.in` would be the following: 

```
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy4'
/
&pairing
/
&evolution
maxiter=300
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
/
&IO
InputFilename='init'
OutputFilename='output.wf'
/
&MomentParam
/
&Cranking
/
```

Please see [this page in the documentation](docs/user-guide/structure.md) to start your journey into the many input options that MOCCa offers.

---

## Documentation

The documentation is available in the `docs/` directory. To build it locally:

```bash
mkdocs serve
```
Then open the provided address in your browser. 

This requires a working [mkdocs](https://www.mkdocs.org/) installation. If you can't install that for some reasons: all documentation files are Markdown and so quite universally readable with your preferred application.

---

## Getting help

- Open an issue on GitHub
- Contact: [your email]

---

## Contributors

- Dr. W. Ryssens [@wryssens](https://github.com/wryssens)
- Dr. M. Bender [@mbipnl](https://github.com/mbipnl)
- Dr. Nikolai Shchechilin [@NikolaiNikolaevic](https://github.com/NikolaiNikolaevic)
- Dr. Luis González-Miret Zaragoza [@luigonzar](https://github.com/luigonzar)
- Dr. Pepijn Demol [@PepijnDemol](https://github.com/PepijnDemol)

---

## Citing MOCCa

When you use MOCCa in your research, we ask that you cite our work. The minimal requirement is to include the DOI of the version you used:

> W. Ryssens and M. Bender, *The MOCCa code*, https://doi.org/10.1140/epja/s10050-021-00365-3

Unfortunately, there is currently no dedicated peer-reviewed paper about MOCCa specifically. The closest reference is the original PhD thesis:

> W. Ryssens, *Symmetry breaking in nuclear mean-field models*, PhD Thesis, Université libre de Bruxelles (2016).

A digital copy is available [here](https://difusion.ulb.ac.be/vufind/Record/ULB-DIPOT:oai:dipot.ulb.ac.be:2013/235692/Holdings).

---

## License

This project is licensed under **AGPL v3** — see the [LICENSE](LICENSE) file for details.
