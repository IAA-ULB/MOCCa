# MOCCa

[![Github](https://img.shields.io/badge/github-MOCCa-blue?logo=github)](https://www.github.com/IAA-nuclear/tantalus_full)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Tantalus CI/CD Workflow](https://github.com/IAA-nuclear/tantalus_full/actions/workflows/workflow.yml/badge.svg?branch=master)](https://github.com/IAA-nuclear/tantalus_full/actions/workflows/workflow.yml)

**Nuclear Mean-field and Linear Response with Skyrme energy density functionals.**

---

## About MOCCa

MOCCa is a solver for the Skyrme mean-field and linear response equations which represents the single-particle wavefunctions on a three-dimensional coordinate space mesh. To save on both computational and human resources, the code is capable of operating in different symmetry modes, i.e. imposing different assumptions on the self-consistent symmetries of the nuclear mean-field many-body wavefunction. Despite this generality, an efficient numerical representation in terms of a [Lagrange mesh](https://www-sciencedirect-com.ezproxy.ulb.ac.be/science/article/pii/S0370157314004086) and clever [algorithms](https://link.springer.com/article/10.1140/epja/i2019-12766-6) make MOCCa sufficiently fast for [applications to the entire nuclear chart](https://link.springer.com/article/10.1140/epja/s10050-021-00642-1).

### Structure of the code and this repository 

MOCCa has evolved into a complex project: you will find that there is no one unique "source code" for the actual solver that covers all cases. Instead, a collection of Fortran90 source code files in `src_orig` gets modified by a complex preprocessor code, `Hephaestos.py` and its modules in `src_heph/`. Depending on your choice for 

- **Functional type**: LO, NLO, BXL, etc.
- **Self-consistent symmetries**: maximally symmetric, broken parity, broken time-reversal, etc.
- **Calculation type:**: static mean-field or linear response
- **System type**:: atomic nuclei or dense matter with periodic boundary conditions

Hephaestos will generate a customized MOCCa version tailored to your needs. For more information, see the [compilation section](docs/user-guide/compilation.md) in the documentation. 


### Some history 

MOCCa is a direct descendant of the Skyrme mean-field codes developed by the Brussels - Orsay - Saclay collaboration of P.-H. Heenen, H. Flocard and P. Bonche in the '80s. The most widely known of these codes is EV8 ([GitHub](https://github.com/wryssens/ev8),[original paper](https://doi.org/10.1016/j.cpc.2005.05.001)).
Today, MOCCa is capable of essentially all calculations that EV8 and its (unpublished) variants were designed for but 
supercedes its predecessors greatly in efficiency and generality. Nevertheless, the overall design and concept remains true to the original and we remain indebted to the original developers for their vision and many years of guidance.

In fact, MOCCa started as the PhD project of W. Ryssens under the guidance of P.-H.-Heenen. At the time, the goal was the study of rotational bands in a context of general symmetry breaking through cranking calculations. Coined at that time, the acronym MOCCa stands for **MO**dular **C**ranking **C**ode, with an extra 'a'. As always, our intentions evolved and so did MOCCa. Today the name is no longer truly suited to the tool; although still capable of cranking calculations, those are no longer the primary goal.

---

## Quickstart

### Prerequisites
- Python 3 with the numpy and scipy libraries, tested for Python 3.12.8
- a Fortran compiler (gfortran, ifort)
- BLAS and LAPACK (or SCALAPACK for MPI) installations
- HDF5 (optional, if you want to use it)

### Compilation
```bash
git clone https://github.com/IAA-nuclear/MOCCa.git
cd MOCCa
cp make_include/make.inc.gnu-serial make.inc
make mf 
```
> [!WARNING]
> Make sure that the `make.inc` file matches your system specifics!

The last command will 

 1. run Hephaestos for the default configuration file which is suitable for NLO EDFs,
 2. compile the resulting MOCCa source code with gfortran, 
 3. place the executable `MOCCa.default.exe` in the `exec/` directory.

Note that compilation may take several minutes, depending on your choices.

See the [compilation section](docs/user-guide/compilation.md) of the full documentation for advanced options.


### Running MOCCa

MOCCa takes its main user input from STDIN, e.g. you can type 

```bash
./exec/MOCCa.default.exe < MOCCa.in 
```

which will start a calculation. Note that essentially all relevant output will be printed to STDOUT, so you might want to save that.

A minimal example of contents for `MOCCa.in` would be: 

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

> [!WARNING]
> The explicit documentation of MOCCa is work-in-progress. If you don't find what you are looking for, check the comments in the source code.


The documentation is available in the `docs/` directory. To build it locally:

```bash
mkdocs serve
```
Then open the provided address in your browser. 

This requires a working [mkdocs](https://www.mkdocs.org/) installation. If you can't install that for some reasons: all documentation files are Markdown and so quite universally readable with your preferred application.

---

## Getting help

If the documentation is not sufficient to help you solve your problems, you can either

- Open an issue on GitHub; tagging one of the authors will help you get noticed.
- Contact Dr. W. Ryssens directly via email. You can find his adress easily on the internet. 

> [!WARNING]
> We welcome questions and will respond as time permits. However, we do not guarantee anything: our response times can be long and we might not provide a solution to any given issue, particularly if it concerns the implementation of new functionality.

---

## Contributors

- W. Ryssens [@wryssens](https://github.com/wryssens)
- M. Bender [@mbipnl](https://github.com/mbipnl)
- N. Shchechilin [@NikolaiNikolaevic](https://github.com/NikolaiNikolaevic)
- L. González-Miret Zaragoza [@luigonzar](https://github.com/luigonzar)
- P. Demol [@PepijnDemol](https://github.com/PepijnDemol)

---

## Citing MOCCa

When you use MOCCa in your research, we ask that you cite our work. The minimal requirement is to include the DOI of the version you used:

> W. Ryssens and M. Bender, *The MOCCa code*, [INSERT LINK]

Unfortunately, there is currently no dedicated peer-reviewed paper about MOCCa specifically. The closest reference is the original PhD thesis:

> W. Ryssens, *Symmetry breaking in nuclear mean-field models*, PhD Thesis, Université libre de Bruxelles (2016).

A digital copy is available [here](https://difusion.ulb.ac.be/vufind/Record/ULB-DIPOT:oai:dipot.ulb.ac.be:2013/235692/Holdings).

---

## License

This project is licensed under **AGPL v3** — see the [LICENSE](LICENSE) file for details.
