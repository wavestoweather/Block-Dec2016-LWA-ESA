## Python Library: hn2016_falwa (v0.5.0)

Compute from gridded climate data the Finite-amplitude Local Wave Activity (FALWA) and flux terms introduced in:

- [Huang and Nakamura (2016, JAS)](http://dx.doi.org/10.1175/JAS-D-15-0194.1)
- [Huang and Nakamura (2017, GRL)](http://onlinelibrary.wiley.com/doi/10.1002/2017GL073760/full).
- [Nakamura and Huang (2018, Science)](https://doi.org/10.1126/science.aat0721) *Atmospheric Blocking as a Traffic Jam in the Jet Stream*.

## Package Installation

This current version works in both Python 2.7 and 3.6. Note that since v0.3.0, some functions have backend in Fortran. You will need a fortran compiler (e.g. [gfortran](http://hpc.sourceforge.net/)) to implement the installation.

Since the package is still being actively developed, please use the *develop* mode for installation:
```
git clone https://github.com/csyhuang/hn2016_falwa.git
cd hn2016_falwa
python setup.py develop
```

To incorporate updates, pull the new version of the code from GitHub by:
```
git pull
```

If there are updates in the Fortran modules in the new commits, please re-compile them:
```
python setup develop -u
python setup develop
pytest # to check if the package can be run properly
```

## Quick start

The jupyter notebook in `examples/nh2018_science` demonstrates how to compute wave activity and reference states presented in Nakamura and Huang (2018). 
To make sure the package is properly installed in your environment, run through the notebook after installation to see if there is error.

THe conda environment for running the notebook can be found in `environment.yml`. To create the conda environment, execute:

```bash
conda env create -f environment.yml
```

## Interfaces

There are two interfaces for this library. One is the **developer interface**; the other is the **object-oriented 
interface**, which is a wrapper for the basis functions in the developer interface and also compiled fortran modules. 
Users are strongly adviced to use only the object-oriented interface.

### Object-oriented interface

The **object-oriented interface** is an easy-to-use interface that takes in the climate field and coordinates as the attributes of an object, and implement the wrapper functions above as methods.

There are two classes with object-oriented interface: *QGField* and *BarotropicField*. Please refer to the example/ directory:

Sample Script | Description
------------- | -------------
nh2018_science/demo_script_for_nh2018.ipynb | Compute wave activity and flux terms in the QG framework presented in Nakamura and Huang (2018, Science). Sample data can be retrieved with `download_example.py` in the same directory.
simple/oopinterface_example_BarotropicField.ipynb | It reads in a sample datasets "barotropic_vorticity.nc", which contains absolute vorticity field snapsnots from a barotropic decay model (Held and Phillips 1987). It computes both the equivalent-latitude relationship (e.g. Nakamura 1996) and local wave activity (Huang and Nakamura 2016) in a global domain.


### Developer Interface

The **developer interface**  contains separate functions that users can alter the inputs more flexibly. Functions 
are added upon users' request on new functionalities to test hypotheses (also see the *test* branch). The 
**developer interface** consists of 4 types of functions:  

- The **basis functions** are smallest unit of functions that make up the **wrapper functions** and **object-oriented interface**.  

- The **wrapper functions** implement particular analysis tasks for published work/manuscripts in preparation  

- The **utility functions** compute general quantities, such as static stability or quasi-geostrophic potential vorticity that are not specific to the finite-amplitude wave theory.   


## Inquiries / Issues reporting

- If you are interested in using the package, please leave your contact [here](https://goo.gl/forms/5L8fv0mUordugq6v2) or email me(csyhuang at uchicago.edu) such that I can keep you updated of any changes made.
- If you encounter *coding issues/bugs* when using the package, please create an [Issue ticket](https://github.com/csyhuang/hn2016_falwa/issues).
- If you have scientific questions, please contact Clare S. Y. Huang via email(csyhuang at uchicago.edu).