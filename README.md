# Deproj

A simple tool to deproject direct imaging observations of debris disks.

In principle it could work for any kind of disks that you are interested in, but since the disk is assumed to be vertically flat, it should be more accurate for debris disks rather than proto-planetary disks.


## Installation

Simply clone the repository and run

```
python3 setup.py develop
```

and you can then import the class with

```
from deproj import Deproj
```

![HR4796](screenshots/HR4796.png)

Work in progress, more later on

## To be done

- [x] Install script
- [ ] Additional parameters for the plotting, `vmin`, `vmax`, `xlim`, `cmap` etc
- [ ] Plotting options, only `polar`, only `cartesian`, `full`, etc
- [ ] Provide either a file name, or a 2D frame directly
- [ ] Automatic determination of the inclination and position angle? It would increase the dependencies required though.


