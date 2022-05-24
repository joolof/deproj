# Deproj

A simple tool to deproject direct imaging observations of debris disks.

In principle it could work for any kind of disks that you are interested in, but since the disk is assumed to be vertically flat, it should be more accurate for debris disks rather than proto-planetary disks.


## Installation

Simply clone the repository and run

```python
python3 setup.py develop
```

and you can then import the class with

```python
from deproj import Deproj
```

## First example

There are very few parameters that you need to provide to get the deprojected map. Here's a quick use example

```python
disk = Deproj('data_example/HR4796_Qphi_400.fits', nr = 30, nt = 60, pixscale = 0.0072)
```

You need to provide at least the name of the fits file with the observations, and then there are the following optional parameters that you can pass:

```python
nr = 30           # Number of radial bins
nt = 60           # Number of azimuthal bins
pixscale = 0.0072 # Size of one pixel in arcsec
```

Afterwards, you can call the main method as

```python
disk.go(amin = 0.7, amax = 1.3, incl = 77.72, pa = -151.59)
```

where the parameters are

```python
amin = 0.7        # The inner radius of the mask
amax = 1.3        # The outer radius of the mask
incl = 77.72      # The inclination of the disk
pa = -151.59      # The position angle of the disk
```

and this will produce the feollowing output:

![HR4796](screenshots/HR4796.png)

To help you better visualize how this works you can also call the following method

```
disk.debug(amin = 0.7, amax = 1.3, incl = 77.72, pa = -151.59)
```

which will produce the output below

![debug](screenshots/debug.png)

where you can see on the left the distance in the midplane, with the contours showing the inner and outer radii of the mask, while the right panel shows the azimuthal angle in the midplane.


## To be done

- [x] Install script
- [ ] Additional parameters for the plotting, `vmin`, `vmax`, `xlim`, `cmap` etc
- [ ] Plotting options, only `polar`, only `cartesian`, `full`, etc
- [ ] Provide either a file name, or a 2D frame directly
- [ ] Automatic determination of the inclination and position angle? It would increase the dependencies required though.


