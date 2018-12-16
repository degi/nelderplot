# A Computer-aided Tool for Nelder Plot Design

## Overview

Nelder plot (Nelder 1962) is an experimental design for testing a multiple tree
spacings in a single plot. The shape of plot is circular with some number of concentric wheels or arcs radiating outward and the spokes connecting the center. The intersections of spokes and arcs are the tree positions. The tree position on different wheels represent a different the tree spacing distance or tree density. This design eliminates the need for separate experimental plots for each tree density. The difficulties with Nelder plot design is their uncommon layout and complex parameterization (Parrott et al. 2012). A computer-aided tool for designing the Nelder plot can be very useful for the practical implementation of this experimental design.  

## Software Specification

### Input

#### Practical Parameters

* Minimum planting density
* Maximum planting density
* Number of density range
* Rectangularity
  
#### Standard Parameters

* Plot dimension (Width, Height)
* Initial radius (r_0)
* Constant increment of arcs radius (α)
* Division angle of spokes (θ)

#### Additional Parameters

* Center of rotation (x, y)
* Number of Species

### Output

* Plot design
* Individual position coordinates and variables

### Features

* Standalone application
* Can be used in R package library
* Saving plot design output as image and data file

## Nelder Plot Parameters

Plot design type A1 require a parameter to define the initial radius (r_0) and a constant increment value (α). The wheels radius is then calculated with formula:

r_n = r_0 α^n

The individual growing space is calculated based on formula:

A_n=r_n^2 θ(α-α^(-1))/2

where θ is degrees between the wheel spokes (Nelder 1962).

## Nelder Plot Parameters for practitioner

The parameterization with practical application use the formula proposed by Parrot et al. (2012) to convert the variable into standard Nelder plot parameters. The calculation is as follow:  

α = exp((log(d_1) - log(d_2))/(2 * n - 2))

where α is constant increment, d_1 is minimum planting density, d_1 is maximum planting density, and n is number of density range.

θ = t * (α^0.5 - 1/α^0.5)

where θ is division angle of spokes and t is rectangularity.

r_0 = (20000/(d_1 θ (α**3 - α)))^0.5

where r_0 is initial wheel radius.

## References

* Nelder JA (1962) New kinds of systematic designs for spacing experiments. Biometrics 283–307
* Parrott DL, Brinks JS, Lhotka JM (2012) Designing Nelder wheel plots for tree density experiments. New Forests 43:245–254. doi: 10.1007/s11056-011-9278-4
