# A Computer-aided Tool for Nelder Plot Design

Author: Degi Harja Asmara and Steve Cumming

## Overview

Nelder plot (Nelder 1962) is an experimental design for testing a multiple tree
spacings in a single plot. The shape of plot is circular with some number of concentric wheels or arcs radiating outward and the spokes connecting the center. The intersections of spokes and arcs are the tree positions. The tree position on different wheels represent a different the tree spacing distance or tree density. This design eliminates the need for separate experimental plots for each tree density. The difficulties with Nelder plot design is their uncommon layout and complex parameterization (Parrott et al. 2012). This computer-aided tool for designing the Nelder plot can be very useful for the practical implementation of this experimental design.

The binary distribution of this software can be found at:
<https://sourceforge.net/projects/nelderplot.ecomodels.p/>

The source code location:
<https://github.com/degi/nelderplot>

## Software Specification

The software is designed to be used by either forestry practitioner or those who one to explore the Nelder plot parameterization. Two type of input parameters option is available: practical and standard parametrization. The plot design will be automatically displayed as soon as the required parameters are satisfied.  

## Installation

Download the source code, and go to the app root folder and type: 
```
python setup.py install
```
It will install the required files. After, run the application by typing:
```
python -m nelderplot.app
```
*NOTE: This app was developed using python 3*

### Input

#### Practical Parameters

* Minimum planting density (trees/ha)
* Maximum planting density (trees/ha)
* Number of density range
* Rectangularity (1 - for 'square' approximation)
  
#### Standard Parameters

* Plot dimension (Width, Height)
* Initial radius (r_0)
* Arcs radius rate of change (α)
* Division angle of spokes (θ)

#### Additional Parameters

* Center of rotation (x, y)
* Number of Species (max 3)

### Output

The software provides various output including:

* Text file with all the parameters, plot variables and data points
* Plot design that can be exported into various image format
* KML file which includes individual plant position coordinates and plot design information

Below are some screenshot of the software:

* Main interfaces of the software

![Main interfaces](screenshot/nelderplot.png?raw=true "Main interfaces")

* The option for half wheel Nelder by selecting the center of rotation

![Alt text](screenshot/nelderplot_half.png?raw=true "Main interfaces")

* Another option with a quarter wheel

![Alt text](screenshot/nelderplot_quarter.png?raw=true "Main interfaces")

* Standard parameterization may use a rectangle plot

![Alt text](screenshot/nelderplot_rec_plot.png?raw=true "Main interfaces")

* The rectangle plot with center of rotation on the plot border
  
![Alt text](screenshot/nelderplot_rec_plot_half.png?raw=true "Main interfaces")

* Option for number of species

![Alt text](screenshot/nelderplot_full_3sp.png?raw=true "Main interfaces")

* Two type of species arrangement: Alternate spokes and  alternate wheels combination

![Alt text](screenshot/nelderplot_full_3sp_altspokes.png?raw=true "Main interfaces")

* The variable information of Nelder plot wheels

![Alt text](screenshot/nelderplot_variables.png?raw=true "Main interfaces")

* Data point information for relative plant coordinates
  
![Alt text](screenshot/nelderplot_datapoints.png?raw=true "Main interfaces")

* The KML file can be imported to google map as map layer

![Alt text](screenshot/gmap1.png?raw=true "Main interfaces")

* The KML file includes the information of Nelder plot parametrization and design

![Alt text](screenshot/gmap3.png?raw=true "Main interfaces")

## Nelder Plot Parameters

Plot design type A1 require a parameter to define the initial radius (r_0) and a constant increment value (α). The wheels radius is then calculated with formula:

r_n = r_0 α^n

The individual growing space is calculated based on formula:

A_n=r_n^2 θ(α-α^(-1))/2

where θ is degrees between the wheel spokes (Nelder 1962).

### Nelder Plot Parameters for practitioner

The parameterization with practical application use the formula proposed by Parrot et al. (2012) to convert the variable into standard Nelder plot parameters. The calculation is as follow:  

α = exp((log(d_1) - log(d_2))/(2 * n - 2))

where α is constant increment, d_1 is minimum planting density, d_1 is maximum planting density, and n is number of density range.

θ = t * (α^0.5 - 1/α^0.5)

where θ is division angle of spokes and t is rectangularity.

r_0 = (20000/(d_1 θ (α^3 - α)))^0.5

where r_0 is initial wheel radius.

## References

* Nelder JA (1962) New kinds of systematic designs for spacing experiments. Biometrics 283–307
* Parrott DL, Brinks JS, Lhotka JM (2012) Designing Nelder wheel plots for tree density experiments. New Forests 43:245–254. doi: 10.1007/s11056-011-9278-4
