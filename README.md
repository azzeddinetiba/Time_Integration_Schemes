# Time Integration Schemes
This is a finite element research project, where the goal was to study the damping effect of different numerical integration schemes for the time domain in a transient analysis, for beam elements. Especially for an impact load on the longitudinal degrees of freedom. The main studied integration schemes were Backward-Euler, the Newmark Method,  Centered scheme, some Runge-kutta schemes and the generalized alpha scheme.

The carried simulations were compared to the commercial code Abaqus, and different properties were investigated, e.g accuracy, stability, robustness ...
These methods and the main finite element method were implemented and programmed in Matlab.
A graphical user interface was also programmed to allow the user to choose easily the numerical integration method and view the displacement results.

[![Product Name Screen Shot][product-screenshot]]()


The GUI show the results of an 3500 N impact simulation on the end of a 2m beam, for a time period of 1e-4 seconds. 
Theses parameters can be changed on the 
```matlab
Newmark.m
```

and
```matlab
dynamics.m
```



ENJOY !






[product-screenshot]: GUI_screen.PNG


