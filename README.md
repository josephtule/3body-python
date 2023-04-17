# 3body-python

The intent of this project is to create a 3-body simulator that will eventually generate an optimal trajectory from an orbit about the Earth to an orbit about the Moon, then generate an optimal control sequence in order to stay in the reference trajectory.


## Equations of Motion
Using the Earth as the central body in an ECI reference frame, the equations of motion can be represented as the following:

$$ \ddot{\vec{r}}_{ES} = -\frac{\mu_{E}}{||\vec{r}_{ES}||^3} \vec{r}_{ES} -\frac{\mu_{M}}{||\vec{r}_{MS}||^3} \vec{r}_{MS} $$

$$ = -\frac{\mu_{E}}{||\vec{r}_{ES}||^3} \vec{r}_{ES} -\frac{\mu_{M}}{||\vec{r}_{ES} - \vec{r}_{EM}||^3} (\vec{r}_{ES} - \vec{r}_{EM})$$ 

and 

$$ \ddot{\vec{r}}_{EM} = -\frac{\mu_{E}}{||\vec{r}_{EM}||^3} \vec{r}_{EM} $$ 


Where the subscripts represent a vector from the first element (E - Earth, M - Moon, S - Spacecraft) to the second element.

## Perturbations

### Gravitational
### Aerodynamic
### 


## Numerical Integration
### RK4
### RK45

## Optimial Transfer Trajectory

## Optimal Control
