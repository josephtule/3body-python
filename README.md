# 3body-python

The intent of this project is to create a 3-body simulator that will eventually generate an optimal trajectory from an orbit about the Earth to an orbit about the Moon, then generate an optimal control sequence in order to stay in the reference trajectory.

## Assumptions
To start, we will assume the following:

- The initial orbit Earth-Moon system is circular
- The Moon orbits around the Earth's equator
- There are no other bodies outside of the 3 body system that may affect the orbits


These assumptions can be removed at later phases with additional modeling.


## Equations of Motion
Using the Earth as the central body in an ECI reference frame, the equations of motion can be represented as the following:

$$\begin{equation}
\ddot{\vec{r}}_{ES} = -\frac{\mu_{E}}{||\vec{r}_{ES}||^3} \vec{r}_{ES} -\frac{\mu_{M}}{||\vec{r}_{MS}||^3} \vec{r}_{MS}
= -\frac{\mu_{E}}{||\vec{r}_{ES}||^3} \vec{r}_{ES} -\frac{\mu_{M}}{||\vec{r}_{ES} - \vec{r}_{EM}||^3} (\vec{r}_{ES} - \vec{r}_{EM})
\end{equation}$$ 

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
### Model Predictive Control


## References
https://epubs.siam.org/doi/pdf/10.1137/16M1062569
https://www.youtube.com/watch?v=yMJ_VU3jt7c
