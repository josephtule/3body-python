# 3body-python

The intent of this project is to create a 3-body simulator that will eventually generate an optimal trajectory from an orbit about the Earth to an orbit about the Moon, then generate an optimal control sequence in order to stay in the reference trajectory.

## Assumptions
To start, we will assume the following:

- The initial orbit Earth-Moon system is circular
- The Moon orbits around the Earth's equator
- There are no other bodies outside of the 3 body system that may affect the orbits


These assumptions can be removed at later phases with additional modeling.


## Equations of Motion
Using an arbitrary central body (C) with an inertial reference frame as well as an arbitrary number of external "planets" (O), the equations of motion can be represented as the following:

$$\begin{equation}
\ddot{\vec{r}}_{CS,i} = -\frac{\mu_{C}}{||\vec{r}_{CS,i}||^3} \vec{r}_{CS,i} - \sum_{j} \frac{\mu_{O,j}}{||\vec{r}_{OS,i,j}||^3} \vec{r}_{OS,i,j}
= -\frac{\mu_{C}}{||\vec{r}_{CS,i}||^3} \vec{r}_{CS,i} - \sum_{j} \frac{\mu_{O,j}}{||\vec{r}_{CS,i} - \vec{r}_{CM,j}||^3} (\vec{r}_{CS,i} - \vec{r}_{CM,j})
\end{equation}$$ 

Where i is the ith body of interest and j is the jth external planet. This equation will be the governing equation for each body other than the central body (including satellites and all "planets"). Since satellites are low in mass, we will ignore their pull on the other bodies but the external planets will affect each other. 

### State Space

To allow for trajectory optimization, we will take to account the position of the moon at all times, this will result in size 12 state vector, with the state being:

$$
\vec{x} =  \left\lbrack \begin{array}{c}
\vec{r}_{CS,i} \\
\dot{\vec{r}}_{CS,i} \\
\end{array}\right\rbrack = 
\left\lbrack \begin{array}{c}
\vec{x}_{1:3} \\ 
\vec{x}_{4:6} \\ 
\end{array}\right\rbrack
$$

After the optimal trajectory has been generated, the state can be reduced to a state vector of size 6, and including the dynamics for the moon along with the generated trajectory since the orbit of the moon is assumed not to be affected by the spacecraft.

Taking the derivative of the state vector yields

$$
\dot{\vec{x}}_i = 
\left\lbrack \begin{array}{c}
\dot{\vec{x}}_{i,1:3} \\ 
\dot{\vec{x}}_{i,4:6} \\ 
\end{array}\right\rbrack$$ 

$$ =
\left\lbrack \begin{array}{c}
\vec{x}_{4:6}\\ 
-\frac{\mu_C}{(\vec{x}_{i,1:3}^T\vec{x}_{i,1:3})^{3/2}} \vec{x}_{i,1:3} - \frac{\mu_{O,j}}{(\vec{x}_{i,1:3}-\vec{x}_{i,4:6})^T(\vec{x}_{1:3}-\vec{x}_{4:6})^{3/2}} (\vec{x}_{1:3}-  \vec{x}_{4:6}) + u_{1:3}\\ 
\end{array}\right\rbrack
$$

We add $u_{1:3}$ here to imply a thrust control on the satellite.

## Perturbations
In the generation of the optimimal trajectory, perturbations will be ignored, but in the optimal control problem, external perturnbations will be included in order to work the controller. In this project, only gravitational perturbations and aerodynamic perturbations will be considered. 
### Gravitational
Gravitational perturbations come in the form of 3rd body forces (in a 2 body problem) and the geometry of the planetary bodies that a spacecraft may be orbiting.
3rd body forces will be generated through the dynamics of the system. The geometric perturbations will be represented by the $J_i$ perturbations, especially $J_2$ as the other coefficients will have much less of an effect on the spacecraft. 
### Aerodynamic
### 


## Numerical Integration
### RK4
### RK45

## Optimial Transfer Trajectory
### Objective Function

### Constraints
- System Dyanmics
- Path Constraints
- Boundary Constraints
- Path and Control Bounds
- Time and Time-Dependent State Bounds

## Optimal Control
### Model Predictive Control

## Future Plans
- Allow rendezvous with other spacecraft
- Allow central body changes
- Change code to reflect above changes

## References
https://epubs.siam.org/doi/pdf/10.1137/16M1062569  
https://www.youtube.com/watch?v=yMJ_VU3jt7c
