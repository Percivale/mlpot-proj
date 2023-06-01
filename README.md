# A Machine Learning Potential for Amorphous Alumina in LAMMPS

# The Matsui Potential
This potential was developed by M. Matsui and published in 1994, "A transferable interatomic potential model for crystals and melts in the system CaO-MgO-$Al_2O_3$-$SiO_2$. The potential is,

$
\begin{equation}
    V(\Vec{r_{ij}}) = \frac{q_i q_j}{\Vec{r_{ij}}} - 
\end{equation}
$

Maybe plot the potential here...

# Amorphous aluminum oxide
The structures were simulated in LAMMPS with the born/coul/long potential which is of the same form as the Matsui potential. When producing the structures the system starts out in an initial state that is heated far over the melting point of the material. This results in the atoms experiencing large velocities which originally resulted in clusters of Aluminum atoms due to the fact that they moved so fast that the potential had no time to repulse the atoms. To bypass this issue, an additional strong short term repulsive potential was added to make sure that no unphysical clustering of atoms would take place. The term scaled with $\frac{1}{r^{24}}$. 

Plot potential with 1/r24 term...

When the system has turned liquid, the material is cooled down to achieve the amorphous structure. If the system was cooled down slow enough it would have time to rearrange itself into the crystal. The faster the cooling is, the more unordered the structure should end up as. 

Radial distribution functions from the simulations....
-,,- around melting point....

# Neural Networks
The hope is that the machine learning models can be used to decrease the computational requirements for molecular dynamics simulations. 
Downsides: 
The machine learning models will only be trained on amorphous alumina oxide, which means that it will only be able to do predictions on these sorts of structures. 
## Linear Layers
## Convolutional Layers
## Activation Function
## Loss Function
## Optimization


