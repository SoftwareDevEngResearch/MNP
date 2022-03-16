# MNP
Package for handling Multiphysics Neutronics systems. Requires numpy and numba.
Run the software with  

>$python wrapper.py material.npz input.in

material.npz specifies the material parameter used for the simulation. This .npz file has the following numpy arrays:
<ol>
<li> E,             Energy grid boundaries,         [1xG+1] </li>
<li> SigmaT,        Total cross section,            [1xG] </li>
<li> SigmaF,        Fission cross section,          [1xG] </li>
<li> SigmaS,        Scattering cross section        [GxG] </li>
<li> v,             Speed in energy group           [1xG] </li>
<li> nu_prompt,     Neutrons per prompt fission,    [1xG] </li>
<li> nu_delayed,    Neutrons per delayed fission,   [1xG] </li>
<li> chi_prompt,    Probability of prompt fission,  [GxG] </li>
<li> chi_delayed    Probability of delayed fission, [GxJ] </li>
<li> decay,         Decay constant for precursor    [1xJ] </li>
<li> beta_frac,     Probability of prompt fission,  [1xJ] </li>
</ol>

Input.in specifies simulation parameters. It is self documenting.

Two example problems have been included, a subcritical problem and a delayed supercritical problem. These can be run with the previously mentioned command and varying the input file between the two choices. A rudementary plotting script has been included in the Outputs directory to visualize the output on a very basic level. 
