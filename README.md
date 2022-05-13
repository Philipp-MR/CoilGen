<div id="top"></div>


 <img src="./Documentation/GoilGen_Logo.png" width="300">

<!-- ABOUT THE PROJECT -->

The CoilGen Project is a community-based tool for the generation of coil Layouts within the MRI/NMR environment. It is based on a boundary element method and generates an interconnected non-overlapping wire-tracks on 3D support structures. The focus of this work is post processing.

The user must specify a target field (e.g., bz(x,y,z)=y for a constant gradient in the y-direction) and a surface mesh geometry (e.g., a cylinder defined in an .stl file). The code then generates a coil layout in the form of a non-overlapping, interconnected wire trace to achieve the desired field.

Up to now, the code is written in MATLAB, but future migration to python might be advantageous, especially since it does not need proprietary software licenses. The author is very willing to collaborate with anyone who wants do the translation.



<!-- GETTING STARTED -->
## Getting Started

Check the documentation to get started



### Installation

1. Download and extract the file of the CoilGen repository
2. Run one of the examples in the folder "Examples"

The project requires MATLAB and optionally FastHenry2 for calculation of the coil inductance.  The MATLAB version should not be older than 2020A. Important: It also requires the MATLAB MAPPING-toolbox since the functions polyxpoly and inpolyon are used. Non-proprietary versions of these functions are also very welcome.


### Exemplary Images

<!-- ![plot](./Documentation/flow_chart_algorithm_revised.png)  -->
![plot](./Documentation/Results_CoilGen_YGradient.png)





<!-- LICENSE -->
## License

 See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Philipp Amrein

Project Link: [https://github.com/Philipp-MR/CoilGen]

<p align="right">(<a href="#top">back to top</a>)</p>


## Citation

For citation of this work, please refer to the following publication:
https://onlinelibrary.wiley.com/doi/10.1002/mrm.29294
https://doi.org/10.1002/mrm.29294

<p align="right">(<a href="#top">back to top</a>)</p>
