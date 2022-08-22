# 🌤 Solecito

A collection of MATLAB functions for the detailed simulation of photovoltaic plants.

The project began in 2013 as an extension of PVLIB (Matlab), answering the need of independent consultant Hugo Capdevila for a PV simulation software that could account for the effects of electrical layout and complex geometry in shading and mismatch losses. Since then, several commercial tools have caught up to consider string layout information in some way or another. But our library still offers unique features that might be of value to the PV modelling community. Namely:

- View factor calculation for arbitrary tessellations of the sky dome, allowing the calculation of anisotropic diffuse shading. Our custom polygonal projection approach accounts also for the cosine-response of the modules for each diffuse irradiance component.

- Spherical interpolation of diffuse 3D view-factors, and adaptive resolution of direct-shading calculations, reducing the computational expense of run-time electrical simulations. 

- Performance optimized C++ library for One-Diode-Model current-voltage (IV) curve calculations, and series-parallel IV curve addition (by Huan Zhou).

- Functions for probabilistic sensor-fusion of irradiation data, including quality control, irradiance separation, and transposition (See [matlab-pvmeteo](https://github.com/martinherrerias/matlab-pvmeteo) and [radiance_RBE](https://github.com/martinherrerias/radiance_RBE)).

## Acknowledgements

Since 2019, the code has been reworked in the framework of the **HyForPV** project, funded by the German Federal Ministry of Economic Affairs and Climate Action on the basis of a decision by the German Bundestag. The project was comissioned and supported privately for years by Hugo Capdevila -- I trust he would be happy to see it bear fruit, even if not in the way he envisioned. Thanks to the financial support of the Mexican National Council of Science and Technology (CONACYT) and the German Academic Exchange Service (DAAD), in 2011-2013.

## References

[1] Capdevila, H., Herrerias A., M., Marola, A., 2014. Anisotropic Diffuse Shading Model for Sun-tracking Photovoltaic Systems. Energy Procedia, 2013 ISES Solar World Congress 57, 144–151. https://doi.org/10.1016/j.egypro.2014.10.018

[2] Capdevila, H., Marola, A., Herrerias A., M., 2013. High resolution shading modeling and performance simulation of sun-tracking photovoltaic systems. 9th International Conference on Concentrator Photovoltaic Syste,s: CPV-9, Miyazaki, Japan, pp. 201–204. https://doi.org/10.1063/1.4822231

[3] Herrerias A., M., Capdevila, H., 2018. Detailed Calculation of Electrical Mismatch Losses for Central- and String-Inverter Configurations on Utility-Scale PV Arrays, in: 35th European Photovoltaic Solar Energy Conference and Exhibition. pp. 1973–1978. https://doi.org/10.4229/35thEUPVSEC20182018-6CV.2.36

[4] Herrerias A., M., Capdevila, H., Zhou, H., Hammer, A., 2020. Simulation of Large PV Plants Using a Continuous Radiance Distribution Model and Cell-Resolution Mismatch Calculation. EU PVSEC 2020 37th European Photovoltaic Solar Energy Conference and Exhibition, online, pp. 1311–1316.

[5] Zhou, H., Niethammer, C., Herrerias A., M., 2021. Usage Experiences of Performance Tools for Modern C++Code Analysis and Optimization, in: Tools for High Performance Computing 2018 / 2019. Springer International Publishing, Cham, pp. 103–121. https://doi.org/10.1007/978-3-030-66057-4_5

[6] Herrerias A., M., Hammer, A., 2021. Probabilistic Estimation of Irradiance Components: A State Space Framework to Integrate Empirical Models, Forecasts, and Imperfect Sensor Measurements, ISES SWC 2021. Under review, [see preprint](https://www.researchgate.net/publication/356491420_Probabilistic_Estimation_of_Irradiance_Components_A_State_Space_Framework_to_Integrate_Empirical_Models_Forecasts_and_Imperfect_Sensor_Measurements).

## Use

#### 0. Try running the `startup.m` function from Matlab. It should add all the required folders to the path, and compile all libraries. If it worked, congratulations! jump to step 4... otherwise grab a coffee, and go to step 1.

#### 1. Clone the repository (*with submodules!*):

`git clone --recurse-submodules git@github.com:martinherrerias/solecito.git`

#### 2. Add code to the matlab path (note some folders are removed, to avoid path cluttering and name conflicts):

```
BASEDIR = '~/solecito';  % your cloned repository
cd(BASEDIR);
addpath(genpath('.'));
rmpath(genpath('./.git'));
rmpath(genpath('./tests/ODM'));

rmpath(genpath('./pvCplusplus'));
addpath('./pvCplusplus/mexMdpwl');

rmpath('./matlab-pvmeteo/PVLIB/Example Data');
rmpath('./external/clipper/cpp');
rmpath(genpath('./external/jigsaw/jigsaw-matlab'));
```

#### 3. Compile CPP libraries

- Compile C++ One-Diode-Model and Piece-Wise-Linear electrical solver
```
cd(fullfile(BASEDIR, 'pvCplusplus/mexMdpw'));
commonBuiltScript()
addSeriesBuiltScript()
addParallelBuildScript()
ODMpwlapproxBuildScript()
```

- Compile CLIPPER polygon-clipping library/mex wrapper (see `/external/clipper/matlab/README.txt`)
```
cd(fullfile(BASEDIR,'external/clipper/matlab'));
cd('./external/clipper/matlab');
mex '-D__int64=__int64_t' '-I../cpp' '../cpp/clipper.cpp' 'mexclipper.cpp'
```

- Compile JIGSAW mesh refining utility (see `/external/jigsaw/jigsaw-matlab/README.md`)
```
cd(fullfile(BASEDIR,'external/jigsaw'));
initjigsaw()
compile()
endjigsaw()
```

#### 4. (Optional) add user accounts

User accounts are required to access [SoDa](http://www.soda-pro.com/) and [geonames](https://www.geonames.org/). In most cases these provide optional features (look for the option `useweb = false`). To use the services, modify the file `./Remote/useraccounts.txt` with your credentials and save as `useraccounts.m`:

```
copyfile('./Remote/useraccounts.txt','./Remote/useraccounts.m');
% (update `useraccounts.m` with your credentials)
```

To avoid repeated web requests (or to allow manual downloads), _mcclear*.csv_ and _merra2*.csv_ files can also be retrieved from local folders in the search-path. Files are identified by name (see `getremote.searchexisting`). The functions `merra2_healthcheck` and `mcclear_healthcheck` can be used periodically to auto-rename the files and speed-up the search:

```
addpath('~/Documents/CAMS/MERRA2');
copyfile('./Remote/resources/merra2_healthcheck.m','~/Documents/CAMS/MERRA2');

addpath('~/Documents/CAMS/McClear');
copyfile('./Remote/resources/mcclear_healthcheck.m','~/Documents/CAMS/McClear');
```

#### 5. Where to start?

The program was originally designed to be a backend for a web-front application, so the user interface is limited. A simple placeholder UI is in `main/GUI/splitGUI.m`, which is just a wrapper for a series of simulation "steps" `GUIsetup`, `GUIphystrck`, ..., `GUIsolver`. To use, just `cd` into a new project directory and launch `splitGUI()`, or continue working on an existing project with `splitGUI('my_project.ssp')`.

A batch version, designed for unattended operation (e.g. in a cluster), is in `main/GUI/splitscript.m`. By default, `splitscript('my_project.ssp')` will load an existing project (or create it new) and attempt all steps in order, but will crash as soon as a required file is not found or more than one option exists (see note on extensions below).

The folder `./tests/common` contains example files to run simple test cases:

- `*_simoptions.json`: Global options configuration file, see `main/options/DefaultOptions.m` for documentation.
- `*.meteo`: delimited text file with meteorological data (and meta-data). See `/matlab-pvmeteo/@MeteoData/getMeteoData.m` for specs. and alternatives.
- `*.sensors`: sensor configuration YAML file, used to parse the `*.meteo` file, and to define sensor meta-data (location, orientation and uncertainty).
- `*.samlib`: a module model definition file (plain text file) compatible with the old CEC database format, but with additional fields for the PVSyst model, and for reverse-bias characteristics (see `models/module/checkODM.m')`. A spreadsheet tool to test and generate these models can be found in `/models/module/ODMSpec_LibreOffice.ods`.
- `*.ondcp` (or `*.ond`, `*.samlib`): an inverter model definition file (text file). See `main/GUI/GUImodels.m` and `models/inverter/ONDread.m` for details.
- `*.tiff` (and/or `*.xyz`,`*.hor`): DEM raster files (and/or a point cloud text file, and/or a far-horizon-profile text file) that define terrain and surrounding obstacle geometry (see `main/GUI/GUIterrain.m`, and `terrain/readterrain.m`).

To run a test case, `cd` to this directory and run `testbase` to generate additional required files (plant geometry and electrical configuration):

```
cd(fullfile(BASEDIR,'tests'));
testbase('types','0a')
```

Running `testbase('0a',..)` will call `layout/samplesystem.m` which creates rectangular arrays with simple, repetitive geometry and electrical layout (`0a` stands for zero-axis-trackers, i.e. a fixed-mount plant). This will create a folder `./tests/0a` with additional files required for the simulation:
 
- `*.mpoly`: plain text file that defines internal module geometry (use `layout/stdmodulepoly.m` for common configurations).
- `*.tpoly`: plain text file that defines "mount (table or tracker) geometry and specs (see `layout/importpolygonsfile.m` and `main/GUI/GUIphystrck.m`)
- `*.mounts`: plain text file that defines geometrical layout, i.e. mount (axis) coordinates and orientations. See `layout/checkcoordsystem.m` and `layout/mountrotations.m` for notes on coordinate systems and tracker specifications.
- `*.arrdef`: plain text file that defines the electrical topology of the plant. Contains a bijective list of module's physical and electrical coordinates, i.e. which module at which mount is connected to which string at which MPPT, inverter, and combiner box. See `main/GUI/GUIarrdef.m` and `connectivity/readarraydefinition.m`.

Output files from a full simulation include the following:

- `*.ssp`: Solar-Simulation-Project files are just `*.mat` files with a different extension. You can see their contents directly with `load('myproject.ssp','-mat')`.
- `*.modint`: module model interpolant (binary `*.mat` file) is a module model that has been parsed into a set of IV curves that can be used for fast interpolation (see `models/module/ModuleInterpolant.m` for details).
- `horprof_*.mat`: near-horizon profiles, result of a simplified view-shed calculation at each mount's center point (see `terrain/gethorizons.m`).
- `*.shdres`: shading results (binary `*.mat` file). Contains an an object of the class `shading/ShadingResults.m` that encodes direct- and diffuse shading calculation results for a given plant geometry. These are independent of electrical layout, and, if performed over a grid of solar positions, also independent of time-step resolution (see `main/GUI/GUIshading.m` and `shading/ShadingAnalysis.m`).
- `solres_*.mat`: electrical simulation results (binary `*.mat` file). Contains a structure with production results at string level and for every simulation time-step, including information like voltage and clipping losses. See `main/GUI/GUIsolver.m` and `pvArraySolver.m`.

#### 6. Caveats and limitations

- Plant layout is currently limited to one mount type, and modules and inverters of the same model. As a work-around, define different inverters as separate projects, and add the rest of the mounts as static obstacles using an `*.xyz` point cloud file.

- Point-clouds are currently parsed through a Delaunay triangulation algorithm and merged into a single TIN together with DEM raster information. Unlike neighboring mounts, obstacles defined in this way do not project sharp shadows at the sub-module level (they are used only for near-horizon profile generation).

- There is still no support for bifacial modules, and ground-shading is still experimental.

- Shading analysis is slow. Migration of `ShadingAnalysis` using a graphics engine should be considered, before spending additional effort into the issues described above.

- There is still no account of ohmic wiring losses or power conditioning (trafo and MV transmission line) losses.

- The `MeteoData` and `MeteoSensor` classes are still not implemented beyond importing and simple quality-control. The Recursive Bayesian Estimation algorithm (`RBE`) is still a stand-alone library. The end goal is that sensor-shading information, and data-assimilated anisotropic state vector information are propagated into the irradiance transposition step (`poairradiance` and `diffusecomponents`).

- The ported and optimized classes `OneDiodeModel`, `mexODMpwlapprox`, and `mexVecODMpwlapproxAddseries` are still not fully merged into the electrical solver. The end goal is to port the complete `pvArraySolver` to avoid the current mex-interface overhead.

A (more) complete list of known issues and limitations is included in `todo.xlsx`. Volunteer contributions are very welcome!



