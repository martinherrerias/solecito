# PV4K-matlab
Simulate PV plants with so much detail that you can see their pores.

Capdevila, H., Herrerias A., M., Marola, A., 2014. Anisotropic Diffuse Shading Model for Sun-tracking Photovoltaic Systems. Energy Procedia, 2013 ISES Solar World Congress 57, 144–151. https://doi.org/10.1016/j.egypro.2014.10.018
Capdevila, H., Marola, A., Herrerias A., M., 2013. High resolution shading modeling and performance simulation of sun-tracking photovoltaic systems. 9th International Conference on Concentrator Photovoltaic Syste,s: CPV-9, Miyazaki, Japan, pp. 201–204. https://doi.org/10.1063/1.4822231
Herrerias A., M., Capdevila, H., 2018. Detailed Calculation of Electrical Mismatch Losses for Central- and String-Inverter Configurations on Utility-Scale PV Arrays, in: 35th European Photovoltaic Solar Energy Conference and Exhibition. pp. 1973–1978. https://doi.org/10.4229/35thEUPVSEC20182018-6CV.2.36
Herrerias A., M., Capdevila, H., Zhou, H., Hammer, A., 2020. Simulation of Large PV Plants Using a Continuous Radiance Distribution Model and Cell-Resolution Mismatch Calculation. EU PVSEC 2020 37th European Photovoltaic Solar Energy Conference and Exhibition, online, pp. 1311–1316.
Zhou, H., Niethammer, C., Herrerias A., M., 2021. Usage Experiences of Performance Tools for Modern C++Code Analysis and Optimization, in: Tools for High Performance Computing 2018 / 2019. Springer International Publishing, Cham, pp. 103–121. https://doi.org/10.1007/978-3-030-66057-4_5


## Use
#### 1. Clone the repository (*with submodules!*):

`git clone --recurse-submodules git@github.com:martinherrerias/PV4K-matlab.git`

#### 2. Add code to the matlab path (note some folders are removed, to avoid path cluttering and name conflicts):

```
BASEDIR = '~/PV4K-matlab';  % your cloned repository
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

#### 5. Test



## Acknowledgements

Since 2019, the code has been reworked in the framework of the HyForPV project, funded by the Federal Ministry of Economic Affairs and Energy (Germany) on the basis of a decision by the German
Bundestag.

This project was comissioned and supported privately for years by Hugo Capdevila. I trust he would be happy to see it bear fruit, even if not in the way he envisioned.

Thanks to the financial support of the Mexican National Council of Science and Technology (CONACYT) and the German Academic Exchange Service (DAAD), in 2011-2013.

