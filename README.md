# EigenOPF
Copyright (C) 2020 Chao Duan

Generation and Demand Co-optimization for Small-signal Stability Enhancement.


This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


The full text of the GNU General Public License can be found in the file "LICENSE.txt".


# Dependencies

MATLAB 2017b or later.

Matpower 5.1 https://matpower.org/download/

IPOPT Solver or Knitro Solver, a free version of IPOPT can be found in https://www.inverseproblem.co.nz/OPTI/index.php/DL/DownloadOPTI

# Usage
The main function TX_opt maximizes the system damping ratio by simultaneously adjusting generation and demand.

[success,optratio,improve,Vps1] = TX_opt(ps,Vps0,dispratio,load_level)

The function takes as inputs a data structure specifying the transmission system (ps), initial power flow solution (Vps0), relative increase/decrease allowed for the loads (dispratio), and load level (load_level). It returns the indicator for successful optimization (success), the optimized damping ratio (optratio), percentages improvement of the damping ratio (improve), and the power flow solution after optimization (Vps1).

See call_OptimizeRatio.m (and OptimizeRatio.m) for an example of using TX_opt.

The detailed power system model and techiques to derive the linearized system associated derivative information are presented in [appendices.](https://github.com/cduan2020/EigenOPF/blob/master/appendices.pdf)



