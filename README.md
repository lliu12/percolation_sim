# Explosive rigidity percolation in kirigami

<img src = "https://github.com/lliu12/percolation_sim/blob/main/explosive_kirigami.png" width="750" height="450" />

The code contains several methods for rigidifying a kirigami system by incrementally changing either the connectivity or the rigidity of individual components, which allow us to control the explosive rigidity percolation in kirigami.

Any comments and suggestions are welcome. 

If you use this code in your own work, please cite the following paper:

G. P. T. Choi, L. Liu, and L. Mahadevan, "[Explosive rigidity percolation in kirigami.](https://arxiv.org/abs/2211.15073)" Preprint, arXiv:2211.15073, 2022.

============================================================

Main programs:
* `kirigami_explosive_percolation_connection.m`: The connection-based method.
* `kirigami_explosive_percolation_angle.m`: The angle-based method.
* `kirigami_explosive_percolation_coordinate.m` and `simulation.py`: The coordinate-based method.

Results:
* `result_connection`: The connection-based method.
* `result_angle`: The angle-based method.
* `result_coordinates`: The coordinate-based method.
