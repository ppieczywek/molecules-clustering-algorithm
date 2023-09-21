
load('my_test_data.mat')
%  Sample data comes from MD simulation of water molecules in a small
%  simulation box (edge length = 31.036 A, cubic shape, PBC assumed). 
%  Only the heavy atoms are considered during the analysis, therefore the
%  data contain only the positions of the oxygen atoms. 
%  The assumed cg_ratio is equal to 4.  

cg_ratio = 4;
aggregation_tolerance = 2;
grid_cell_size = 7.0;%
x_limits = [0.0 31.036];
y_limits = [0.0 31.036];
z_limits = [0.0 31.036];

bead_position = aggregate(md_trajectories, cg_ratio, aggregation_tolerance, x_limits, y_limits, z_limits, grid_cell_size );

hold on
bead_id = 10;
plot3(squeeze(bead_position(bead_id,1,1:20)), ...
      squeeze(bead_position(bead_id,2,1:20)), ...
      squeeze(bead_position(bead_id,3,1:20)), ...
      '-r');
plot3(squeeze(bead_position(bead_id,1,1:20)), ...
      squeeze(bead_position(bead_id,2,1:20)), ...
      squeeze(bead_position(bead_id,3,1:20)), ...
      'ob');

axis('equal');  
