# Water molecules clustering algorithm

Algorithm was designed to cluster water particles from MD simulations based on their coordinates into equally sized groups. It is used to aggregate non-bounded MD (water) molecules in order to map their parameters into the coarse-grained model (such as based on dissipative particle dynamics). See the publication below for a full description of the procedure:

Pieczywek, P.M., Płaziński, W. & Zdunek, A. Dissipative particle dynamics model of homogalacturonan based on molecular dynamics simulations. Sci Rep 10, 14691 (2020). https://doi.org/10.1038/s41598-020-71820-2


The input data should be a matrix containing the coordinates of the molecules for the successive MD simulation frames. Sample data is included in the file "test_data.mat". Script “example.mat” clarifies how the functions works.

Cite As:
Pieczywek, P.M., Płaziński, W. & Zdunek, A. Dissipative particle dynamics model of homogalacturonan based on molecular dynamics simulations. Sci Rep 10, 14691 (2020). https://doi.org/10.1038/s41598-020-71820-2