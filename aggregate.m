%    Algorithm was designed to cluster particles based on their coordinates
%    in space into equally sized groups. It is used to aggregate non-bounded
%    MD (water) molecules in order to map their parameters into the
%    coarse-grained model. See the publication below for a full description
%    of the procedure:
%    Pieczywek, P.M., P³aziñski, W. & Zdunek, A. Dissipative particle
%    dynamics model of homogalacturonan based on molecular dynamics
%    simulations. Sci Rep 10, 14691 (2020).
%    https://doi.org/10.1038/s41598-020-71820-2
%
%    Copyright (C) 2021  Piotr Mariusz Pieczywek

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
function [ bead_position ] = aggregate(md_trajectories, cg_ratio, aggregation_tolerance, x_limits, y_limits, z_limits, grid_cell_size )
    
%   INPUT PARAMETERS 
%
%   md_trajectories        - [N by 3 by F] matrix, which holds x,y,z
%                            coordinates of MD atoms from MD simulation.
%                            N is the number of atoms, F is the number of data
%                            frames from the MD simulation. In case of water 
%                            particles data usually contains only coordinates
%                            of oxygen atoms.  
%   cg_ratio               - pisitive integer which defines how many atoms
%                            will be aggregated into one CG bead/particle
%   aggregation_tolerance  - positive integer, number of CG beads, that
%                            should be left over; the actual number of 
%                            aggregates is equal to:
%                            (atoms_num / cg_ratio) - aggregation_tolerance; 
%                            leaving a few atoms unassigned to any cluster 
%                            improves the stability of the algorithm
%   x_limits, 
%   y_limits, 
%   z_limits               - limits of MD simulation box in x,y,z directions,
%                            written as 2-elemen vectors [min max] this depends
%                            on yours data
%   grid_cell_size         - edge length of rectangular grid used by the
%                            nearest neighbour search algorithm; In general, 
%                            the bigger the unit cell is, the easier to find
%                            the neigbouring atoms. However, large unit cells
%                            significantly slows down calculations. 
%
%    Final output data is written to "bead_position" matrix, which holds
%    CG bead coordinates for each processed simulation frame.

    h = waitbar(0.0,'calculations in progress');
    bead_position = [];
    
    if isempty(md_trajectories)
        disp('md_trajectories is empty!')
        return;
    end
    
    [atoms_number atomic_coordinates max_frames] = size(md_trajectories);
    
    if atomic_coordinates ~= 3 
        disp('md_trajectories expected to cointain x,y,z atomic coordinates!')
        return;
    end
    
    if ~ismatrix(md_trajectories(:,:,1))
        disp('md_trajectories is expected to be a matrix!')
    	return;
    end
        
    if cg_ratio <= 1
        disp('cg_ratio is expected to be a positive integer larger than 1!')
        return;
    else
        if floor(cg_ratio) ~= cg_ratio
            disp('cg_ratio is expected to be a positive integer!')
            return;
        end
    end
    
    if aggregation_tolerance < 0
        disp('aggregation_tolerance is expected to be a non-negative integer!')
        return;
    else
        if floor(aggregation_tolerance) ~= aggregation_tolerance
            disp('aggregation_tolerance is expected to be a integer!')
            return;
        end
    end
    
    if isvector(x_limits)
        if length(x_limits) ~= 2
            disp('x_limits is expected to be a 2-element vactor!')
            return;
        else
            if x_limits(1) > x_limits(2)
                disp('lower limit of x_limits is expected to be smaller than upper limit - x_limits(1) << x_limits(2) !')
                return;
            end
        end
    else
        disp('x_limits is expected to be a 2-element vactor!')
        return;
    end
    
    if isvector(y_limits)
        if length(y_limits) ~= 2
            disp('y_limits is expected to be a 2-element vactor!')
            return;
        else
            if y_limits(1) > y_limits(2)
                disp('lower limit of y_limits is expected to be smaller than upper limit - y_limits(1) << y_limits(2) !')
                return;
            end
        end
    else
        disp('y_limits is expected to be a 2-element vactor!')
        return;
    end
    
    if isvector(z_limits)
        if length(z_limits) ~= 2
            disp('z_limits is expected to be a 2-element vactor!')
            return;
        else
            if z_limits(1) > z_limits(2)
                disp('lower limit of z_limits is expected to be smaller than upper limit - z_limits(1) << xz_limits(2) !')
                return;
            end
        end
    else
        disp('z_limits is expected to be a 2-element vactor!')
        return;
    end    
    
    if grid_cell_size < 0.0
        disp('grid_cell_size is expected to be a positive real number!')
        return;
    end
    
    start_frame = 1;
    bead_num = floor(atoms_number / cg_ratio)  - aggregation_tolerance;

    % define the grid for the nearest neighbour search algorithm
    x_bins = floor(x_limits(2) / grid_cell_size);
    y_bins = floor(y_limits(2) / grid_cell_size);
    z_bins = floor(z_limits(2) / grid_cell_size);

    x_span = (x_limits(2) - x_limits(1)) / x_bins; 
    y_span = (y_limits(2) - y_limits(1)) / y_bins; 
    z_span = (z_limits(2) - z_limits(1)) / z_bins; 


    % Initial positions of DPD beads are generated from randomly
    % picked atoms from MD model.
    start_idx = randperm(atoms_number);
    start_idx = start_idx(1:bead_num);

    % An output matrix with DPD bead positions
    bead_position = zeros(bead_num,3,max_frames);
    bead_position(:,:,1) = md_trajectories(start_idx,:,1);

    for frame_id=(start_frame+1):1:(max_frames) 

        md_atoms = md_trajectories(:,:,frame_id);

        nearest_atom   = zeros(bead_num,5) - 1;
        atom_has_bead  = zeros(size(md_atoms,1),1);
        bead_has_atom  = zeros(bead_num,1);
        atoms_assigned = zeros(bead_num,1);

        grid = cell(x_bins, y_bins, z_bins);

        for ii=1:1:size(md_atoms,1)

            x_cell = ceil(md_atoms(ii,1)/x_span);
            y_cell = ceil(md_atoms(ii,2)/y_span);
            z_cell = ceil(md_atoms(ii,3)/z_span);

            if x_cell > x_bins
                x_cell = x_bins;
            end

            if y_cell > y_bins
                y_cell = y_bins;
            end

            if z_cell > z_bins
                z_cell = z_bins;
            end
            grid{x_cell, y_cell, z_cell} = [grid{x_cell, y_cell, z_cell}; ii];
        end

        old_bead_position = bead_position(:,:,(frame_id-start_frame)); 
        new_bead_position = zeros(size(old_bead_position));

        for cg_level=1:1:cg_ratio

            for itetarion=1:1:20

                atom_nh = cell(size(md_atoms,1),1);

                % In the following loop script assigns DPD bead to the grid cell
                % and looks for the closest atom in this cell and the sourounding
                % cells. Algorithm assumes that periodic boundary conditions
                % are used. If the DPD bead is on the edge of the box, atoms will
                % also be searched on the opposite side of the simulation area.

                for ii=1:1:bead_num

                    if bead_has_atom(ii) == 0

                        x_cell = ceil(old_bead_position(ii,1)/x_span);
                        y_cell = ceil(old_bead_position(ii,2)/y_span);
                        z_cell = ceil(old_bead_position(ii,3)/z_span);

                        shift_x = 0.0;
                        shift_y = 0.0;
                        shift_z = 0.0;

                        xx = 0;
                        yy = 0;
                        zz = 0;

                        for pz = (z_cell-1):1:(z_cell+1)

                            shift_z = 0.0;
                            zz = pz;
                            if pz <= 0 
                                shift_z = -z_limits(2);
                                zz = z_bins;
                            end
                            if pz > z_bins
                                shift_z = z_limits(2);
                                zz = 1;
                            end

                            for py = (y_cell-1):1:(y_cell+1)

                                shift_y = 0.0;
                                yy = py;
                                if py <= 0 
                                    shift_y = -y_limits(2);
                                    yy = y_bins;
                                end
                                if py > y_bins
                                    shift_y = y_limits(2);
                                    yy = 1;
                                end

                                for px = (x_cell-1):1:(x_cell+1)

                                    shift_x = 0.0;
                                    xx = px;
                                    if px <= 0 
                                        shift_x = -x_limits(2);
                                        xx = x_bins;
                                    end
                                    if px > x_bins
                                        shift_x = x_limits(2);
                                        xx = 1;
                                    end

                                    if ~isempty(grid{xx,yy,zz})

                                        atoms_in_cell = grid{xx,yy,zz};
                                        for kk=1:1:length(atoms_in_cell)

                                            if atom_has_bead(atoms_in_cell(kk)) == 0

                                                atom_position = md_atoms(atoms_in_cell(kk),:);
                                                atom_position(1) = atom_position(1) + shift_x;
                                                atom_position(2) = atom_position(2) + shift_y;
                                                atom_position(3) = atom_position(3) + shift_z;

                                                distance = old_bead_position(ii,:) - atom_position;
                                                distance = sqrt(sum(distance.^2));

                                                if nearest_atom(ii,1) == -1
                                                    nearest_atom(ii,1) = atoms_in_cell(kk);
                                                    nearest_atom(ii,2) = distance;
                                                    nearest_atom(ii,3:5) = atom_position;    
                                                else    
                                                    if distance < nearest_atom(ii,2)
                                                        nearest_atom(ii,1) = atoms_in_cell(kk);
                                                        nearest_atom(ii,2) = distance;
                                                        nearest_atom(ii,3:5) = atom_position;
                                                    end 
                                                end

                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                % Following loop runs through all MD atoms and makes a list
                % of the DPD beads to which the atom has been assigned. Initially 
                % one atom can be defined as the closest one to a few DPD beads.
                % Following lines of code ensure that only one - the closest one -
                % atom will be assigned to one DPD molecule.
                for ii=1:1:bead_num
                    if nearest_atom(ii,1) > 0
                        atom_nh{nearest_atom(ii,1)} = [ atom_nh{nearest_atom(ii,1)}; [ii nearest_atom(ii,2:5)]];
                    end
                end

                atom_cnt = 0;
                for ii=1:1:size(md_atoms,1)

                    nh = atom_nh{ii};
                    if ~isempty(nh) 

                         % Sorts list of atoms based on distance from their
                         % beads, from the closest to the farthest.
                         nh = sortrows(nh,2);

                         bead_has_atom(nh(1,1)) = 1;
                         atom_cnt = atom_cnt + 1; 
                         atom_has_bead(ii) = 1;

                         if size(nh,1) > 1
                            for qq=2:1:size(nh,1)
                                nearest_atom(nh(qq,1), : ) = -1;
                            end
                         end
                    end
                end

                % The loop breaks and code proceedes to next cg level iteration 
                % if all beads have one atom assigned or the number of tries 
                % exceeds 20 itetarions.
                if atom_cnt == bead_num
                    break
                end

            end

            % This loop accumulates positions of assigned atoms,  
            % flags assigned atoms and increase the assigned atoms counter for
            % each bead.
            for ii=1:1:bead_num
                if nearest_atom(ii,1) > 0
                    new_bead_position(ii,:) =  new_bead_position(ii,:) + nearest_atom(ii,3:5);
                    atom_has_bead(nearest_atom(ii,1)) = 1;
                    atoms_assigned(ii) = atoms_assigned(ii) + 1;  
                end
            end
            bead_has_atom  = zeros(bead_num,1);
            nearest_atom   = zeros(bead_num,5) - 1;
        end

        % Calculates the averaged position of asigned atoms.
        new_bead_position =  new_bead_position ./ repmat(atoms_assigned,1,3);

        % This loop translates the position of the DPD bead across pbc edges if
        % necessary.
        for ii=1:1:bead_num
            if new_bead_position(ii,1) < 0
               new_bead_position(ii,1) =  x_limits(2) + new_bead_position(ii,1);
            end
            if new_bead_position(ii,2) < 0
               new_bead_position(ii,2) =  y_limits(2) + new_bead_position(ii,2);
            end
            if new_bead_position(ii,3) < 0
               new_bead_position(ii,3) =  z_limits(2) + new_bead_position(ii,3);
            end
            if new_bead_position(ii,1) > x_limits(2)
               new_bead_position(ii,1) = new_bead_position(ii,1) - x_limits(2);
            end
            if new_bead_position(ii,2) > y_limits(2)
               new_bead_position(ii,2) = new_bead_position(ii,2) - y_limits(2);
            end
            if new_bead_position(ii,3) > z_limits(2)
               new_bead_position(ii,3) = new_bead_position(ii,3) - z_limits(2);
            end
        end

        bead_position(:,:,(frame_id)) = new_bead_position;
        waitbar((frame_id-start_frame)/max_frames ,h,'calculations in progress');

    end

    delete(h)

end

