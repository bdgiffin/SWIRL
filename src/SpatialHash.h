#ifndef SPATIAL_HASH_H
#define SPATIAL_HASH_H

#include <vector>
#include <algorithm>
#include <iostream>

class SpatialHash {
public:

  // --------------- Declare public member functions --------------- //
  
  // Default constructor
  SpatialHash(void) { }

  // Set the dimensions of the spatial hash
  void set_dimensions(double xmin, double ymin, double zmin,
		      double xmax, double ymax, double zmax,
		      int Nx_in, int Ny_in, int Nz_in) {
    // set the grid's base point
    x0 = xmin;
    y0 = ymin;
    z0 = zmin;

    // set the number of grid cells per axis
    Nx = Nx_in;
    Ny = Ny_in;
    Nz = Nz_in;

    // set the grid cell dimensions
    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;
    dz = (zmax - zmin)/Nz;

    // allocate and initialize the grid
    hashed_segment_ids.clear();
    hashed_segment_ids.resize(Nx*Ny*Nz);
    
  } // set_dimensions()

  // Insert a new line segment into the hash
  void insert_segment(int segment_id, double x1, double y1, double z1, double x2, double y2, double z2) {
    // find the starting and ending normalized grid coordinates
    double xg1 = (x1-x0)/dx;
    double yg1 = (y1-y0)/dy;
    double zg1 = (z1-z0)/dz;
    double xg2 = (x2-x0)/dx;
    double yg2 = (y2-y0)/dy;
    double zg2 = (z2-z0)/dz;
    double xg12 = xg2 - xg1;
    double yg12 = yg2 - yg1;
    double zg12 = zg2 - zg1;

    // determine the starting and ending grid indicies
    int i1 = xg1;
    int j1 = yg1;
    int k1 = zg1;
    int i2 = xg2;
    int j2 = yg2;
    int k2 = zg2;

    // determine the oriented number of grid planes crossed in each spatial direction
    int ni = i2 - i1;
    int nj = j2 - j1;
    int nk = k2 - k1;
    int si = (ni < 0) ? -1 : +1;
    int sj = (nj < 0) ? -1 : +1;
    int sk = (nk < 0) ? -1 : +1;
    ni *= si;
    nj *= sj;
    nk *= sk;

    // keep track of the grid indicies crossed by the current segment
    std::vector<int> grid_indicies;

    // get the index of the first grid cell
    grid_indicies.push_back(Nx*(Ny*k1 + j1) + i1);

    // for each spatial direction, loop over intersected grid planes
    for (int i = 0; i < ni; i++) {
      // determine the x coordinate of the current intersected yz grid plane
      double xi = i1 + si*i + 0.5*(1+si);
      
      // determine the parametric coordinate along the line segment at the intersection point
      double s = (xi - xg1)/xg12;

      // determine the y,z coordinates at the intersection point
      double yi = yg1 + s*yg12;
      double zi = zg1 + s*zg12;

      // determine the hashed index of the newly entered grid cell
      int index = Nx*(Ny*int(zi) + int(yi)) + int(xi+0.5*si);

      // keep track of the currently indexed grid cell
      grid_indicies.push_back(index);
    }
    for (int j = 0; j < nj; j++) {
      // determine the y coordinate of the current intersected xz grid plane
      double yj = j1 + sj*j + 0.5*(1+sj);
      
      // determine the parametric coordinate along the line segment at the intersection point
      double s = (yj - yg1)/yg12;

      // determine the x,z coordinates at the intersection point
      double xj = xg1 + s*xg12;
      double zj = zg1 + s*zg12;

      // determine the hashed index of the newly entered grid cell
      int index = Nx*(Ny*int(zj) + int(yj+0.5*sj)) + int(xj);

      // keep track of the currently indexed grid cell
      grid_indicies.push_back(index);
    }
    for (int k = 0; k < nk; k++) {
      // determine the z coordinate of the current intersected xy grid plane
      double zk = k1 + sk*k + 0.5*(1+sk);
      
      // determine the parametric coordinate along the line segment at the intersection point
      double s = (zk - zg1)/zg12;

      // determine the x,y coordinates at the intersection point
      double xk = xg1 + s*xg12;
      double yk = yg1 + s*yg12;

      // determine the hashed index of the newly entered grid cell
      int index = Nx*(Ny*int(zk+0.5*sk) + int(yk)) + int(xk);

      // keep track of the currently indexed grid cell
      grid_indicies.push_back(index);
    }

    // loop over all intersected grid indicies
    for (int i = 0; i < grid_indicies.size(); i++) {
      // insert the current segment ID into the current intersected grid index of the spatial hash
      hashed_segment_ids[grid_indicies[i]].push_back(segment_id);
    }
    
  } // insert_segment()

  // Find the list of all candidate segments in proximity to the specified point
  int find_nearby_segments(double x, double y, double z, int*& nearby_segment_ids) {
    // determine the hashed index of the current grid cell
    int index = hash_index(x,y,z);

    // return the list of segments associated with the indexed grid cell
    if (index < 0) {
      nearby_segment_ids = nullptr;
      return 0;
    } else {
      std::vector<int>& indexed_cell = hashed_segment_ids[index];
      nearby_segment_ids = indexed_cell.data();
      return indexed_cell.size();
    }
  } // find_nearby_segments()

  // Find the list of all segments belonging to the indicated grid index
  int find_segments_in_grid_cell(int i, int j, int k, int*& nearby_segment_ids) {
    // determine the hashed index of the requested grid cell
    int index = Nx*(Ny*k + j) + i;

    // return the list of segments associated with the indexed grid cell
    std::vector<int>& indexed_cell = hashed_segment_ids[index];
    nearby_segment_ids = indexed_cell.data();
    return indexed_cell.size();
  } // find_segments_in_grid_cell()

  // Find the range of grid cells that overlap with the indicated bounding box domain
  bool find_grid_cells_overlapping_bounding_box(double x1, double y1, double z1, double x2, double y2, double z2,
						int& i1, int& j1, int& k1, int& i2, int& j2, int& k2) {
    // determine the spatial grid indicies associated with the corners of the bounding box
    i1 = std::max(int((x1-x0)/dx),0);
    j1 = std::max(int((y1-y0)/dy),0);
    k1 = std::max(int((z1-z0)/dz),0);
    i2 = std::min(int((x2-x0)/dx),Nx-1);
    j2 = std::min(int((y2-y0)/dy),Ny-1);
    k2 = std::min(int((z2-z0)/dz),Nz-1);

    // return true if a non-null range of grid cells was found
    return (((i2-i1) >= 0) && ((j2-j1) >= 0) && ((k2-k1) >= 0));
    
  } // find_grid_cells_overlapping_bounding_box()

private:

  // --------------- Declare private member functions --------------- //

  // determine the hashed spatial index for a given point
  int hash_index(double x, double y, double z) {
    int i = (x-x0)/dx;
    int j = (y-y0)/dy;
    int k = (z-z0)/dz;
    if ((i < 0) || (j < 0) || (k < 0) || (i >= Nx) || (j >= Ny) || (k >= Nz)) {
      return -1; // this point does not lie within the bounds of the spatial grid
    } else {
      return Nx*(Ny*k + j) + i;
    }
  } // hash_index

  // --------------- Declare private data members --------------- //

  // Spatial hash size and dimensions
  double x0, y0, z0;
  double dx, dy, dz;
  int Nx, Ny, Nz;

  // Spatial hash of segment IDs
  std::vector<std::vector<int> > hashed_segment_ids;

}; // SpatialHash

#endif /* SPATIAL_HASH_H */
