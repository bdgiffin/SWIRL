#include "gtest/gtest.h"
#include <math.h>
#include "SpatialHash.h"

TEST(SpatialHash_tests, find_nearby_segments) {
  // Define spatial hash object
  SpatialHash hash;

  // Set the dimensions of the spatial hash
  double xmin = 0.0;
  double ymin = 0.0;
  double zmin = 0.0;
  double xmax = 1.0;
  double ymax = 1.0;
  double zmax = 1.0;
  int      Nx = 10;
  int      Ny = 10;
  int      Nz = 10;
  hash.set_dimensions(xmin,ymin,zmin,
		      xmax,ymax,zmax,
		        Nx,  Ny,  Nz);

  // Insert new segment spanning the diagonal of the spatial hash
  int segment_id = 1;
  double xyz1[3] = { 0.0, 0.0, 0.0 };
  double xyz2[3] = { 1.0, 1.0, 1.0 };
  hash.insert_segment(segment_id,xyz1[0],xyz1[1],xyz1[2],xyz2[0],xyz2[1],xyz2[2]);

  // Find the range of grid cells overlapping with a bounding box enclosing a spherical point
  double xpoint[3] = { 1.0, 1.0, 0.0 };
  double radius = std::sqrt(3.0)/2.0 + 1.0e-5;
  int i1, j1, k1, i2, j2, k2;
  bool found = hash.find_grid_cells_overlapping_bounding_box(xpoint[0]-radius,xpoint[1]-radius,xpoint[2]-radius,
	 					             xpoint[0]+radius,xpoint[1]+radius,xpoint[2]+radius,
						             i1,j1,k1,i2,j2,k2);

  // Find nearby segments
  int* nearby_segment_ids = nullptr;
  int Nfound = hash.find_nearby_segments((xmax-xmin)/Nx,(ymax-ymin)/Ny,(zmax-zmin)/Nz,nearby_segment_ids);
  //ASSERT_EQ(Nfound,1);
  for (int i=0; i<Nfound; i++) {
    ASSERT_EQ(nearby_segment_ids[i],segment_id);
  }
  
} /* TEST(SpatialHash_tests, find_nearby_segments) */
