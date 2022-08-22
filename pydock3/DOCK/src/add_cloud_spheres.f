      subroutine add_cloud_spheres(cluster, rigid_heavy_coords, 
     &                    rigid_colors, rigid_heavy_count, db2lig)

      use db2type
      use spheres

      implicit none

      integer cluster !which cluster we're adding this time
      real rigid_heavy_coords(3, MAXPTS) !coordinates of heavy atoms in rigid
      integer rigid_colors(MAXPTS) !actually the only used ligand colors
      integer rigid_heavy_count  !how many heavy atoms, used in matching
      type(db2) :: db2lig

      integer sphcount !temporary counter
      integer xyzcount !temporary counter for coordinates

      do sphcount = db2lig%cluster_match_start(cluster), 
     &                            db2lig%cluster_match_end(cluster)
        rigid_heavy_count = rigid_heavy_count + 1
        rigid_colors(rigid_heavy_count) = 
     &      db2lig%addmatch_color(sphcount)
        do xyzcount = 1, 3
          rigid_heavy_coords(xyzcount, rigid_heavy_count) = 
     &        db2lig%addmatch_coord(xyzcount, sphcount)
        enddo
      enddo

      return
      end

