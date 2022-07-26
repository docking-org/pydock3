
c distance cube parameters are defined here
       integer cb(3)
c cb: integer x,y,z cube dimensions
       real cmin(3),cpa
c cmin: minimum edge x,y,z
c cpa: cubes per anstrom
       common /cubic/ cb,cmin,cpa

       pointer (i_cblower, cblower)
       pointer (i_cbupper, cbupper)
       pointer (i_cbatom, cbatom)
       common /cpntr/ i_cblower,i_cbupper,i_cbatom
