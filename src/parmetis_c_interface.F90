module parmetis_c_interface

  use iso_c_binding, only: c_int32_t, c_int64_t, c_float, c_double
  use mpi_f08, only: MPI_Comm

  use metis_c_interface, only: idx_t, real_t

  implicit none
  private

  integer(idx_t), parameter, public :: PARMETIS_MAJOR_VERSION = 4
  integer(idx_t), parameter, public :: PARMETIS_MINOR_VERSION = 0
  integer(idx_t), parameter, public :: PARMETIS_SUBMINOR_VERSION = 3

  ! Operation type codes ! pmoptype_et
  integer(idx_t), parameter, public :: PARMETIS_OP_KMETIS = 0
  integer(idx_t), parameter, public :: PARMETIS_OP_GKMETIS = 1
  integer(idx_t), parameter, public :: PARMETIS_OP_GMETIS = 2
  integer(idx_t), parameter, public :: PARMETIS_OP_RMETIS = 3
  integer(idx_t), parameter, public :: PARMETIS_OP_AMETIS = 4
  integer(idx_t), parameter, public :: PARMETIS_OP_OMETIS = 5
  integer(idx_t), parameter, public :: PARMETIS_OP_M2DUAL = 6
  integer(idx_t), parameter, public :: PARMETIS_OP_MKMETIS = 7

  ! ***********************************************************************
  ! Various constants used for the different parameters
  ! ***********************************************************************
  ! Matching types
  !> Restrict matching to within processor vertices
  integer(idx_t), parameter, public :: PARMETIS_MTYPE_LOCAL  = 1
  !> Remote vertices can be matched
  integer(idx_t), parameter, public :: PARMETIS_MTYPE_GLOBAL = 2

  ! Separator refinement types
  !> Vertices are visted from highest to lowest gain
  integer(idx_t), parameter, public :: PARMETIS_SRTYPE_GREEDY = 1
  !> Separators are refined in a two-phase fashion using `PARMETIS_SRTYPE_GREEDY` for the 2nd phase
  integer(idx_t), parameter, public :: PARMETIS_SRTYPE_2PHASE = 2

  ! Coupling types for ParMETIS_V3_RefineKway & ParMETIS_V3_AdaptiveRepart
  !> # of partitions == # of processors
  integer(idx_t), parameter, public :: PARMETIS_PSR_COUPLED = 1
  !> # of partitions != # of processors
  integer(idx_t), parameter, public :: PARMETIS_PSR_UNCOUPLED = 2

  ! Debug levels (fields should be ORed)
  !> Perform timing analysis
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_TIME = 1
  !> Perform timing analysis
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_INFO = 2
  !> Show the coarsening progress
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_PROGRESS = 4
  !> Show info on communication during folding
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_REFINEINFO = 8
  !> Show info on matching
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_MATCHINFO = 16
  !> Show info on communication during folding
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_RMOVEINFO = 32
  !> Determines if remapping will take place
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_REMAP = 64
  !> Performs a 2-hop matching
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_TWOHOP= 128
  !> Drop edges during coarsening
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_DROPEDGES= 256
  !> Reduces #trials for various steps
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_FAST = 512
  !> Saves non-active graphs to disk
  integer(idx_t), parameter, public :: PARMETIS_DBGLVL_ONDISK= 1024

  ! TODO: check which arguments are optional
  interface
    !> This routine is used to compute a k-way partitioning of a graph on
    !> processors using the multilevel k-way multi-constraint partitioning
    !> algorithm.
    function ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, &
                                  adjwgt, wgtflag, numflag, ncon, nparts, &
                                  tpwgts, ubvec, options, edgecut, part, &
                                  comm) &
                                  bind(C, name="ParMETIS_V3_PartKway") result(stat)

      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the vertices of the graph are
      !> distributed among the processors. (See discussion in Section 4.2.1).
      !> Its contents are identical for every processor.
      integer(idx_t), intent(in) :: vtxdist(*)
      !> adjacency structure of the graph at each processor.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: xadj(*)
      !> adjacency structure of the graph at each processor.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjncy(*)
      !> store the weights of the vertices and edges.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: vwgt(*)
      !> store the weights of the vertices and edges.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjwgt(*)
      !> This is used to indicate if the graph is weighted.
      !> `wgtflag` can take one of four values:
      !>
      !> - `0`: No weights (`vwgt` and `adjwgt` are both `NULL`).
      !> - `1`: Weights on the edges only (`vwgt` is `NULL`).
      !> - `2`: Weights on the vertices only (`adjwgt` is `NULL`).
      !> - `3`: Weights on both the vertices and edges.
      integer(idx_t), intent(in) :: wgtflag
      !> This is used to indicate the numbering scheme that is used
      !> for the `vtxdist`,`xadj`,`adjncy`, and `part` arrays. `numflag` can
      !> take one of two values:
      !>
      !> - `0`: C-style numbering that starts from `0`.
      !> - `1`: Fortran-style numbering that starts from `1`.
      integer(idx_t), intent(in) :: numflag
      !> This is used to specify the number of weights that each
      !> vertex has. It is also the number of balance constraints that must
      !> be satisfied.
      integer(idx_t), intent(in) :: ncon
      !> This is used to specify the number of subdomains that are
      !> desired. Note that the number of subdomains is independent of the
      !> number of processors that call this routine.
      integer(idx_t), intent(in) :: nparts
      !> An array of size `ncon×nparts` that is used to specify the
      !> fraction of vertex weight that should be distributed to each
      !> subdomain for each balance constraint. If all of the subdomains are
      !> to be of the same size for every vertex weight, then each of the
      !> `ncon×nparts` elements should be set to a value of `1/nparts`. If
      !> `ncon` is greater than `1`, the target subdomain weights for each
      !> subdomain are stored contiguously (similar to the `vwgt` array).
      !> Note that the sum of all of the `tpwgts` for a give vertex weight
      !> should be one.
      real(real_t), intent(in) :: tpwgts(*) ! ncon*nparts
      !> An array of size `ncon` that is used to specify the
      !> imbalance tolerance for each vertex weight, with `1` being perfect
      !> balance and `nparts` being perfect imbalance. A value of `1.05` for
      !> each of the `ncon` weights is recommended.
      real(real_t), intent(in) :: ubvec(ncon)
      !> This is an array of integers that is used to pass
      !> additional parameters for the routine. The first element
      !> (i.e.,`options[0]`) can take either the value of `0` or `1`. If it
      !> is `0`, then the default values are used, otherwise the remaining
      !> two elements of`options`are interpreted as follows:
      !>
      !> - `options[1]`: This specifies the level of information to be
      !>   returned during the execution of the algorithm. Timing information
      !>   can be obtained by setting this to `1`. Additional options for
      !>   this parameter can be obtained by looking at `parmetis.h`. The
      !>   numerical values there should be added to obtain the correct
      !>   value. The default value is `0`.
      !> - `options[2]`: This is the random number seed for the routine.
      integer(idx_t), intent(in) :: options(*) ! of len 3, could be optional
      !>  Upon successful completion, the number of edges that are
      !> cut by the partitioning is written to this parameter.
      integer(idx_t), intent(out) :: edgecut
      !> This is an array of size equal to the number of
      !> locally-stored vertices. Upon successful completion the partition
      !> vector of the locally-stored vertices is written to this array. (See
      !> discussion in Section 4.2.4).
      integer(idx_t), intent(out) :: part(*)
      !> This is a pointer to the MPI communicator of the processes
      !> that call `PARMETIS`. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK`: Indicates that the function returned normally.
      !> - `METIS_ERROR`: Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V3_PartKway

    !> This routine is used to compute a k-way partitioning of a graph on
    !> processors by combining the coordinate- based and multi-constraint
    !> k-way partitioning schemes.
    !>
    !> ### Note
    !>
    !> The quality of the partitionings computed by
    !> `ParMETIS_V3_PartGeomKway` are comparable to those produced by
    !> `ParMETIS_V3_PartKway`. However, the run time of the routine may be up
    !> to twice as fast.
    function ParMETIS_V3_PartGeomKway(vtxdist, xadj, adjncy, vwgt, &
                                      adjwgt, wgtflag, numflag, ndims, xyz, &
                                      ncon, nparts, tpwgts, ubvec, options, &
                                      edgecut, part, comm) &
                                      bind(C, name="ParMETIS_V3_PartGeomKway") result(stat)
      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the vertices of the graph are
      !> distributed among the processors. (See discussion in Section 4.2.1).
      !> Its contents are identical for every processor.
      integer(idx_t), intent(in) :: vtxdist(*)
      !> store the (local) adjacency structure of the
      !> graph at each processor. (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: xadj(*)
      !> store the (local) adjacency structure of the
      !> graph at each processor. (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjncy(*)
      !> store the weights of the vertices and edges.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: vwgt(*)
      !> store the weights of the vertices and edges.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjwgt(*)
      !> This is used to indicate if the graph is weighted.
      !> `wgtflag` can take one of four values:
      !>
      !> - `0`: No weights (`vwgt` and `adjwgt` are both `NULL`).
      !> - `1`: Weights on the edges only (`vwgt` is `NULL`).
      !> - `2`: Weights on the vertices only (`adjwgt` is `NULL`).
      !> - `3`: Weights on both the vertices and edges.
      integer(idx_t), intent(in) :: wgtflag
      !> This is used to indicate the numbering scheme that is used
      !> for the `vtxdist`, `xadj`, `adjncy`, and `part` arrays. `numflag`
      !> can take one of two values:
      !>
      !> - `0`: C-style numbering that starts from `0`.
      !> - `1`: Fortran-style numbering that starts from `1`.
      integer(idx_t), intent(in) :: numflag
      !> The number of dimensions of the space in which the graph is
      !> embedded.
      integer(idx_t), intent(in) :: ndims
      !> The array storing the coordinates of the vertices (described
      !> in Section 4.2.2).
      real(real_t), intent(in) :: xyz(*)
      !> This is used to specify the number of weights that each
      !> vertex has. It is also the number of balance constraints that must
      !> be satisfied.
      integer(idx_t), intent(in) :: ncon
      !> This is used to specify the number of subdomains that are
      !> desired. Note that the number of subdomains is independent of the
      !> number of processors that call this routine.
      integer(idx_t), intent(in) :: nparts
      !> An array of size `ncon×nparts` that is used to specify the
      !> fraction of vertex weight that should be distributed to each
      !> subdomain for each balance constraint. If all of the subdomains are
      !> to be of the same size for every vertex weight, then each of
      !> the`ncon×nparts` elements should be set to a value of `1/nparts`. If
      !> `ncon` is greater than one, the target subdomain weights for each
      !> subdomain are stored contiguously (similar to the `vwgt` array).
      !> Note that the sum of all of the `tpwgts` for a give vertex weight
      !> should be one.
      real(real_t), intent(in) :: tpwgts(*) ! ncon*nparts
      !> An array of size `ncon` that is used to specify the
      !> imbalance tolerance for each vertex weight, with `1` being perfect
      !> balance and `nparts` being perfect imbalance. A value of `1.05` for
      !> each of the `ncon` weights is recommended.
      real(real_t), intent(in) :: ubvec(ncon)
      !> This is an array of integers that is used to pass
      !> parameters to the routine. Their meanings are identical to those of
      !> `ParMETIS_V3_PartKway`.
      integer(idx_t), intent(in) :: options(*) ! of len 3, could be optional
      !> Upon successful completion, the number of edges that are
      !> cut by the partitioning is written to this parameter.
      integer(idx_t), intent(out) :: edgecut
      !> This is an array of size equal to the number of
      !> locally-stored vertices. Upon successful completion the partition
      !> vector of the locally-stored vertices is written to this array. (See
      !> discussion in Section 4.2.4).
      integer(idx_t), intent(out) :: part(*)
      !> This is a pointer to the MPI communicator of the processes
      !> that call `PARMETIS`. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK`: Indicates that the function returned normally.
      !> - `METIS_ERROR`: Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V3_PartGeomKway

    !> This routine is used to compute ap-way partitioning of a graph on
    !> processors using a coordinate-based space- filling curves method.
    !>
    !> ### Note
    !>
    !> The quality of the partitionings computed by `ParMETIS_V3_PartGeom`
    !> are significantly worse than those produced by `ParMETIS_V3_PartKway`
    !> and `ParMETIS_V3_PartGeomKway`.
    function ParMETIS_V3_PartGeom(vtxdist, ndims, xyz, part, comm) &
                                  bind(C, name="ParMETIS_V3_PartGeom") result(stat)
      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the vertices of the graph are
      !> distributed among the processors. (See discussion in Section 4.2.1).
      !> Its contents are identical for every processor.
      integer(idx_t), intent(in) :: vtxdist(*)
      !> The number of dimensions of the space in which the graph is
      !> embedded.
      integer(idx_t), intent(in) :: ndims
      !> The array storing the coordinates of the vertices (described
      !> in Section 4.2.2).
      real(real_t), intent(in) :: xyz(*)
      !> This is an array of size equal to the number of locally
      !> stored vertices. Upon successful completion stores the partition
      !> vector of the locally stored graph (described in Section 4.2.4).
      integer(idx_t), intent(out) :: part(*)
      !> This is a pointer to the MPI communicator of the processes
      !> that call PARMETIS. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK` Indicates that the function returned normally.
      !> - `METIS_ERROR` Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V3_PartGeom

    !> This routine is used to compute a k-way partitioning of a mesh on
    !> processors. The mesh can contain elements of different types.
    function ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt, wgtflag, &
                                      numflag, ncon, ncommonnodes, nparts, &
                                      tpwgts, ubvec, options, edgecut, part, comm) &
                                      bind(C, name="ParMETIS_V3_PartMeshKway") result(stat)
      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the elements of the mesh are
      !> distributed among the processors. It is analogous to the `vtxdist`
      !> array. Its contents are identical for every processor. (See
      !> discussion in Section 4.2.3).
      integer(idx_t), intent(in) :: elmdist(*)
      !> specifies the elements that are stored
      !> locally at each processor. (See discussion in Section 4.2.3).
      integer(idx_t), intent(in) :: eptr(*)
      !> specifies the elements that are stored
      !> locally at each processor. (See discussion in Section 4.2.3).
      integer(idx_t), intent(in) :: eind(*)
      !> This array stores the weights of the elements. (See
      !> discussion in Section 4.2.3).
      integer(idx_t), intent(in) :: elmwgt(*)
      !> This is used to indicate if the elements of the mesh have
      !> weights associated with them. The `wgtflag` can take two values:
      !> - 0 No weights (`elmwgt` is `NULL`).
      !> - 2 Weights on the vertices only.
      integer(idx_t), intent(in) :: wgtflag
      !> This is used to indicate the numbering scheme that is used
      !> for the `elmdist`, `elements`, and `part` arrays. `numflag` can take
      !> one of two values:
      !> - `0`: C-style numbering that starts from `0`.
      !> - `1`: Fortran-style numbering that starts from `1`.
      integer(idx_t), intent(in) :: numflag
      !> This is used to specify the number of weights that each
      !> vertex has. It is also the number of balance constraints that must
      !> be satisfied.
      integer(idx_t), intent(in) :: ncon
      !> This parameter determines the degree of connectivity
      !> among the vertices in the dual graph. Specifically, an edge is
      !> placed between any two elements if and only if they share at least
      !> this many nodes. This value should be greater than `0`, and for most
      !> meshes a value of two will create reasonable dual graphs. However,
      !> depending on the type of elements in the mesh, values greater than
      !> `2` may also be valid choices. For example, for meshes containing
      !> only triangular, tetrahedral, hexahedral, or rectangular elements,
      !> this parameter can be set to two, three, four, or two, respectively.
      !> Note that setting this parameter to a small value will increase the
      !> number of edges in the resulting dual graph and the corresponding
      !> partitioning time.
      integer(idx_t), intent(in) :: ncommonnodes
      !> This is used to specify the number of subdomains that are
      !> desired. Note that the number of subdomains is independent of the
      !> number of processors that call this routine.
      integer(idx_t), intent(in) :: nparts
      !> An array of size `ncon×nparts` that is used to specify the
      !> fraction of vertex weight that should be distributed to each
      !> subdomain for each balance constraint. If all of the subdomains are
      !> to be of the same size for every vertex weight, then each of the
      !> `ncon×nparts` elements should be set to a value of `1/nparts`. If
      !> `ncon` is greater than 1, the target subdomain weights for each
      !> subdomain are stored contiguously (similar to the `vwgt` array).
      !> Note that the sum of all of the `tpwgts` for a give vertex weight
      !> should be one.
      real(real_t), intent(in) :: tpwgts(*) ! ncon*nparts
      !> An array of size `ncon` that is used to specify the
      !> imbalance tolerance for each vertex weight, with 1 being perfect
      !> balance and `nparts` being perfect imbalance. A value of 1.05 for
      !> each of the `ncon` weights is recommended.
      real(real_t), intent(in) :: ubvec(ncon)
      !> This is an array of integers that is used to pass
      !> parameters to the routine. Their meanings are identical to those of
      !> `ParMETIS_V3_PartKway`.
      integer(idx_t), intent(in) :: options(*) ! of len 3, could be optional
      !> Upon successful completion, the number of edges that are
      !> cut by the partitioning is written to this parameter.
      integer(idx_t), intent(out) :: edgecut
      !> This is an array of size equal to the number of
      !> locally-stored vertices. Upon successful completion the partition
      !> vector of the locally-stored vertices is written to this array. (See
      !> discussion in Section 4.2.4).
      integer(idx_t), intent(out) :: part(*)
      !> This is a pointer to the MPI communicator of the processes
      !> that call `PARMETIS`. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK` Indicates that the function returned normally.
      !> - `METIS_ERROR` Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V3_PartMeshKway

    !> This routine is used to balance the work load of a graph that
    !> corresponds to an adaptively refined mesh.
    function ParMETIS_V3_AdaptiveRepart(vtxdist, xadj, adjncy, vwgt, &
                                        vsize, adjwgt, wgtflag, numflag, ncon, &
                                        nparts, tpwgts, ubvec, ipc2redist, &
                                        options, edgecut, part, comm) &
                                        bind(C, name="ParMETIS_V3_AdaptiveRepart") result(stat)
      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the vertices of the graph are
      !> distributed among the processors. (See discussion in Section 4.2.1).
      !> Its contents are identical for every processor.
      integer(idx_t), intent(in) :: vtxdist(*)
      !> store the (local) adjacency structure of the
      !> graph at each processor. (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: xadj(*)
      !> store the (local) adjacency structure of the
      !> graph at each processor. (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjncy(*)
      !> store the weights of the vertices and edges.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: vwgt(*)
      !> This array stores the size of the vertices with respect to
      !> redistribution costs. Hence, vertices associated with mesh elements
      !> that require a lot of memory will have larger corresponding entries
      !> in this array. Otherwise, this array is similar to the `vwgt` array.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: vsize(*)
      !> store the weights of the vertices and edges.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjwgt(*)
      !> This is used to indicate if the graph is weighted.
      !> `wgtflag` can take one of four values:
      !>
      !> - `0`: No weights (`vwgt` and adjwgt are both `NULL`).
      !> - `1`: Weights on the edges only (`vwgt` is `NULL`).
      !> - `2`: Weights on the vertices only (`adjwgt` is `NULL`).
      !> - `3`: Weights on both the vertices and edges.
      integer(idx_t), intent(in) :: wgtflag
      !> This is used to indicate the numbering scheme that is used
      !> for the `vtxdist`, `xadj`, `adjncy`, and `part` arrays. `numflag`
      !> can take one of two values:
      !>
      !> - `0`: C-style numbering that starts from `0`.
      !> - `1`: Fortran-style numbering that starts from `1`.
      integer(idx_t), intent(in) :: numflag
      !> This is used to specify the number of weights that each
      !> vertex has. It is also the number of balance constraints that must
      !> be satisfied.
      integer(idx_t), intent(in) :: ncon
      !> This is used to specify the number of subdomains that are
      !> desired. Note that the number of subdomains is independent of the
      !> number of processors that call this routine.
      integer(idx_t), intent(in) :: nparts
      !> An array of size `ncon×nparts` that is used to specify the
      !> fraction of vertex weight that should be distributed to each
      !> subdomain for each balance constraint. If all of the subdomains are
      !> to be of the same size for every vertex weight, then each of the
      !> `ncon×nparts` elements should be set to a value of `1/nparts`. If
      !> `ncon` is greater than one, the target subdomain weights for each
      !> subdomain are stored contiguously (similar to the`vwgt`array). Note
      !> that the sum of all of the `tpwgts` for a give vertex weight should
      !> be one.
      real(real_t), intent(in) :: tpwgts(*) ! ncon*nparts
      !> An array of size `ncon` that is used to specify the
      !> imbalance tolerance for each vertex weight, with 1 being perfect
      !> balance and `nparts` being perfect imbalance. A value of `1.05` for
      !> each of the `ncon` weights is recommended.
      real(real_t), intent(in) :: ubvec(ncon)
      !> This parameter describes the ratio of inter-processor
      !> communication time compared to data redistribution time. It should
      !> be set between `0.000001` and `1000000.0`. If `ipc2redist` is set high, a
      !> repartitioning with a low edge-cut will be computed. If it is set
      !> low, a repartitioning that requires little data redistribution will
      !> be computed. Good values for this parameter can be obtained by
      !> dividing inter-processor communication time by data redistribution
      !> time. Otherwise, a value of `1000.0` is recommended.
      real(real_t), intent(in) :: ipc2redist(*)
      !> This is an array of integers that is used to pass
      !> additional parameters for the routine. The first element
      !> (i.e.,`options[0]`) can take either the value of `0` or `1`. If it
      !> is `0`, then the default values are used, otherwise the remaining
      !> three elements of `options` are interpreted as follows:
      !>
      !> - `options[1]` This specifies the level of information to be
      !>   returned during the execution of the algorithm. Timing information
      !>   can be obtained by setting this to 1. Additional options for this
      !>   parameter can be obtained by looking at `parmetis.h`. The
      !>   numerical values there should be added to obtain the correct
      !>   value. The default value is `0`.
      !> - `options[2]` This is the random number seed for the routine.
      !> - `options[3]` This specifies whether the subdomains and processors
      !>   are coupled or un-coupled. If the number of subdomains desired
      !>   (i.e.,`nparts`) and the number of processors that are being used
      !>   is not the same, then these must be un-coupled. However, if
      !>   `nparts` equals the number of processors, these can either be
      !>   coupled or decoupled. If subdomains and processors are coupled,
      !>   then the initial partitioning will be obtained implicitly from the
      !>   graph distribution. However, if subdomains are un-coupled from
      !>   processors, then the initial partitioning needs to be obtained
      !>   from the initial values assigned to the `part` array. A value of
      !>   `PARMETIS_PSR_COUPLED` indicates that subdomains and processors
      !>   are coupled and a value of `PARMETIS_PSR_UNCOUPLED` indicates that
      !>   these are de-coupled. The default value is `PARMETIS_PSR_COUPLED`
      !>   if `nparts` equals the number of processors and
      !>   `PARMETIS_PSR_UNCOUPLED`(un-coupled) otherwise. These constants
      !>   are defined in `parmetis.h`.
      integer(idx_t), intent(in) :: options(*) ! of len 3, could be optional
      !> Upon successful completion, the number of edges that are
      !> cut by the partitioning is written to this parameter.
      integer(idx_t), intent(out) :: edgecut
      !> This is an array of size equal to the number of
      !> locally-stored vertices. Upon successful completion the partition
      !> vector of the locally-stored vertices is written to this array. (See
      !> discussion in Section 4.2.4). If the number of processors is not
      !> equal to the number of subdomains and/or `options[3]` is set to
      !> `PARMETIS_PSR_UNCOUPLED`, then the previously computed partitioning
      !> must be passed to the routine as a parameter via this array.
      integer(idx_t), intent(out) :: part(*)
      !> This is a pointer to the MPI communicator of the processes
      !> that call `PARMETIS`. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK` Indicates that the function returned normally.
      !> - `METIS_ERROR` Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V3_AdaptiveRepart

    !> This routine is used to improve the quality of an existing a k-way
    !> partitioning on processors using the multilevel k-way refinement
    !> algorithm.
    function ParMETIS_V3_RefineKway(vtxdist, xadj, adjncy, vwgt, &
                                    adjwgt, wgtflag, numflag, ncon, nparts, &
                                    tpwgts, ubvec, options, edgecut, &
                                    part, comm) &
                                    bind(C, name="ParMETIS_V3_RefineKway") result(stat)
      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the vertices of the graph are
      !> distributed among the processors. (See discussion in Section 4.2.1).
      !> Its contents are identical for every processor.
      integer(idx_t), intent(in) :: vtxdist(*)
      !> Store the (local) adjacency structure of the
      !> graph at each processor. (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: xadj(*)
      !> Store the (local) adjacency structure of the
      !> graph at each processor. (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjncy(*)
      !> Store the weights of the vertices and edges.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: vwgt(*)
      !> Store the weights of the vertices and edges.
      !> (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjwgt(*)
      !> This is used to indicate if the graph is weighted.
      !> `wgtflag` can take one of four values:
      !>
      !> - `0`: No weights (`vwgt` and `adjwgt` are both `NULL`).
      !> - `1`: Weights on the edges only (`vwgt` is `NULL`).
      !> - `2`: Weights on the vertices only (`adjwgt` is `NULL`).
      !> - `3`: Weights on both the vertices and edges.
      integer(idx_t), intent(in) :: wgtflag
      !> This is used to indicate the numbering scheme that is used
      !> for the `vtxdist`, `xadj`, `adjncy`, and `part` arrays. `numflag`
      !> can take the following two values:
      !>
      !> - `0`: C-style numbering is assumed that starts from `0`
      !> - `1`: Fortran-style numbering is assumed that starts from `1`
      integer(idx_t), intent(in) :: numflag
      !> This is used to specify the number of weights that each
      !> vertex has. It is also the number of balance constraints that must
      !> be satisfied.
      integer(idx_t), intent(in) :: ncon
      !> This is used to specify the number of subdomains that are
      !> desired. Note that the number of subdomains is independent of the
      !> number of processors that call this routine.
      integer(idx_t), intent(in) :: nparts
      !> An array of size `ncon×nparts` that is used to specify the
      !> fraction of vertex weight that should be distributed to each
      !> subdomain for each balance constraint. If all of the subdomains are
      !> to be of the same size for every vertex weight, then each of the
      !> `ncon×nparts` elements should be set to a value of `1/nparts`. If
      !> `ncon` is greater than `1`, the target subdomain weights for each
      !> subdomain are stored contiguously (similar to the `vwgt` array).
      !> Note that the sum of all of the `tpwgts` for a give vertex weight
      !> should be one.
      real(real_t), intent(in) :: tpwgts(*) ! ncon*nparts
      !> An array of size `ncon` that is used to specify the
      !> imbalance tolerance for each vertex weight, with 1 being perfect
      !> balance and `nparts` being perfect imbalance. A value of 1.05 for
      !> each of the `ncon` weights is recommended.
      real(real_t), intent(in) :: ubvec(ncon)
      !> This is an array of integers that is used to pass
      !> parameters to the routine. Their meanings are identical to those of
      !> `ParMETIS_V3_AdaptiveRepart`.
      integer(idx_t), intent(in) :: options(*) ! of len 3, could be optional
      !> Upon successful completion, the number of edges that are
      !> cut by the partitioning is written to this parameter.
      integer(idx_t), intent(out) :: edgecut
      !> This is an array of size equal to the number of
      !> locally-stored vertices. Upon successful completion the partition
      !> vector of the locally-stored vertices is written to this array. (See
      !> discussion in Section 4.2.4). If the number of processors is not
      !> equal to the number of subdomains and/or `options[3]` is set to
      !> `PARMETIS_PSR_UNCOUPLED`, then the previously computed partitioning
      !> must be passed to the routine as a parameter via this array.
      integer(idx_t), intent(out) :: part(*)
      !> This is a pointer to the MPI communicator of the processes
      !> that call `PARMETIS`. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK` Indicates that the function returned normally.
      !> - `METIS_ERROR` Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V3_RefineKway

    !> This routine is used to compute a fill-reducing ordering of a sparse
    !> matrix using multilevel nested dissection.
    function ParMETIS_V3_NodeND(vtxdist, xadj, adjncy, &
                                numflag, options, order, sizes, comm) &
                                bind(C, name="ParMETIS_V3_NodeND") result(stat)
      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the vertices of the graph are
      !> distributed among the processors. (See discussion in Section 4.2.1).
      !> Its contents are identical for every processor.
      integer(idx_t), intent(in) :: vtxdist(*)
      !> Store the (local) adjacency structure of the
      !> graph at each processor (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: xadj(*)
      !> Store the (local) adjacency structure of the
      !> graph at each processor (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjncy(*)
      !> This is used to indicate the numbering scheme that is used
      !> for the `vtxdist`, `xadj`, `adjncy`, and `order` arrays. `numflag`
      !> can take the following two values:
      !> - `0`: C-style numbering is assumed that starts from `0`
      !> - `1`: Fortran-style numbering is assumed that starts from `1`
      integer(idx_t), intent(in) :: numflag
      !> This is an array of integers that is used to pass
      !>  parameters to the routine. Their meanings are identical to those
      !>  of `ParMETIS_V3_PartKway`.
      integer(idx_t), intent(in) :: options(*) ! of len 3, could be optional
      !> This array returns the result of the ordering (described in
      !> Section 4.2.4).
      integer(idx_t), intent(out) :: order(*)
      !> This array returns the number of nodes for each subdomain
      !> and each separator (described in Section 4.2.4).
      integer(idx_t), intent(out) :: sizes(*)
      !> This is a pointer to the MPI communicator of the processes
      !> that call `PARMETIS`. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK`: Indicates that the function returned normally.
      !> - `METIS_ERROR`: Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V3_NodeND

    !> This routine is used to compute a fill-reducing ordering of a sparse
    !> matrix using multilevel nested dissection.
    function ParMETIS_V32_NodeND(vtxdist, xadj, adjncy, vwgt, &
                                  numflag, mtype, rtype, p_nseps, s_nseps, &
                                  ubfrac, seed, dbglvl, order, &
                                  sizes, comm) &
                                  bind(C, name="ParMETIS_V32_NodeND") result(stat)
      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the vertices of the graph are
      !> distributed among the processors. (See discussion in Section 4.2.1).
      !> Its contents are identical for every processor.
      integer(idx_t), intent(in) :: vtxdist(*)
      !> Store the (local) adjacency structure of the
      !> graph at each processor (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: xadj(*)
      !> Store the (local) adjacency structure of the
      !> graph at each processor (See discussion in Section 4.2.1).
      integer(idx_t), intent(in) :: adjncy(*)
      !> These store the weights of the vertices. A value of `NULL`
      !> indicates that each vertex has unit weight. (See discussion in
      !> Section 4.2.1).
      integer(idx_t), intent(in) :: vwgt(*)
      !> This is used to indicate the numbering scheme that is used
      !> for the`vtxdist`, `xadj`, `adjncy`, and `order` arrays. The possible
      !> values are:
      !>
      !> - `0`: C-style numbering is assumed that starts from 0
      !> - `1`: Fortran-style numbering is assumed that starts from 1
      integer(idx_t), intent(in) :: numflag
      !> This is used to indicate the scheme to be used for computing
      !> the matching. The possible values, defined in `parmetis.h`, are:
      !>
      !> - `PARMETIS_MTYPE_LOCAL` A local matching scheme is used in which
      !>   each pair of matched vertices reside on the same processor.
      !> - `PARMETIS_MTYPE_GLOBAL` A global matching scheme is used in which
      !>   the pairs of matched vertices can reside on different processors.
      !>   This is the default value if a `NULL` value is passed.
      integer(idx_t), intent(in) :: mtype
      !> This is used to indicate the separator refinement scheme
      !> that will be used. The possible values, defined in `parmetis.h`,
      !> are:
      !>
      !> - `PARMETIS_SRTYPE_GREEDY`: Uses a simple greedy refinement
      !>   algorithm.
      !> - `PARMETIS_SRTYPE_2PHASE`: Uses a higher quality refinement
      !>   algorithm, which is somewhat slower. This is the default value if
      !>   a `NULL` value is passed.
      integer(idx_t), intent(in) :: rtype
      !> Specifies the number of different separators that will be
      !> computed during each bisection at the first
      integer(idx_t), intent(in) :: p_nseps
      !> Specifies the number of different separators that will be
      !> computed during each of the bisections levels of the remaining
      !> levels of the nested dissection (when the matrix has been divided
      !> among the processors and each processor proceeds independently to
      !> order its portion of the matrix). The bisections that achieve the
      !> smallest separator are selected. The default value is `1` (when
      !> `NULL` is supplied), but values greater than `1` can lead to better
      !> quality orderings. However, this is a time- quality trade-off.
      integer(idx_t), intent(in) :: s_nseps
      !> This value indicates how unbalanced the two partitions are
      !> allowed to get during each bisection level. The default value (when
      !> `NULL` is supplied) is `1.05`, but higher values (typical ranges
      !> `1.05–1.25`) can lead to smaller separators.
      real(real_t), intent(in) :: ubfrac
      !> This is the seed for the random number generator. When `NULL`
      !> is supplied, a default seed is used.
      integer(idx_t), intent(in) :: seed
      !> This specifies the level of information to be returned
      !> during the execution of the algorithm. This is identical to the
      !> `options[2]` parameter of the other routines. When `NULL` is
      !> supplied, a value of `0` is used.
      integer(idx_t), intent(in) :: dbglvl
      !> This array returns the result of the ordering (described in
      !> Section 4.2.4).
      integer(idx_t), intent(out) :: order(*)
      !> This array returns the number of nodes for each subdomain
      !> and each separator (described in Section 4.2.4).
      integer(idx_t), intent(out) :: sizes(*)
      !> This is a pointer to the MPI communicator of the processes
      !> that call `PARMETIS`. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK`: Indicates that the function returned normally.
      !> - `METIS_ERROR`: Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V32_NodeND

    !> This routine is used to construct a distributed graph given a
    !> distributed mesh. It can be used in conjunction with other routines in
    !> the `PARMETIS` library. The mesh can contain elements of different
    !> types.
    function ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, numflag, &
                                  ncommonnodes, xadj, adjncy, comm) &
                                  bind(C, name="ParMETIS_V3_Mesh2Dual") result(stat)
      import :: idx_t, real_t, MPI_Comm
      !> This array describes how the elements of the mesh are
      !> distributed among the processors. It is analogous to the `vtxdist`
      !> array. Its contents are identical for every processor. (See
      !> discussion in Section 4.2.3).
      integer(idx_t), intent(in) :: elmdist(*)
      !> Specifies the elements that are stored
      !> locally at each processor. (See discussion in Section 4.2.3).
      integer(idx_t), intent(in) :: eptr(*)
      !> Specifies the elements that are stored
      !> locally at each processor. (See discussion in Section 4.2.3).
      integer(idx_t), intent(in) :: eind(*)
      !> This is used to indicate the numbering scheme that is used
      !> for the `elmdist`, `elements`, `xadj`, and `adjncy` arrays.
      !> `numflag` can take one of two values:
      !>
      !> - `0`: C-style numbering that starts from `0`.
      !> - `1`: Fortran-style numbering that starts from `1`.
      integer(idx_t), intent(in) :: numflag
      !> This parameter determines the degree of connectivity
      !> among the vertices in the dual graph. Specifically, an edge is
      !> placed between any two elements if and only if they share at least
      !> this many nodes. This value should be greater than 0, and for most
      !> meshes a value of two will create reasonable dual graphs. However,
      !> depending on the type of elements in the mesh, values greater than 2
      !> may also be valid choices. For example, for meshes containing only
      !> triangular, tetrahedral, hexahedral, or rectangular elements, this
      !> parameter can be set to two, three, four, or two, respectively. Note
      !> that setting this parameter to a small value will increase the
      !> number of edges in the resulting dual graph and the corresponding
      !> partitioning time.
      integer(idx_t), intent(in) :: ncommonnodes
      !> Upon the successful completion of the routine,
      !> pointers to the constructed `xadj` and `adjncy` arrays will be
      !> written to these parameters. (See discussion in Section 4.2.1). The
      !> calling program is responsible for freeing this memory by calling
      !> the `METISFree()` routine described in `METIS` manual.
      integer(idx_t), intent(out) :: xadj(*)    ! TODO: is pointer of pointers
      !> Upon the successful completion of the routine,
      !> pointers to the constructed `xadj` and `adjncy` arrays will be
      !> written to these parameters. (See discussion in Section 4.2.1). The
      !> calling program is responsible for freeing this memory by calling
      !> the `METISFree()` routine described in `METIS` manual.
      integer(idx_t), intent(out) :: adjncy(*)  ! TODO: is pointer of pointers
      !> This is a pointer to the MPI communicator of the processes
      !> that call `PARMETIS`. For most programs this will point to
      !> `MPI_COMM_WORLD`.
      type(MPI_Comm), intent(in) :: comm
      !> - `METIS_OK`: Indicates that the function returned normally.
      !> - `METIS_ERROR`: Indicates some other type of error.
      integer(idx_t) :: stat
    end function ParMETIS_V3_Mesh2Dual
  end interface

end module parmetis_c_interface
