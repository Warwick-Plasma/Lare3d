MODULE output

  USE shared_data
  USE iocommon

  IMPLICIT NONE

  SAVE

  PRIVATE

  PUBLIC :: cfd_open_clobber, cfd_write_block_header, cfd_write_meshtype_header
  PUBLIC :: cfd_safe_write_string, cfd_write_snapshot_data, &
      cfd_write_stitched_vector
  PUBLIC :: cfd_write_stitched_magnitude, cfd_write_real_constant
  PUBLIC :: cfd_write_visit_expression

CONTAINS

  SUBROUTINE cfd_open_clobber(filename)

    CHARACTER(len = *), INTENT(IN) :: filename

    ! Set the block header
    block_header_size = max_string_len * 2 + 4 + 2 * 8

    ! Delete file and wait
    IF (cfd_rank == default_rank) &
        CALL MPI_FILE_DELETE(TRIM(filename), MPI_INFO_NULL, cfd_errcode)

    CALL MPI_BARRIER(cfd_comm, cfd_errcode)
    CALL MPI_FILE_OPEN(cfd_comm, TRIM(filename), cfd_mode, MPI_INFO_NULL, &
        cfd_filehandle, cfd_errcode)

    IF (cfd_rank == default_rank) THEN
      ! Write the header
      CALL MPI_FILE_WRITE(cfd_filehandle, "CFD", 3, MPI_CHARACTER, cfd_status, &
          cfd_errcode)

      ! This goes next so that stuff can be added to the global header
      ! without breaking
      ! Everything
      CALL MPI_FILE_WRITE(cfd_filehandle, header_offset, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, block_header_size, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, cfd_version, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, cfd_revision, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, max_string_len, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)

      ! This is where the nblocks variable will go, put a zero for now
      CALL MPI_FILE_WRITE(cfd_filehandle, 0, 1, MPI_INTEGER, cfd_status, &
          cfd_errcode)
    END IF

    ! Currently no blocks written
    nblocks = 0
    ! Current displacement is just the header
    current_displacement = header_offset

  END SUBROUTINE cfd_open_clobber



  SUBROUTINE cfd_safe_write_string(string)

    CHARACTER(LEN = *), INTENT(IN) :: string
    CHARACTER(LEN = max_string_len) :: output
    INTEGER :: len_s

    len_s = LEN(string)

    ! This subroutine expects that the record marker is in place and that
    ! the view is set correctly. Call it only on the node which is doing the
    ! writing
    ! You still have to advance the file pointer yourself on all nodes

    output(1:MIN(max_string_len, len_s)) = string(1:MIN(max_string_len, len_s))

    ! If this isn't the full string length then tag in a ACHAR(0) to help
    ! With C + + string handling
    IF (len_s + 1 < max_string_len) output(len_s + 1:max_string_len) = ACHAR(0)

    CALL MPI_FILE_WRITE(cfd_filehandle, output, max_string_len, MPI_CHARACTER, &
        cfd_status, cfd_errcode)

  END SUBROUTINE cfd_safe_write_string



  SUBROUTINE cfd_write_block_header(blockname, blockclass, blocktype, &
      blocklength, blockmetadatalength, rank_write)

    CHARACTER(len = *), INTENT(IN) :: blockname, blockclass
    INTEGER, INTENT(IN) :: blocktype, rank_write
    INTEGER(KIND = 8), INTENT(IN) :: blocklength, blockmetadatalength
    INTEGER :: len_bn, len_bc

    len_bn = LEN(blockname)
    len_bc = LEN(blockclass)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL cfd_safe_write_string(blockname)
      CALL cfd_safe_write_string(blockclass)
    END IF
    current_displacement = current_displacement + 2 * max_string_len

    ! Write the block type
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, blocktype, 1, MPI_INTEGER, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + 4

    ! Write the block skip and metadata skip data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER8, &
        MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, blockmetadatalength, 1, &
          MPI_INTEGER8, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, blocklength, 1, MPI_INTEGER8, &
          cfd_status, cfd_errcode)
    END IF

    current_displacement = current_displacement + 2 * 8

    nblocks = nblocks + 1

  END SUBROUTINE cfd_write_block_header



  SUBROUTINE cfd_write_meshtype_header(type, dim, sof, rank_write)

    ! MeshTypes (Meshes, fluid variables, multimat blocks etc)
    ! All have a common header, this is what writes that (although the content
    ! Of type will depend on what meshtype you're using)

    INTEGER, INTENT(IN) :: type, dim, rank_write, sof

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, type, 1, MPI_INTEGER, cfd_status, &
          cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, dim, 1, MPI_INTEGER, cfd_status, &
          cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, sof, 1, MPI_INTEGER, cfd_status, &
          cfd_errcode)
    END IF

    current_displacement = current_displacement + meshtype_header_offset

  END SUBROUTINE cfd_write_meshtype_header



  SUBROUTINE cfd_write_snapshot_data(time, CYCLE, rank_write)

    INTEGER, INTENT(IN) :: rank_write, CYCLE
    INTEGER(8) :: mdlength
    REAL(8), INTENT(IN) :: time

    mdlength = soi + num

    CALL cfd_write_block_header("Snapshot", "Snapshot", TYPE_SNAPSHOT, &
        mdlength, mdlength, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, CYCLE, 1, MPI_INTEGER, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, "native", MPI_INFO_NULL, &
        cfd_errcode)

    IF (cfd_rank == rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, time, 1, MPI_DOUBLE_PRECISION, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + 8

  END SUBROUTINE cfd_write_snapshot_data



  SUBROUTINE cfd_write_stitched_vector(vector_name, vector_class, mesh_name, &
      mesh_class, name, class, rank_write)

    CHARACTER(len = *), DIMENSION(:), INTENT(IN) :: name, class
    CHARACTER(len = *), INTENT(IN) :: vector_name, vector_class, mesh_name, &
        mesh_class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: n_dims, mdlength, blocklength
    INTEGER :: i_loop

    n_dims = SIZE(name)

    mdlength = 2 * max_string_len + soi
    blocklength = mdlength + n_dims * 2 * max_string_len

    CALL cfd_write_block_header(vector_name, vector_class, &
        TYPE_STITCHED_VECTOR, blocklength, mdlength, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    END IF
    current_displacement = current_displacement + 2 * max_string_len

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, n_dims, 1, MPI_INTEGER, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      DO i_loop = 1, n_dims
        CALL cfd_safe_write_string(name(i_loop))
        CALL cfd_safe_write_string(class(i_loop))
      END DO
    END IF
    current_displacement = current_displacement + 2 * n_dims* max_string_len

  END SUBROUTINE cfd_write_stitched_vector



  SUBROUTINE cfd_write_stitched_magnitude(magn_name, magn_class, mesh_name, &
      mesh_class, name, class, rank_write)

    CHARACTER(len = *), DIMENSION(:), INTENT(IN) :: name, class
    CHARACTER(len = *), INTENT(IN) :: magn_name, magn_class, mesh_name, &
        mesh_class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: n_dims, mdlength, blocklength
    INTEGER :: i_loop

    n_dims = SIZE(name)

    mdlength = 2 * max_string_len + soi
    blocklength = mdlength + n_dims * 2 * max_string_len

    CALL cfd_write_block_header(magn_name, magn_class, &
        TYPE_STITCHED_MAGNITUDE, blocklength, mdlength, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    END IF
    current_displacement = current_displacement + 2 * max_string_len

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, n_dims, 1, MPI_INTEGER, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      DO i_loop = 1, n_dims
        CALL cfd_safe_write_string(name(i_loop))
        CALL cfd_safe_write_string(class(i_loop))
      END DO
    END IF
    current_displacement = current_displacement + 2 * n_dims* max_string_len

  END SUBROUTINE cfd_write_stitched_magnitude



  SUBROUTINE cfd_write_real_constant(name, class, value, rank_write)

    CHARACTER(len = *), INTENT(IN) :: name, class
    REAL(num), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: mdlength

    mdlength = num

    CALL cfd_write_block_header(name, class, TYPE_CONSTANT, mdlength, &
        mdlength, rank_write)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, value, 1, mpireal, cfd_status, &
          cfd_errcode)
    END IF
    current_displacement = current_displacement + num

  END SUBROUTINE cfd_write_real_constant



  SUBROUTINE cfd_write_1d_integer_array(name, class, values, rank_write)

    CHARACTER(len = *), INTENT(IN) :: name, class
    INTEGER, DIMENSION(:), INTENT(IN) :: values
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: mdlength

    mdlength = 2 * soi

    CALL cfd_write_block_header(name, class, TYPE_INTEGERARRAY, mdlength, &
        mdlength, rank_write)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      ! 1D
      CALL MPI_FILE_WRITE(cfd_filehandle, 1, 1, MPI_INTEGER, cfd_status, &
          cfd_errcode)
      ! Size of array
      CALL MPI_FILE_WRITE(cfd_filehandle, 1, SIZE(values), MPI_INTEGER, &
          cfd_status, cfd_errcode)
      ! Actual Array
      CALL MPI_FILE_WRITE(cfd_filehandle, values, SIZE(values), MPI_INTEGER, &
          cfd_status, cfd_errcode)
    END IF
    current_displacement = current_displacement + mdlength

  END SUBROUTINE cfd_write_1d_integer_array



  SUBROUTINE cfd_write_visit_expression(expression_name, expression_class, &
      expression)

    CHARACTER(LEN = *), DIMENSION(:), INTENT(IN) :: expression_name, &
        expression_class, expression

    PRINT *, LEN(expression(1)), LEN(expression(2))

  END SUBROUTINE cfd_write_visit_expression

END MODULE output
