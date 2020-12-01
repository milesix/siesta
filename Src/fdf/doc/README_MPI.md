For MPI operation, the "reading node" has to call the fdf initialization routine,
serialize the fdf data structure, and broadcast it. Other nodes receive the
structure and de-serialize it. In the following snippet we assume that 0 is the
reading node:

```
      if (Node .eq. 0) then
         call fdf_init(filein, fileout)
      else
         !                                                                                                                                   
      endif
#ifdef MPI      
      call broadcast_fdf_struct(0,mpi_comm_world)
      !                                                                                                                                      
      ! Example of multi-node logging (optional)                                                                                                         
      !                                                                                                                                      
      if ( Node .ne. 0 .and.
     $     fdf_get('fdf-log-in-all-nodes', .false.)) then

         debug_level = fdf_get('fdf-debug', 0)
         write(prefix,"(a,i0)") "fdf-node-", Node
         call fdf_setdebug(debug_level,trim(prefix)//".debug")
         output_level = fdf_get('fdf-output', 1)
         call fdf_setoutput(output_level,trim(prefix)//".log")
         if (debug_level >= 2) call fdf_print_struct()
      endif
#endif
```

The routine `broadcast_fdf_struct` can be found in file `broadcast_fdf_struct.F90`
can be found in this directory. It should be general enough for most codes.

