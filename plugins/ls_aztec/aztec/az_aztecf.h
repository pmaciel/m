












#define AZ_MSG_TYPE      1234
#define AZ_NUM_MSGS        20



#ifndef AZ_MAX_MEMORY_SIZE
#define AZ_MAX_MEMORY_SIZE   16000000  
                                       
#endif
#ifndef AZ_MAX_MSG_BUFF_SIZE
#define AZ_MAX_MSG_BUFF_SIZE 100000    
#endif
#define AZ_MAX_NEIGHBORS     250
#define AZ_MAX_MESSAGE_SIZE  (AZ_MAX_MSG_BUFF_SIZE / (2*AZ_MAX_NEIGHBORS))
#define AZ_FALSE               0
#define AZ_TRUE                1
#define AZ_MAX_POLY_ORDER     10 
#define AZ_default           -10 



#define AZ_OPTIONS_SIZE       15
#define AZ_PARAMS_SIZE         5
#define AZ_PROC_SIZE           3
#define AZ_STATUS_SIZE         7
#define AZ_COMM_SIZE          (10 + 3*AZ_MAX_NEIGHBORS)



#define AZ_cg               0 
#define AZ_gmres            1 
#define AZ_cgs              2 
#define AZ_tfqmr            3 
#define AZ_bicgstab         4 
#define AZ_slu              5 
#define AZ_symmlq           6 
#define AZ_lu               8 
                              
                              
                              





#define AZ_BJacobi          2 
#define AZ_row_sum          3 
#define AZ_sym_diag         4 
#define AZ_sym_row_sum      5 
#define AZ_sym_BJacobi      6 
                              
                              



#define AZ_none             0 
                              
#define AZ_Jacobi           1 
                              
#define AZ_sym_GS           2 
#define AZ_Neumann          3 
#define AZ_ls               4 
#define AZ_ilu              6 
#define AZ_bilu             7 

/*************************************************************/
/* my modification to include the block milu preconditioner. */
/*************************************************************/

#define AZ_bmilu            9
#define AZ_icc             10


/*************************************************************/
/*                 end of my modification                    */
/*************************************************************/



#define AZ_r0               0 
#define AZ_rhs              1 
#define AZ_Anorm            2 
#define AZ_sol              3 
#define AZ_weighted         4 
                              
                              



#define AZ_all             -3 
                              
                              

#define AZ_last            -1 
#define AZ_warnings        -2 



#define AZ_input_form       0 
                              
                              
                              
                              
#define AZ_global_mat       1 
                              
                              
                              
                              
                              
                              
                              
                              
#define AZ_explicit         2 
                              
                              
                              
                              
                              
                              



#define AZ_calc             1 
#define AZ_recalc           2 
#define AZ_reuse            3 
#define AZ_sys_reuse        4 
                              
                              




#define AZ_diag             1 
#define AZ_sym_full         2 
#define AZ_full             3 
                              
                              



#define AZ_classic          0
#define AZ_modified         1



#define AZ_resid            0
#define AZ_rand             1



#define AZ_normal           0 
#define AZ_param            1 
#define AZ_breakdown        2 
#define AZ_maxits           3 
#define AZ_loss             4 
#define AZ_ill_cond         5 



#define AZ_solver              0
#define AZ_scaling             1
#define AZ_precond             2
#define AZ_conv                3
#define AZ_output              4
#define AZ_pre_calc            5
#define AZ_max_iter            6
#define AZ_poly_ord            7
#define AZ_overlap             8
#define AZ_kspace              9
#define AZ_orthog              10
#define AZ_aux_vec             11
#define AZ_print_freq          12



#define AZ_tol                 0
#define AZ_drop                2
#define AZ_weights             3



#define AZ_matrix_type         0
#define AZ_N_internal          1
#define AZ_N_border            2
#define AZ_N_external          3
#define AZ_N_int_blk           4
#define AZ_N_bord_blk          5
#define AZ_N_ext_blk           6
#define AZ_N_neigh             7
#define AZ_total_send          8
#define AZ_name                9
#define AZ_neighbors           10
#define AZ_rec_length          (10 +   AZ_MAX_NEIGHBORS)
#define AZ_send_length         (10 + 2*AZ_MAX_NEIGHBORS)
#define AZ_send_list           (10 + 3*AZ_MAX_NEIGHBORS)



#define AZ_its                 0
#define AZ_why                 1
#define AZ_r                   2
#define AZ_rec_r               3
#define AZ_scaled_r            4
#define AZ_first_precond       5     
                                     
                                     
                                     
                                     



#define AZ_node                0
#define AZ_N_procs             1
#define AZ_dim                 2



#define AZ_linear              0
#define AZ_file                1
#define AZ_box                 2



#define AZ_ALLOC               0
#define AZ_CLEAR               1
#define AZ_REALLOC             2
#define AZ_SYS                 -14901
#define AZ_OLD_ADDRESS         0
#define AZ_NEW_ADDRESS         1



#define AZ_MSR_MATRIX          0
#define AZ_VBR_MATRIX          1
#define AZ_USER_MATRIX         2



#define AZ_SCALE_MAT           0
#define AZ_RESCALE_SOL         1



#define AZ_NOT_FIRST           0   
                                   
                                   
                                   
#define AZ_FIRST_TIME          1   
                                   
                                   
                                   
                                   
                                   



#define AZ_QUIT             5
#define AZ_CONTINUE         6



#define AZ_PACK             0
#define AZ_SEND             1




#define AZ_TEST_ELE         3
#define AZ_ALL              1 
#define AZ_EXTERNS          2 
#define AZ_GLOBAL           1 
#define AZ_LOCAL            2 
