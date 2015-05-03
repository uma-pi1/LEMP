# define source files


set(ta_SRCS
# 1. structs

structs/read.cpp


# ap/cmdline.c
# ap/csr.c
# ap/memory.c
# ap/select.c
# ap/sort.c
# ap/util.c
# 
# #classifier/BayesianClassifier.cpp
# 
# 
# #my_mlpack/core/tree/mrkd_statistic.cpp
# # 2. algos
# 
# 
#   my_mlpack/core/arma_extend/arma_extend.hpp
#   my_mlpack/core/arma_extend/fn_ccov.hpp
#   my_mlpack/core/arma_extend/glue_ccov_meat.hpp
#   my_mlpack/core/arma_extend/glue_ccov_proto.hpp
#   my_mlpack/core/arma_extend/hdf5_misc.hpp
#   my_mlpack/core/arma_extend/op_ccov_meat.hpp
#   my_mlpack/core/arma_extend/op_ccov_proto.hpp
#   my_mlpack/core/arma_extend/promote_type.hpp
#   my_mlpack/core/arma_extend/restrictors.hpp
#   my_mlpack/core/arma_extend/traits.hpp
#   my_mlpack/core/arma_extend/typedef.hpp
#   my_mlpack/core/arma_extend/SpMat_extra_bones.hpp
#   my_mlpack/core/arma_extend/SpMat_extra_meat.hpp
# 
#   my_mlpack/core/util/cli.hpp
#   my_mlpack/core/util/cli.cpp
#   my_mlpack/core/util/cli_deleter.hpp
#   my_mlpack/core/util/cli_deleter.cpp
#   my_mlpack/core/util/cli_impl.hpp
#   my_mlpack/core/util/log.hpp
#   my_mlpack/core/util/log.cpp
#   my_mlpack/core/util/nulloutstream.hpp
#   my_mlpack/core/util/option.hpp
#   my_mlpack/core/util/option.cpp
#   my_mlpack/core/util/option_impl.hpp
#   my_mlpack/core/util/prefixedoutstream.hpp
#   my_mlpack/core/util/prefixedoutstream.cpp
#   my_mlpack/core/util/prefixedoutstream_impl.hpp
#   my_mlpack/core/util/sfinae_utility.hpp
#   my_mlpack/core/util/string_util.hpp
#   my_mlpack/core/util/string_util.cpp
#   my_mlpack/core/util/timers.hpp
#   my_mlpack/core/util/timers.cpp
#   my_mlpack/core/util/version.hpp
#   my_mlpack/core/util/version.cpp
#    
#   my_mlpack/core/kernels/cosine_distance.hpp
#   my_mlpack/core/kernels/cosine_distance_impl.hpp
#   my_mlpack/core/kernels/example_kernel.hpp
#   my_mlpack/core/kernels/kernel_traits.hpp
#   my_mlpack/core/kernels/linear_kernel.hpp  
#   
#   my_mlpack/core/math/clamp.hpp
#   my_mlpack/core/math/random.hpp
#   my_mlpack/core/math/random.cpp
#   my_mlpack/core/math/range.hpp
#   my_mlpack/core/math/range_impl.hpp
#   my_mlpack/core/math/round.hpp
#   
#   my_mlpack/core/metrics/ip_metric.hpp
#   my_mlpack/core/metrics/ip_metric_impl.hpp
#   my_mlpack/core/metrics/lmetric.hpp
#   my_mlpack/core/metrics/lmetric_impl.hpp  
#   my_mlpack/core/tree/ballbound.hpp
#   my_mlpack/core/tree/ballbound_impl.hpp
#   my_mlpack/core/tree/bounds.hpp
#   my_mlpack/core/tree/cover_tree/cover_tree.hpp
#   my_mlpack/core/tree/cover_tree/cover_tree_impl.hpp
#   my_mlpack/core/tree/cover_tree/first_point_is_root.hpp
#   my_mlpack/core/tree/cover_tree/single_tree_traverser.hpp
#   my_mlpack/core/tree/cover_tree/single_tree_traverser_impl.hpp
#   my_mlpack/core/tree/cover_tree/traits.hpp
#   my_mlpack/core/tree/hrectbound.hpp
#   my_mlpack/core/tree/hrectbound_impl.hpp
#   my_mlpack/core/tree/periodichrectbound.hpp
#   my_mlpack/core/tree/periodichrectbound_impl.hpp
#   my_mlpack/core/tree/statistic.hpp
#   my_mlpack/core/tree/tree_traits.hpp  
)

# define main header files
set(ta_HDRS
#         ap/includes.h  
 	taLib.h

)



# combine to all headers 
set(ta_all_HDRS
	${ta_HDRS}
)

# add library targets
add_library(ta SHARED ${ta_all_HDRS} ${ta_SRCS}) 

# link with boost
target_link_libraries( ta rg ${Boost_LIBRARIES})
# target_link_libraries( ta rg armadillo ${Boost_LIBRARIES} blas lapack ${GSL_LIBRARIES} GKlib)

# install library and main headers
install(TARGETS ta DESTINATION ${LIB_INSTALL_DIR})
install(FILES ${ta_HDRS} DESTINATION ${INCLUDE_INSTALL_DIR}/taLib)



