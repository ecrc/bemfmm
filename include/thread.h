//#ifndef THREAD_ACOUSTICS_H
//#define THREAD_ACOUSTICS_H
//#if EXAFMM_WITH_TBB
//#define num_threads(E)                tbb::task_scheduler_init init(E)
//#include <tbb/task_scheduler_init.h>
//#elif EXAFMM_WITH_CILK
//#include <cilk/cilk.h>
//#include <cilk/cilk_api.h>
//#define num_threads(E)                char nworkers[32]; sprintf(nworkers,"%d",E); __cilkrts_set_param("nworkers",nworkers)
//#elif EXAFMM_WITH_MTHREAD
//#define num_threads(E)                mtbb::task_scheduler_init init(E)
//#include <myth/myth.h>
//#include <mtbb/task_scheduler_init.h>
//#elif EXAFMM_WITH_OMP
//#include <omp.h>
//#define num_threads(E)                omp_set_num_threads(E);
//#else 
//#define num_threads(E) 
//#endif
//#endif
