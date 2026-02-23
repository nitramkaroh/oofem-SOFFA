#ifdef DEBUG
#pragma push_macro( "DEBUG" )
#undef DEBUG
#define RESTORE_DEBUG_MACRO
#endif

#define SYM_ILDL_DEBUG

#include <solver.h>

#ifdef RESTORE_DEBUG_MACRO
#pragma pop_macro( "DEBUG" )
#undef RESTORE_DEBUG_MACRO
#endif