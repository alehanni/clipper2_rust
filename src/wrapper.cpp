#include "clipper2/clipper.h"
#include "clipper2/clipper.export.h"

// extra wrappers
namespace Clipper2Lib {
    
EXTERN_DLL_EXPORT CPaths64 SimplifyPaths64(const CPaths64 paths, double epsilon, bool is_open_path=false) {
    Paths64 pp = ConvertCPaths(paths);
    Paths64 result = SimplifyPaths(pp, epsilon, is_open_path);
    return CreateCPaths(result);
}

}