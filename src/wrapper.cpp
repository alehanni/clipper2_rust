#include "clipper2/clipper.h"
#include "clipper2/clipper.export.h"

// extra wrappers
namespace Clipper2Lib {
    
EXTERN_DLL_EXPORT CPaths64 SimplifyPaths64(const CPaths64 paths, double epsilon, bool is_open_path=false) {
    Paths64 pp = ConvertCPaths(paths);
    Paths64 result = SimplifyPaths(pp, epsilon, is_open_path);
    return CreateCPaths(result);
}

EXTERN_DLL_EXPORT CPaths64 TrimCollinear64(const CPaths64 path, bool is_open_path=false) {
    Paths64 pp = ConvertCPaths(path);
    Paths64 result;
    result.push_back(TrimCollinear(pp[0], is_open_path));
    return CreateCPaths(result);
}

EXTERN_DLL_EXPORT CPaths64 MinkowskiSum64(const CPaths64 pattern, const CPaths64 path, bool is_closed) {
    Paths64 ppattern = ConvertCPaths(pattern);
    Paths64 ppath = ConvertCPaths(path);
    Paths64 result = MinkowskiSum(ppattern[0], ppath[0], is_closed);
    return CreateCPaths(result);
}

EXTERN_DLL_EXPORT CPaths64 MinkowskiDiff64(const CPaths64 pattern, const CPaths64 path, bool is_closed) {
    Paths64 ppattern = ConvertCPaths(pattern);
    Paths64 ppath = ConvertCPaths(path);
    Paths64 result = MinkowskiDiff(ppattern[0], ppath[0], is_closed);
    return CreateCPaths(result);
}

}