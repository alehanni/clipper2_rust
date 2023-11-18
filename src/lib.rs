use libc::{c_long, size_t, c_double};
use thiserror::Error;

mod ffi {
    use libc::{c_int, c_long, c_double, c_char};

    pub type CPaths64 = *const c_long;
    pub type CPaths64Mut = *mut c_long;
    pub type CPaths64Ref = *mut *mut c_long;

    #[repr(C)]
    pub struct CRect64{
        pub left: c_long,
        pub top: c_long,
        pub right: c_long,
        pub bottom: c_long,
    }

    impl From<[i64; 4]> for CRect64 {
        fn from(slice: [i64; 4]) -> CRect64 {
            CRect64 { left: slice[0], top: slice[1], right: slice[2], bottom: slice[3] }
        }
    }

    impl From<(i64, i64, i64, i64)> for CRect64 {
        fn from(tup: (i64, i64, i64, i64)) -> CRect64 {
            CRect64 { left: tup.0, top: tup.1, right: tup.2, bottom: tup.3 }
        }
    }

    extern "C" {
        pub fn Version() -> *const c_char;

        pub fn BooleanOp64(
            cliptype: u8,
            fillrule: u8,
            subjects: CPaths64,
            subjects_open: CPaths64,
            clips: CPaths64,
            solution: CPaths64Ref,
            solution_open: CPaths64Ref,
            preserve_collinear: bool,
            reverse_solution: bool,
        ) -> c_int;
        
        pub fn InflatePaths64(
            paths: CPaths64,
            delta: c_double,
            jointype: u8,
            endtype: u8,
            miter_limit: c_double,
            arc_tolerance: c_double,
            reverse_solution: bool,
        ) -> CPaths64Mut;

        pub fn RectClip64(
            rect: *const CRect64,
            paths: CPaths64,
        ) -> CPaths64Mut;

        pub fn RectClipLines64(
            rect: *const CRect64,
            paths: CPaths64,
        ) -> CPaths64Mut;

        pub fn MinkowskiSum64(
            pattern: CPaths64,
            path: CPaths64,
            is_closed: bool,
        ) -> CPaths64Mut;

        pub fn MinkowskiDiff64(
            pattern: CPaths64,
            path: CPaths64,
            is_closed: bool,
        ) -> CPaths64Mut;

        pub fn SimplifyPaths64(
            paths: CPaths64,
            epsilon: c_double,
            is_open_path: bool,
        ) -> CPaths64Mut;

        pub fn TrimCollinear64(
            paths: CPaths64,
            is_open_path: bool,
        ) -> CPaths64Mut;

        pub fn DisposeExportedCPaths64(p: CPaths64Ref);
    }
}

pub enum ClipType {
    Intersection = 1,
    Union,
    Difference,
    Xor,
}

pub enum FillRule {
    EvenOdd,
    NonZero,
    Positive,
    Negative,
}

pub enum JoinType {
    Square,
    Round,
    Miter,
}

pub enum EndType {
    Polygon,
    Joined,
    Butt,
    Square,
    Round,
}

// NOTE:
// in order to interface with the exported clipper functions
// we need to comply with its data structures
//
// Path64 data structure
// __________________________________
// |counter|coord1|coord2|...|coordN|
// |N, 0   |x1, y1|x2, y2|...|xN, yN|
// __________________________________
//
// N: int64_t <=> number of coords
// xN, yN: int64_t <=> said coords
//
// CPaths64 data structure
// _______________________________
// |counter|path1|path2|...|pathN|
// |len, N |     |     |...|pathN|
// _______________________________
//
// len: int64_t <=> array length
// N: int64_t <=> number of paths
//
// here, path1, path2, etc. aren't pointers to the structures
// but rather represent the paths packed one after the other 
// in a single contiguous int64_t array

fn cpaths_from_vec<T>(paths: &Vec<Vec<T>>) -> Vec<c_long>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let n_paths = paths.len();
    let n_coords = paths.iter().map(|x| x.len()).sum::<usize>();
    let buffer_len = 2 * (n_paths + n_coords + 1);

    let mut cpaths_buffer = vec![0 as c_long; buffer_len];

    // fill buffer with appropriate data
    cpaths_buffer[0] = buffer_len as i64;
    cpaths_buffer[1] = n_paths as i64;
    
    let mut i = 2;
    for path in paths {
        cpaths_buffer[i] = path.len() as i64;
        cpaths_buffer[i+1] = 0;
        i += 2;

        for coord in path {
            let coord: [c_long; 2] = (*coord).into();
            cpaths_buffer[i] = coord[0];
            cpaths_buffer[i+1] = coord[1];
            i += 2;
        }
    }

    return cpaths_buffer;
}

// single-path special case of cpaths_from_vec
fn cpaths_from_coords<T>(coords: &Vec<T>) -> Vec<c_long>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let n_coords = coords.len();
    let buffer_len = 2 * (n_coords + 2);

    let mut cpaths_buffer = vec![0 as c_long; buffer_len];

    cpaths_buffer[0] = buffer_len as i64;
    cpaths_buffer[1] = 1; // only 1 path
    cpaths_buffer[2] = n_coords as i64;
    cpaths_buffer[3] = 0;

    let mut i = 4;
    for coord in coords {
        let coord: [c_long; 2] = (*coord).into();
        cpaths_buffer[i] = coord[0];
        cpaths_buffer[i+1] = coord[1];
        i += 2;
    }

    return cpaths_buffer;
}

fn vec_from_cpaths<T>(cpaths_buffer: &Vec<c_long>) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let n_paths = cpaths_buffer[1];

    let mut cpaths_iter = cpaths_buffer.into_iter();
    cpaths_iter.next(); // skip past the metadata
    cpaths_iter.next(); //

    let mut paths64 = Vec::<Vec<T>>::with_capacity(n_paths as size_t);
    for _ in 0..n_paths {
        let path_len = *cpaths_iter.next().unwrap();
        let zero = *cpaths_iter.next().unwrap(); // this should always be 0
        debug_assert_eq!(0, zero);

        let mut path = Vec::<T>::with_capacity(2*path_len as size_t);
        for _ in 0..path_len { // change when https://github.com/rust-lang/rust/issues/98326 has been closed?
            let point = [
                *cpaths_iter.next().unwrap(), // x
                *cpaths_iter.next().unwrap(), // y
            ];
            path.push(T::from(point));
        }
        paths64.push(path);
    }

    return paths64;
}

unsafe fn vec_from_raw_cpaths<T>(cpaths_ptr: *mut *mut c_long) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let buffer_len = **cpaths_ptr as size_t;
    let cpaths_buffer = Vec::from_raw_parts(*cpaths_ptr, buffer_len, buffer_len); // this vec can not be re-alloc'd or freed with rust's allocator

    let result = vec_from_cpaths::<T>(&cpaths_buffer);
    
    std::mem::forget(cpaths_buffer); // must be leaked since it holds a pointer to c++-allocated memory...
    ffi::DisposeExportedCPaths64(cpaths_ptr); // ...which is properly freed here instead

    return result;
}

#[derive(Error, Debug)]
pub enum BooleanOpError {
    #[error("argument \"cliptype\" is not a valid ClipType enum")]
    InvalidClipType,
    #[error("argument \"fillrule\" is not a valid FillRule enum")]
    InvalidFillRule,
    #[error("clipper encountered an error during execute")]
    ClipperExecuteError,
    #[error("unknown clipper error")]
    Unknown,
}

pub fn boolean_op<T>(
    subjects: &Vec<Vec<T>>,
    clips: &Vec<Vec<T>>,
    fillrule: FillRule,
    cliptype: ClipType,
    preserve_collinear: bool,
    reverse_solution: bool,
) -> Result<Vec<Vec<T>>, BooleanOpError>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let subj_buffer = cpaths_from_vec(subjects);
    let clip_buffer = cpaths_from_vec(clips);
    let sol_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();      // the foreign function sets these to point to allocated memory
    let sol_open_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut(); //

    let errcode = unsafe { 
        ffi::BooleanOp64(
            cliptype as u8,
            fillrule as u8,
            subj_buffer.as_ptr(),
            std::ptr::null_mut(), // a null pointer is interpreted as an empty set of paths
            clip_buffer.as_ptr(),
            sol_buffer,
            sol_open_buffer, // unused, but should point to an empty CPaths64 that needs freeing* (on errcode==0)
            preserve_collinear,
            reverse_solution,
        )
    };

    if errcode != 0 {
        // sol(_open)_buffer shouldn't point to anything and thus don't require clean up
        match errcode {
            -4 => return Err(BooleanOpError::InvalidClipType),
            -3 => return Err(BooleanOpError::InvalidFillRule),
            -1 => return Err(BooleanOpError::ClipperExecuteError),
            _ => return Err(BooleanOpError::Unknown),
        }
    }

    let solution = unsafe { vec_from_raw_cpaths(sol_buffer) };
    unsafe { ffi::DisposeExportedCPaths64(sol_open_buffer); } // *freed here

    return Ok(solution);
}

pub fn intersect<T>(subjects: &Vec<Vec<T>>, clips: &Vec<Vec<T>>, fillrule: FillRule) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    match boolean_op(subjects, clips, fillrule, ClipType::Intersection, true, false) {
        Ok(res) => res,
        Err(e) => panic!("intersect returned an error: {}", e),
    }
}

pub fn union<T>(subjects: &Vec<Vec<T>>, clips: &Vec<Vec<T>>, fillrule: FillRule) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    match boolean_op(subjects, clips, fillrule, ClipType::Union, true, false) {
        Ok(res) => res,
        Err(e) => panic!("union returned an error: {}", e),
    }
}

pub fn difference<T>(subjects: &Vec<Vec<T>>, clips: &Vec<Vec<T>>, fillrule: FillRule) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    match boolean_op(subjects, clips, fillrule, ClipType::Difference, true, false) {
        Ok(res) => res,
        Err(e) => panic!("difference returned an error: {}", e),
    }
}

pub fn xor<T>(subjects: &Vec<Vec<T>>, clips: &Vec<Vec<T>>, fillrule: FillRule) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    match boolean_op(subjects, clips, fillrule, ClipType::Xor, true, false) {
        Ok(res) => res,
        Err(e) => panic!("xor returned an error: {}", e),
    }
}

pub fn inflate_paths_ext<T>(
    paths: &Vec<Vec<T>>,
    delta: f64,
    jointype: JoinType,
    endtype: EndType,
    miter_limit: f64,
    arc_tolerance: f64,
    reverse_solution: bool
) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let input_buffer = cpaths_from_vec(paths);
    let result_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();

    let result = unsafe {
        *result_buffer = ffi::InflatePaths64(
            input_buffer.as_ptr(),
            delta as c_double,
            jointype as u8,
            endtype as u8,
            miter_limit as c_double,
            arc_tolerance as c_double,
            reverse_solution,
        );
        vec_from_raw_cpaths(result_buffer)
    };

    return result;
}

pub fn inflate_paths<T>(paths: &Vec<Vec<T>>, delta: f64, jointype: JoinType, endtype: EndType) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    inflate_paths_ext(paths, delta, jointype, endtype, 2.0, 0.0, false)
}

#[derive(Error, Debug)]
pub enum RectClipError {
    #[error("the provided rectangle has 0 or negative area")]
    CRectIsEmpty,
    #[error("expected paths, found null")]
    NullPaths,
}

pub fn rect_clip_ext<R, T>(rect: &R, paths: &Vec<Vec<T>>) -> Result<Vec<Vec<T>>, RectClipError>
where
    R: Into<[c_long; 4]> + From<[c_long; 4]> + Copy,
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let crect = ffi::CRect64::from((*rect).into());

    let input_buffer = cpaths_from_vec(paths);
    let result_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();

    let result = unsafe {
        *result_buffer = ffi::RectClip64(
            &crect,
            input_buffer.as_ptr(),
        );

        // a nullptr is returned if an error occurred
        if (*result_buffer).is_null() {
            if (crect.right <= crect.left) || (crect.bottom <= crect.top) {
                return Err(RectClipError::CRectIsEmpty);
            } else {
                return Err(RectClipError::NullPaths); // this may not even be reachable
            }
        }

        vec_from_raw_cpaths(result_buffer)
    };

    return Ok(result);
}

pub fn rect_clip<R, T>(rect: &R, paths: &Vec<Vec<T>>) -> Vec<Vec<T>>
where
    R: Into<[c_long; 4]> + From<[c_long; 4]> + Copy,
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    match rect_clip_ext(rect, paths) {
        Ok(res) => res,
        Err(e) => panic!("rect_clip_ext returned an error: {}", e),
    }
}

pub fn rect_clip_lines_ext<R, T>(rect: &R, paths: &Vec<Vec<T>>) -> Result<Vec<Vec<T>>, RectClipError>
where
    R: Into<[c_long; 4]> + From<[c_long; 4]> + Copy,
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let crect = ffi::CRect64::from((*rect).into());

    let input_buffer = cpaths_from_vec(paths);
    let result_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();

    let result = unsafe {
        *result_buffer = ffi::RectClipLines64(
            &crect,
            input_buffer.as_ptr(),
        );

        // a nullptr is returned if an error occurred
        if (*result_buffer).is_null() {
            if (crect.right <= crect.left) || (crect.bottom <= crect.top) {
                return Err(RectClipError::CRectIsEmpty);
            } else {
                return Err(RectClipError::NullPaths);
            }
        }

        vec_from_raw_cpaths(result_buffer)
    };

    return Ok(result);
}

pub fn rect_clip_lines<R, T>(rect: &R, paths: &Vec<Vec<T>>) -> Vec<Vec<T>>
where
    R: Into<[c_long; 4]> + From<[c_long; 4]> + Copy,
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    match rect_clip_lines_ext(rect, paths) {
        Ok(res) => res,
        Err(e) => panic!("rect_clip_lines_ext returned an error: {}", e),
    }
}

pub fn minkowski_sum<T>(pattern: &Vec<T>, path: &Vec<T>, is_closed: bool) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let pattern_buffer = cpaths_from_coords(pattern); // the original library api uses a single pattern/path
    let path_buffer = cpaths_from_coords(path);       // but the exported functions only deal with paths (plural)
    let result_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();

    let result = unsafe {
        *result_buffer = ffi::MinkowskiSum64(
            pattern_buffer.as_ptr(),
            path_buffer.as_ptr(),
            is_closed
        );
        vec_from_raw_cpaths(result_buffer)
    };

    return result;
}

pub fn minkowski_diff<T>(pattern: &Vec<T>, path: &Vec<T>, is_closed: bool) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let pattern_buffer = cpaths_from_coords(pattern);
    let path_buffer = cpaths_from_coords(path);
    let result_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();

    let result = unsafe {
        *result_buffer = ffi::MinkowskiDiff64(
            pattern_buffer.as_ptr(),
            path_buffer.as_ptr(),
            is_closed
        );
        vec_from_raw_cpaths(result_buffer)
    };

    return result;
}

pub fn simplify_paths_ext<T>(paths: &Vec<Vec<T>>, epsilon: f64, is_open_path: bool) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let input_buffer = cpaths_from_vec(paths);
    let result_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();

    let result = unsafe {
        *result_buffer = ffi::SimplifyPaths64(
            input_buffer.as_ptr(),
            epsilon,
            is_open_path,
        );
        vec_from_raw_cpaths(result_buffer)
    };

    return result;
}

pub fn simplify_paths<T>(paths: &Vec<Vec<T>>, epsilon: f64) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    simplify_paths_ext(paths, epsilon, false)
}

pub fn trim_collinear_ext<T>(path: &Vec<T>, is_open_path: bool) -> Vec<T>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let input_buffer = cpaths_from_coords(path);
    let result_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();

    let mut result = unsafe {
        *result_buffer = ffi::TrimCollinear64(
            input_buffer.as_ptr(),
            is_open_path,
        );
        vec_from_raw_cpaths(result_buffer)
    };

    return result.pop().unwrap();
}

pub fn trim_collinear<T>(path: &Vec<T>) -> Vec<T>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    trim_collinear_ext(path, false)
}

#[cfg(test)]
mod tests {

    use std::f64::consts::PI;
    use crate::*;
    use glam::I64Vec2;

    #[test]
    fn boolean_op_test() {
        let box1_subj = vec![vec![
            I64Vec2::new(48, 48),
            I64Vec2::new(48, -16),
            I64Vec2::new(-16, -16),
            I64Vec2::new(-16, 48),
        ]];

        let box2_clip = vec![vec![
            I64Vec2::new(-48, -48),
            I64Vec2::new(-48, 16),
            I64Vec2::new(16, 16),
            I64Vec2::new(16, -48),
        ]];

        let res = intersect(&box1_subj, &box2_clip, FillRule::EvenOdd);
        assert_eq!(res, [[I64Vec2::new(16, -16), I64Vec2::new(16, 16), I64Vec2::new(-16, 16), I64Vec2::new(-16, -16)]]);
    }

    #[test]
    fn invalid_cliptype_test() {
        let box1_subj = vec![vec![I64Vec2::new(48, 48), I64Vec2::new(48, -16), I64Vec2::new(-16, -16), I64Vec2::new(-16, 48)]];
        let box2_clip = vec![vec![I64Vec2::new(-48, -48), I64Vec2::new(-48, 16), I64Vec2::new(16, 16), I64Vec2::new(16, -48)]];
        let invalid_cliptype: ClipType  = unsafe { std::mem::transmute(100 as u8) };
        assert!(boolean_op(&box1_subj, &box2_clip, FillRule::EvenOdd, invalid_cliptype, true, false).is_err());
    }

    #[test]
    fn invalid_fillrule_test() {
        let box1_subj = vec![vec![I64Vec2::new(48, 48), I64Vec2::new(48, -16), I64Vec2::new(-16, -16), I64Vec2::new(-16, 48)]];
        let box2_clip = vec![vec![I64Vec2::new(-48, -48), I64Vec2::new(-48, 16), I64Vec2::new(16, 16), I64Vec2::new(16, -48)]];
        let invalid_fillrule: FillRule  = unsafe { std::mem::transmute(100 as u8) };
        assert!(boolean_op(&box1_subj, &box2_clip, invalid_fillrule, ClipType::Intersection, true, false).is_err());
    }

    #[test]
    fn inflate_paths_test() {
        let paths = vec![vec![
            I64Vec2::new(100, 100),
            I64Vec2::new(1500, 100),
            I64Vec2::new(100, 1500),
            I64Vec2::new(1500, 1500),
        ]];

        let mut res = inflate_paths(&paths, 200.0, JoinType::Miter, EndType::Square);
        res = simplify_paths(&res, 100.0);
        assert_eq!(res[0].len(), 10);
    }

    #[test]
    fn rect_clip_test() {
        let radius = 100.0 * 2.0_f64.sqrt();
        let n_points = 32;
        let circle_x = |i| radius * (i as f64/n_points as f64 * 2.0*PI).cos();
        let circle_y = |i| radius * (i as f64/n_points as f64 * 2.0*PI).sin();

        let circle_vertices: Vec<Vec<I64Vec2>> = vec![(0..n_points).map(|i| I64Vec2::new(circle_x(i) as i64, circle_y(i) as i64)).collect()];
        let rect = [-50, -50, 50, 50];

        let res = rect_clip(&rect, &circle_vertices);
        assert_eq!(res[0].len(), 4);
    }

    #[test]
    fn rect_clip_empty_test() {
        let subj = vec![vec![I64Vec2::new(48, 48), I64Vec2::new(48, -16), I64Vec2::new(-16, -16), I64Vec2::new(-16, 48)]];
        let rect = [0, 0, 0, 0];
        assert!(rect_clip_ext(&rect, &subj).is_err());
    }

    #[test]
    fn rect_clip_lines_empty_test() {
        let subj = vec![vec![I64Vec2::new(48, 48), I64Vec2::new(48, -16), I64Vec2::new(-16, -16), I64Vec2::new(-16, 48)]];
        let rect = [0, 0, 0, 0];
        assert!(rect_clip_lines_ext(&rect, &subj).is_err());
    }

    #[test]
    fn minkowski_sum_closed_test() {
        let box1_pattern = vec![I64Vec2::new(20, 20), I64Vec2::new(20, 0), I64Vec2::new(0, 0), I64Vec2::new(0, 20)];
        let box2_path = vec![I64Vec2::new(10, 10), I64Vec2::new(10, 40), I64Vec2::new(40, 40), I64Vec2::new(40, 10)];
        let res = minkowski_sum(&box1_pattern, &box2_path, true);
        assert_eq!(trim_collinear(&res[0]), [I64Vec2::new(60, 10), I64Vec2::new(60, 60), I64Vec2::new(10, 60), I64Vec2::new(10, 10)]);
        assert_eq!(trim_collinear(&res[1]), [I64Vec2::new(30, 40), I64Vec2::new(40, 40), I64Vec2::new(40, 30), I64Vec2::new(30, 30)]);
    }

    #[test]
    fn minkowski_diff_closed_test() {
        let box1_pattern = vec![I64Vec2::new(20, 20), I64Vec2::new(20, 0), I64Vec2::new(0, 0), I64Vec2::new(0, 20)];
        let box2_path = vec![I64Vec2::new(10, 10), I64Vec2::new(10, 40), I64Vec2::new(40, 40), I64Vec2::new(40, 10)];
        let res = minkowski_diff(&box1_pattern, &box2_path, true);
        assert_eq!(trim_collinear(&res[0]), [I64Vec2::new(40, -10), I64Vec2::new(40, 40), I64Vec2::new(-10, 40), I64Vec2::new(-10, -10)]);
        assert_eq!(trim_collinear(&res[1]), [I64Vec2::new(10, 20), I64Vec2::new(20, 20), I64Vec2::new(20, 10), I64Vec2::new(10, 10)]);
    }
}