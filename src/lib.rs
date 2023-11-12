#![allow(unused_must_use, dead_code)]
use libc::{c_long, size_t};
use std::vec::Vec;

mod ffi {
    use libc::{c_int, c_long, c_double, c_char};

    // NOTES:
    // removed c_double variants due to them only adding a step that converts to int64_t for the actual computations

    pub type CPaths64 = *const c_long;
    pub type CPaths64Ref = *mut *mut c_long;

    #[repr(C)]
    pub struct CRect64{
        left: c_long,
        top: c_long,
        right: c_long,
        bottom: c_long,
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
        ) -> CPaths64;

        pub fn RectClip64(
            rect: &CRect64,
            paths: CPaths64,
        ) -> CPaths64;

        pub fn RectClipLines64(
            rect: &CRect64,
            paths: CPaths64,
        ) -> CPaths64;

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
// N: size_t <=> number of coords
// xN, yN: int64_t <=> said coords
//
// CPaths64 data structure
// _______________________________
// |counter|path1|path2|...|pathN|
// |len, N |     |     |...|pathN|
// _______________________________
//
// len: size_t <=> array length
// N: size_t <=> number of paths
//
// here, path1, path2, etc. aren't pointers to the structures
// but rather represent the paths packed one after the other 
// in a single contiguous int64_t array

fn cpaths_from_vec<T>(paths: Vec<Vec<T>>) -> Vec<c_long>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy,
{
    let n_paths = paths.len();
    let n_coords = paths.iter().map(|x| x.len()).sum::<usize>();
    let buffer_len = 2 * (n_paths + n_coords + 1);

    let mut cpaths_buffer = vec![0 as c_long; buffer_len];

    // fill buffer with data
    cpaths_buffer[0] = buffer_len as i64;
    cpaths_buffer[1] = n_paths as i64;
    
    let mut i = 2;
    for path in paths {
        cpaths_buffer[i] = path.len() as i64;
        cpaths_buffer[i+1] = 0;
        i += 2;

        for coord in path {
            let coord: [c_long; 2] = coord.into();
            cpaths_buffer[i] = coord[0];
            cpaths_buffer[i+1] = coord[1];
            i += 2;
        }
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
        debug_assert_eq!(0, *cpaths_iter.next().unwrap()); // this should always be 0

        let mut path = Vec::<T>::with_capacity(2*path_len as size_t);
        for _ in 0..path_len { // change when https://github.com/rust-lang/rust/issues/98326 has been closed
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

fn boolean_op<T>(subjects: Vec<Vec<T>>, clips: Vec<Vec<T>>, fillrule: FillRule, cliptype: ClipType) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy + std::fmt::Debug,
{
    let subj_buffer = cpaths_from_vec(subjects);
    let clip_buffer = cpaths_from_vec(clips);
    let sol_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut();      // set to point to allocated memory by the foreign function
    let sol_open_buffer: ffi::CPaths64Ref = &mut std::ptr::null_mut(); //

    unsafe {
        ffi::BooleanOp64(
            cliptype as u8,
            fillrule as u8,
            subj_buffer.as_ptr(),
            std::ptr::null_mut(), // a null pointer is interpreted as an empty set of paths
            clip_buffer.as_ptr(),
            sol_buffer,
            sol_open_buffer, // unused, but should point to an empty CPaths64 that needs freeing
            true,
            false,
        );
    }

    // sol_buffer and sol_open_buffer are set by BooleanOp64 to point to the result
    let buffer_len = unsafe { **sol_buffer } as size_t; // first entry in allocated memory should return length
    let cpaths_buffer = unsafe { Vec::from_raw_parts(*sol_buffer, buffer_len, buffer_len) };
    
    let solution = vec_from_cpaths::<T>(&cpaths_buffer);
    
    // free the buffers
    std::mem::forget(cpaths_buffer); // must be leaked since it holds a pointer to c++-allocated memory
    unsafe {
        ffi::DisposeExportedCPaths64(sol_buffer);
        ffi::DisposeExportedCPaths64(sol_open_buffer);
    }

    return solution;
}

pub fn intersect<T>(subjects: Vec<Vec<T>>, clips: Vec<Vec<T>>, fillrule: FillRule) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy + std::fmt::Debug,
{
    boolean_op(subjects, clips, fillrule, ClipType::Intersection)
}

pub fn union<T>(subjects: Vec<Vec<T>>, clips: Vec<Vec<T>>, fillrule: FillRule) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy + std::fmt::Debug,
{
    boolean_op(subjects, clips, fillrule, ClipType::Union)
}

pub fn difference<T>(subjects: Vec<Vec<T>>, clips: Vec<Vec<T>>, fillrule: FillRule) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy + std::fmt::Debug,
{
    boolean_op(subjects, clips, fillrule, ClipType::Difference)
}

pub fn xor<T>(subjects: Vec<Vec<T>>, clips: Vec<Vec<T>>, fillrule: FillRule) -> Vec<Vec<T>>
where
    T: Into<[c_long; 2]> + From<[c_long; 2]> + Copy + std::fmt::Debug,
{
    boolean_op(subjects, clips, fillrule, ClipType::Xor)
}

#[cfg(test)]
mod tests {

    use crate::intersect;
    use crate::FillRule;
    use glam::I64Vec2;

    #[test]
    fn basic_test() {
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

        let res = intersect(box1_subj, box2_clip, FillRule::EvenOdd);
        //println!("res: {:?}", res);
    }
}