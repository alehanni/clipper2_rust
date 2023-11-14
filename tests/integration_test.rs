use clipper2_rust::{intersect, union, difference, xor, inflate_paths, simplify_paths, JoinType, EndType};
use clipper2_rust::FillRule;
use glam::{I64Vec2, DVec2};
use plotters::prelude::*;
use regex::Regex;
use svg_path_parser;

fn ivec_to_cartesian(input: Vec<I64Vec2>) -> Vec<(f64, f64)> {
    input.into_iter().map(|x| x.as_dvec2().into()).collect()
}

#[test]
fn plot_intersection() {
    let root = SVGBackend::new("test_plots/0_intersection_test.svg", (500, 500)).into_drawing_area();

    root.fill(&TRANSPARENT).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .build_cartesian_2d(0.0..100.0, 0.0..100.0).unwrap();

    let subj_vertices = vec![vec![
        I64Vec2::new(100, 50),
        I64Vec2::new(10, 79),
        I64Vec2::new(65, 2),
        I64Vec2::new(65, 98),
        I64Vec2::new(10, 21),
    ]];

    let clip_vertices = vec![vec![
        I64Vec2::new(98, 63),
        I64Vec2::new(4, 68),
        I64Vec2::new(77, 8),
        I64Vec2::new(52, 100),
        I64Vec2::new(19, 12),
    ]];

    let intersection_vertices = intersect(subj_vertices.clone(), clip_vertices.clone(), FillRule::NonZero);
    println!("{:?}", intersection_vertices);

    chart.draw_series(std::iter::once(Polygon::new(
        ivec_to_cartesian(intersection_vertices[0].clone()),
        RGBColor(255, 0, 255),
    ))).unwrap();

    root.present().unwrap();
}

#[test]
fn plot_union() {
    let root = SVGBackend::new("test_plots/1_union_test.svg", (500, 500)).into_drawing_area();

    root.fill(&TRANSPARENT).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .build_cartesian_2d(0.0..100.0, 0.0..100.0).unwrap();

    let subj_vertices = vec![vec![
        I64Vec2::new(100, 50),
        I64Vec2::new(10, 79),
        I64Vec2::new(65, 2),
        I64Vec2::new(65, 98),
        I64Vec2::new(10, 21),
    ]];

    let clip_vertices = vec![vec![
        I64Vec2::new(98, 63),
        I64Vec2::new(4, 68),
        I64Vec2::new(77, 8),
        I64Vec2::new(52, 100),
        I64Vec2::new(19, 12),
    ]];

    let union_vertices = union(subj_vertices.clone(), clip_vertices.clone(), FillRule::NonZero);
    println!("{:?}", union_vertices);

    chart.draw_series(std::iter::once(Polygon::new(
        ivec_to_cartesian(union_vertices[0].clone()),
        RGBColor(255, 0, 255),
    ))).unwrap();

    root.present().unwrap();
}

#[test]
fn plot_difference() {
    let root = SVGBackend::new("test_plots/2_difference_test.svg", (500, 500)).into_drawing_area();

    root.fill(&TRANSPARENT).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .build_cartesian_2d(0.0..100.0, 0.0..100.0).unwrap();

    let subj_vertices = vec![vec![
        I64Vec2::new(100, 50),
        I64Vec2::new(10, 79),
        I64Vec2::new(65, 2),
        I64Vec2::new(65, 98),
        I64Vec2::new(10, 21),
    ]];

    let clip_vertices = vec![vec![
        I64Vec2::new(98, 63),
        I64Vec2::new(4, 68),
        I64Vec2::new(77, 8),
        I64Vec2::new(52, 100),
        I64Vec2::new(19, 12),
    ]];

    let difference_vertices = difference(subj_vertices.clone(), clip_vertices.clone(), FillRule::NonZero);
    println!("{:?}", difference_vertices);

    for path in difference_vertices {
        chart.draw_series(std::iter::once(Polygon::new(
            ivec_to_cartesian(path.clone()),
            RGBColor(255, 0, 255),
        ))).unwrap();
    }

    root.present().unwrap();
}

#[test]
fn plot_xor() {
    let root = SVGBackend::new("test_plots/3_xor_test.svg", (500, 500)).into_drawing_area();

    root.fill(&TRANSPARENT).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .build_cartesian_2d(0.0..100.0, 0.0..100.0).unwrap();

    let subj_vertices = vec![vec![
        I64Vec2::new(100, 50),
        I64Vec2::new(10, 79),
        I64Vec2::new(65, 2),
        I64Vec2::new(65, 98),
        I64Vec2::new(10, 21),
    ]];

    let clip_vertices = vec![vec![
        I64Vec2::new(98, 63),
        I64Vec2::new(4, 68),
        I64Vec2::new(77, 8),
        I64Vec2::new(52, 100),
        I64Vec2::new(19, 12),
    ]];

    let xor_vertices = xor(subj_vertices.clone(), clip_vertices.clone(), FillRule::NonZero);
    println!("{:?}", xor_vertices);

    for path in xor_vertices {
        chart.draw_series(std::iter::once(Polygon::new(
            ivec_to_cartesian(path.clone()),
            RGBColor(255, 0, 255),
        ))).unwrap();
    }

    root.present().unwrap();
}

#[test]
fn plot_offset_path() {

    // use regex to extract path string
    let svg_str = std::fs::read_to_string("Clipper2/CPP/Examples/Inflate/rabbit.svg").unwrap();
    let re = Regex::new("d=\".*\"").unwrap();

    let caps = &(re.captures(svg_str.as_str()).unwrap())[0];

    let len = caps.len();
    let path_str = &caps[3..len-1];

    // parse as points
    let (_, points) = svg_path_parser::parse(path_str).collect::<Vec<(bool, Vec<(f64, f64)>)>>().pop().unwrap();

    // convert to I64Vec2
    let mut bunny_vertices: Vec<Vec<I64Vec2>> = vec![points.into_iter().map(|x| -DVec2::from(x).as_i64vec2()).collect()];
    bunny_vertices[0] = (&bunny_vertices[0]).into_iter().map(|x| *x + I64Vec2::new(500, 520)).collect(); // translate coords

    // create inflated/offset paths
    let mut p: Vec<Vec<I64Vec2>> = bunny_vertices.clone();
    let mut solution: Vec<Vec<I64Vec2>> = vec![];
    loop {
        p = inflate_paths(p, -2.5, JoinType::Round, EndType::Polygon, 2.0, 0.0, false);
        p = simplify_paths(p, 0.025);

        if p.is_empty() {
            break;
        }

        solution.append(&mut p.clone());
    }

    // draw solution
    let root = SVGBackend::new("test_plots/4_offset_bunny.svg", (500, 500)).into_drawing_area();

    root.fill(&TRANSPARENT).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .build_cartesian_2d(0.0..500.0, 0.0..500.0).unwrap();

    let sol_iter = solution.iter().enumerate();
    for (idx, path) in sol_iter {
        chart.draw_series(std::iter::once(Polygon::new(
            ivec_to_cartesian(path.clone()),
            HSLColor(
                idx as f64 / solution.len() as f64,
                0.9,
                0.5,
            ),
        ))).unwrap();
    }

    root.present().unwrap();
}