use clipper2_rust::{intersect, union, difference, xor};
use clipper2_rust::FillRule;
use glam::I64Vec2;
use plotters::prelude::*;

fn ivec_to_cartesian(input: Vec<I64Vec2>) -> Vec<(f64, f64)> {
    input.into_iter().map(|x| x.as_dvec2().into()).collect()
}

#[test]
fn plot_intersection() {
    let root = BitMapBackend::new("0_intersection_test.png", (256, 256)).into_drawing_area();

    root.fill(&WHITE).unwrap();

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
    let root = BitMapBackend::new("1_union_test.png", (256, 256)).into_drawing_area();

    root.fill(&WHITE).unwrap();

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
    let root = BitMapBackend::new("2_difference_test.png", (256, 256)).into_drawing_area();

    root.fill(&WHITE).unwrap();

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
    let root = BitMapBackend::new("3_xor_test.png", (256, 256)).into_drawing_area();

    root.fill(&WHITE).unwrap();

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