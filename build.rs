use cmake::Config;

fn main() {

    // build Clipper2
    let dst = Config::new("Clipper2/CPP")
                .define("CLIPPER2_UTILS", "OFF")
                .define("CLIPPER2_EXAMPLES", "OFF")
                .define("CLIPPER2_TESTS", "OFF")
                .define("CLIPPER2_USINGZ", "OFF")
                .define("BUILD_SHARED_LIBS", "OFF")
                .build();
    
    // build wrapper
    cc::Build::new()
        .file("src/wrapper.cpp")
        .include("Clipper2/CPP/Clipper2Lib/include")
        .compile("clipper2wrap");

    println!("cargo:rustc-link-search=native={}/lib", dst.display());
    println!("cargo:rustc-link-lib=static=Clipper2");
    println!("cargo:rustc-link-lib=static=clipper2wrap");
    println!("cargo:rustc-link-lib=dylib=stdc++"); // need to change this if targeting other platforms
}