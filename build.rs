use cmake::Config;

fn main() {

    let dst = Config::new("Clipper2/CPP")
                .define("CLIPPER2_UTILS", "OFF")
                .define("CLIPPER2_EXAMPLES", "OFF")
                .define("CLIPPER2_TESTS", "OFF")
                .define("CLIPPER2_USINGZ", "OFF")
                .define("BUILD_SHARED_LIBS", "OFF")
                .build();
    
    println!("cargo:rustc-link-search=native={}/lib", dst.display());
    println!("cargo:rustc-link-lib=static=Clipper2");
    println!("cargo:rustc-link-lib=dylib=stdc++"); // need to change this if targeting other platforms
}