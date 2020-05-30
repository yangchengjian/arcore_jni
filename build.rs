extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rustc-link-lib=dylib=arcore_sdk_c");
    println!("cargo:rustc-link-lib=dylib=jnigraphics");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        .trust_clang_mangling(false)
        .clang_arg("--sysroot=/Users/yangchengjian/captcha/sft/android-ndk-r15c/sysroot")
        .header("include/arcore_c_api.h")
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("arcore_bindings.rs"))
        .expect("Couldn't write bindings!");
}
