use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use android_injected_glue::ffi::__android_log_write;
use android_injected_glue::ffi::ANDROID_LOG_DEBUG;
use android_injected_glue::ffi::ANDROID_LOG_ERROR;
use android_injected_glue::ffi::ANDROID_LOG_INFO;

use android_injected_glue::write_log;

const TAG: &str = "RustAndroidGlueStdouterr";

fn get_c_char(s: &str) -> *const c_char {
    let c_str = CString::new(s).unwrap();
    c_str.as_ptr()
}

/// Log Debug
pub fn d(message: &str) {
    log(ANDROID_LOG_DEBUG, message);
}

/// Log Info
pub fn i(message: &str) {
    log(ANDROID_LOG_INFO, message);
}

/// Log Error
pub fn e(message: &str) {
    log(ANDROID_LOG_ERROR, message);
}

pub fn log(level: i32, message: &str) {
    write_log(message);
//    unsafe {
//        __android_log_write(level, get_c_char(TAG), get_c_char(message));
//    }
}

pub fn print_matrix(tag: &str, mat: &[f32]) {
    d(&format!("{} : {:?}", tag, mat));
}