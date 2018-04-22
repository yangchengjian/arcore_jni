extern crate android_injected_glue;
extern crate jni;
extern crate libc;

use android_injected_glue::ffi::{jobject, JNIEnv};
use self::libc::{c_float, c_int, c_long, c_void, size_t};

#[link(name = "hello_ar_native")]
extern {
    fn create_application() -> c_long;
    fn destroy_application();
    fn on_touched(native_application: c_long, x: c_float, y: c_float);
    fn on_pause(native_application: c_long);
    fn on_resume(env: *mut ::std::os::raw::c_void, context: *mut ::std::os::raw::c_void, native_application: c_long);
    fn on_surface_created(native_application: c_long);
    fn on_display_deometry_changed(native_application: c_long, display_rotation: c_int, width: c_int, height: c_int);
    fn on_draw_frame(native_application: c_long);
    fn has_detected_planes(native_application: c_long) -> c_int;
}

pub fn ar_create_application() -> u64 {
    unsafe {
        let mut app = create_application();
        app as u64
    }
}

pub fn ar_on_resume(env: *mut JNIEnv, context: jobject, native_application: u64) {
    unsafe {
        on_resume(env as *mut ::std::os::raw::c_void, context as *mut ::std::os::raw::c_void, native_application as c_long)
    }
}

pub fn ar_on_surface_created(native_application: u64) {
    unsafe {
        on_surface_created(native_application as c_long)
    }
}

pub fn ar_on_display_geometry_changed(native_application: u64, display_rotation: u32, width: u32, height: u32) {
    unsafe {
        on_display_deometry_changed(native_application as c_long, display_rotation as c_int, width as c_int, height as c_int)
    }
}

pub fn ar_on_draw_frame(native_application: u64) {
    unsafe {
        on_draw_frame(native_application as c_long)
    }
}