#[cfg(target_os = "android")]
extern crate android_injected_glue;
extern crate jni;
extern crate glm;
extern crate gleam;

mod ffi_arcore {
    include!(concat!(env!("OUT_DIR"), "/arcore_bindings.rs"));
}

pub mod background_renderer;
pub mod ffi;
pub mod util;

use std::rc::Rc;

use gleam::gl;

use android_injected_glue::write_log;
use android_injected_glue::ffi::{jobject, jclass, android_app, ANativeActivity, JavaVM, JNIEnv, JNIInvokeInterface, JNINativeInterface};

use glm::mat4;
use self::ffi_arcore::*;

#[cfg(target_os = "android")]
#[allow(non_snake_case)]
#[no_mangle]
pub unsafe extern fn Java_com_swarm_matrix_JniInterface_init(env: JNIEnv, _: jclass, context: jobject) {
    //    let ar_thread = ArThread::new(env, context);
    //    entry_servo(&ar_thread);
}

pub fn init_jni() -> (*mut JNIEnv, jobject) {
    write_log("servo_jni::init_jni");

    use android_injected_glue::Event;
    use android_injected_glue::{add_sender, get_app, write_log};
    use android_injected_glue::ffi::{jobject, android_app, ANativeActivity, JavaVM, JNIEnv, JNIInvokeInterface, JNINativeInterface};

    unsafe {
        let app: &mut android_app = get_app();
        let activity: *const ANativeActivity = (*app).activity;

        let vm: *mut JavaVM = unsafe { (*activity).vm };
        let mut env: *mut JNIEnv = unsafe { (*activity).env };
        let obj: jobject = unsafe { (*activity).clazz };
        //    pub AttachCurrentThread:               extern fn(*mut JavaVM, *mut *mut JNIEnv, *mut c_void) -> jint,

        let jni_invoke_interface: *const JNIInvokeInterface = unsafe { (*vm).functions };
        let jni_native_interface: *const JNINativeInterface = unsafe { (*env).functions };

        let attach = unsafe { (*jni_invoke_interface).AttachCurrentThread };
        attach(vm, &mut env, std::ptr::null_mut());

        //    write_log(&format!("_______status_attach_________{}", status_attach));

        //        let get_object_class = unsafe { (*jni_native_interface).GetObjectClass };
        //        let find_class = unsafe { (*jni_native_interface).FindClass };
        //        let get_method_id = unsafe { (*jni_native_interface).GetMethodID };
        //        let call_object = unsafe { (*jni_native_interface).CallObjectMethod };
        //        let call_int = unsafe { (*jni_native_interface).CallIntMethod };

        (env, obj as jobject)
    }
}

pub fn init_arcore() -> ArCore {
    write_log("servo_jni::init_ar");

    let (env, context) = init_jni();
    ArCore::new(env, context)
}

#[derive(Clone)]
pub struct ArCore {
    pub ar_session: *mut ArSession,
    pub ar_frame: *mut ArFrame,
    install_requested_: bool,
    width_: i32,
    height_: i32,
    display_rotation_: i32,

    first_plane_has_been_found_: bool,
    plane_count_: i32,
    background_renderer_: Option<background_renderer::BackgroundRenderer>,
    //    tracked_obj_set_:Vec<ArAnchor>,

    //    point_cloud_renderer_,
    //    plane_renderer_,
    pub view_mat4x4: [f32; 16],
    pub proj_mat4x4: [f32; 16],
    pub uniform_mvp_mat_: Option<gl::GLuint>,
}

impl ArCore {
    pub fn new(env: *mut JNIEnv, context: jobject) -> ArCore {
        write_log("servo_jni::Ar::new");

        unsafe {
            let mut out_session_pointer: *mut ArSession = ::std::ptr::null_mut();
            let mut ar_status_create: ArStatus = ArSession_create(env as *mut ::std::os::raw::c_void, context as *mut ::std::os::raw::c_void, &mut out_session_pointer);
            if ar_status_create != 0 {
                write_log(&format!("ArSession_create error, ar_status_create = {}", unsafe { ar_status_create }));
            }

            let mut out_config: *mut ArConfig = ::std::ptr::null_mut();
            ArConfig_create(out_session_pointer as *const ArSession, &mut out_config);

            let mut ar_status_check: ArStatus = ArSession_checkSupported(out_session_pointer as *const ArSession, out_config);
            if ar_status_check != 0 {
                write_log(&format!("ArSession_checkSupported error, ar_status_check = {}", unsafe { ar_status_check }));
            }

            let mut ar_status_configure: ArStatus = ArSession_configure(out_session_pointer, out_config);
            if ar_status_configure != 0 {
                write_log(&format!("ArSession_configure error, ar_status_configure = {}", unsafe { ar_status_configure }));
            }
            ArConfig_destroy(out_config);

            let mut out_frame: *mut ArFrame = ::std::ptr::null_mut();
            ArFrame_create(out_session_pointer as *const ArSession, &mut out_frame);

            ArSession_setDisplayGeometry(out_session_pointer, 0, 1, 1);

            let mut ar_status_resume: ArStatus = ArSession_resume(out_session_pointer);
            if ar_status_resume != 0 {
                write_log(&format!("ArSession_resume error, ar_status_resume = {}", unsafe { ar_status_resume }));
            }


            //        cargo apk logcat | grep RustAndroidGlueStdouterr

            let view_mat4x4 = [0.0; 16];
            let proj_mat4x4 = [0.0; 16];

            ArCore {
                ar_session: out_session_pointer,
                ar_frame: out_frame,
                install_requested_: false,
                width_: 1,
                height_: 1,
                display_rotation_: 0,
                first_plane_has_been_found_: false,
                plane_count_: 0,
                background_renderer_: None,
                view_mat4x4: view_mat4x4,
                proj_mat4x4: proj_mat4x4,
                uniform_mvp_mat_: None,
            }
        }
    }

    pub fn on_display_changed(&mut self, gl: &gl::Gl, display_rotation: i32, width: i32, height: i32) {
        write_log(&format!("servo_jni::Ar::on_display_changed display_rotation = {}, width = {}, height = {}", display_rotation, width, height));

        let bgr = ::background_renderer::BackgroundRenderer::initializel_content(gl);
        self.background_renderer_ = Some(bgr.clone());

        self.display_rotation_ = display_rotation;
        self.width_ = width;
        self.height_ = height;
        if self.ar_session != ::std::ptr::null_mut() {
            unsafe { ArSession_setDisplayGeometry(self.ar_session, display_rotation, width, height) };
        }
    }

    pub fn on_draw(&mut self, gl: &gl::Gl) {
        write_log("servo_jni::Ar::on_draw");

        unsafe {
            //            gl.clear_color(0.9, 0.9, 0.9, 1.0);
            //            gl.clear(gl::DEPTH_BUFFER_BIT | gl::COLOR_BUFFER_BIT);
            //
            //            gl.enable(gl::CULL_FACE);
            //            gl.enable(gl::DEPTH_TEST);
            //            gl.enable(gl::BLEND);
            //            gl.blend_func(gl::SRC_ALPHA, gl::ONE_MINUS_SRC_ALPHA);


            write_log(&format!("background_renderer_ : {:?}", self.clone().background_renderer_));

            let mut bgr =
                match self.clone().background_renderer_ {
                    Some(b) => {
                        b
                    }
                    None => {
                        let bgr = ::background_renderer::BackgroundRenderer::initializel_content(gl);
                        self.background_renderer_ = Some(bgr.clone());
                        bgr
                    }
                };


            let mut texture_id = bgr.get_texture_id();

            write_log(&format!("get_texture_id : {:?}", texture_id));

            ArSession_setCameraTextureName(self.ar_session, texture_id);
            let mut ar_status_update: ArStatus = ArSession_update(self.ar_session, self.ar_frame);
            if ar_status_update != 0 {
                write_log(&format!("ArSession_resume error, ar_status_update = {}", unsafe { ar_status_update }));
            }

            let mut out_camera: *mut ArCamera = ::std::ptr::null_mut();
            ArFrame_acquireCamera(self.ar_session, self.ar_frame, &mut out_camera);

            write_log(&format!("out_camera : {:?}", unsafe { *out_camera }));

            ArCamera_getViewMatrix(self.ar_session, out_camera, self.view_mat4x4.as_mut_ptr());

            ArCamera_getProjectionMatrix(self.ar_session,
                                         out_camera as *const ArCamera,
                                         0.1,
                                         100.0,
                                         self.proj_mat4x4.as_mut_ptr());

            write_log(&format!("ArCamera_getViewMatrix : {:?}", self.view_mat4x4));
            write_log(&format!("ArCamera_getProjectionMatrix : {:?}", self.proj_mat4x4));

            let mut x = 0;
            let out_tracking_state: *mut ArTrackingState = &mut x;
            ArCamera_getTrackingState(self.ar_session, out_camera as *const ArCamera, out_tracking_state);
            ArCamera_release(out_camera);

            bgr.draw(gl, self.ar_session, self.ar_frame);
            self.background_renderer_ = Some(bgr);

            //            write_log(&format!("out_tracking_state : {:?}", *out_tracking_state));
        }
    }

    pub fn on_finish(&self) {
        unsafe {
            if self.ar_session != ::std::ptr::null_mut() {
                ArSession_destroy(self.ar_session);
                ArFrame_destroy(self.ar_frame);
            }
        }
    }
}