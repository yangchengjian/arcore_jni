#[cfg(target_os = "android")]
extern crate android_injected_glue;
extern crate jni;
extern crate glm;
extern crate lodepng;
extern crate rgb;
extern crate sparkle;

mod ffi_arcore {
    include!(concat!(env!("OUT_DIR"), "/arcore_bindings.rs"));
}

pub mod background_renderer;
pub mod plane_renderer;
pub mod point_cloud_renderer;
pub mod util;

use std::collections::BTreeMap;
use std::rc::Rc;

use sparkle::gl;

use android_injected_glue::write_log;
use android_injected_glue::ffi::{jlong, jobject, jclass, android_app, ANativeActivity, JavaVM, JNIEnv, JNIInvokeInterface, JNINativeInterface};

use self::ffi_arcore::*;


pub const kMaxNumberOfAndroidsToRender: usize = 20;

//#[cfg(target_os = "android")]
//#[allow(non_snake_case)]
//#[no_mangle]
//pub unsafe extern fn Java_com_swarm_matrix_JniInterface_init(env: JNIEnv, _: jclass, context: jobject) -> jlong {
//    let arcore = ArCore::new(env, context);
//    &arcore  as *const _ as jlong
//}

pub fn log(message: &str) {
    write_log(message);
}

pub fn init_jni() -> (*mut JNIEnv, jobject) {
    write_log("arcore_jni::lib::init_jni");

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
    write_log("arcore_jni::lib::init_arcore");

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
    point_cloud_renderer_: Option<point_cloud_renderer::PointCloudRenderer>,
    plane_renderer_: Option<plane_renderer::PlaneRenderer>,
    tracked_obj_set_: Vec<*mut ArAnchor>,
    plane_color_map_: BTreeMap<*mut ArPlane, ::glm::Vec3>,
    pub view_mat4x4: [f32; 16],
    pub proj_mat4x4: [f32; 16],
    pub mode_mat4x4: [f32; 16],
    pub pv: Vec<f32>,
    pub uniform_mvp_mat_: gl::ffi::types::GLint,
}

impl ArCore {
    pub fn new(env: *mut JNIEnv, context: jobject) -> ArCore {
        write_log("arcore_jni::lib::ArCore::new");

        unsafe {
            let mut out_session_pointer: *mut ArSession = ::std::ptr::null_mut();
            let mut ar_status_create: ArStatus = ArSession_create(env as *mut ::std::os::raw::c_void, context as *mut ::std::os::raw::c_void, &mut out_session_pointer);
            if ar_status_create != 0 {
                write_log(&format!("arcore_jni::lib::new ArSession_create error, ar_status_create = {}", unsafe { ar_status_create }));
            }

            let mut out_config: *mut ArConfig = ::std::ptr::null_mut();
            ArConfig_create(out_session_pointer as *const ArSession, &mut out_config);

            let mut ar_status_check: ArStatus = ArSession_checkSupported(out_session_pointer as *const ArSession, out_config);
            if ar_status_check != 0 {
                write_log(&format!("arcore_jni::lib::new ArSession_checkSupported error, ar_status_check = {}", unsafe { ar_status_check }));
            }

            let mut ar_status_configure: ArStatus = ArSession_configure(out_session_pointer, out_config);
            if ar_status_configure != 0 {
                write_log(&format!("arcore_jni::lib::new ArSession_configure error, ar_status_configure = {}", unsafe { ar_status_configure }));
            }
            ArConfig_destroy(out_config);

            let mut out_frame: *mut ArFrame = ::std::ptr::null_mut();
            ArFrame_create(out_session_pointer as *const ArSession, &mut out_frame);

            ArSession_setDisplayGeometry(out_session_pointer, 0, 1, 1);

            let mut ar_status_resume: ArStatus = ArSession_resume(out_session_pointer);
            if ar_status_resume != 0 {
                write_log(&format!("arcore_jni::lib::new ArSession_resume error, ar_status_resume = {}", unsafe { ar_status_resume }));
            }


            //        cargo apk logcat | grep RustAndroidGlueStdouterr

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
                plane_renderer_: None,
                point_cloud_renderer_: None,
                tracked_obj_set_: Vec::new(),
                plane_color_map_: BTreeMap::new(),
                view_mat4x4: [0.0; 16],
                proj_mat4x4: [0.0; 16],
                mode_mat4x4: [0.0; 16],
                pv: Vec::new(),
                uniform_mvp_mat_: -1,
            }
        }
    }

    pub fn on_display_changed(&mut self, gl: &gl::Gl, display_rotation: i32, width: i32, height: i32) {
        write_log(&format!("arcore_jni::lib::on_display_changed display_rotation = {}, width = {}, height = {}", display_rotation, width, height));

        let bgr = ::background_renderer::BackgroundRenderer::initializel_content(gl);
        self.background_renderer_ = Some(bgr.clone());

        let plr = ::plane_renderer::PlaneRenderer::initializel_content(gl);
        self.plane_renderer_ = Some(plr.clone());

        let pcr = ::point_cloud_renderer::PointCloudRenderer::initializel_content(gl);
        self.point_cloud_renderer_ = Some(pcr.clone());

        //        self.on_touched((width / 2) as f32, (height / 2) as f32);

        self.display_rotation_ = display_rotation;
        self.width_ = width;
        self.height_ = height;
        if self.ar_session != ::std::ptr::null_mut() {
            unsafe { ArSession_setDisplayGeometry(self.ar_session, display_rotation, width, height) };
        }
    }

    pub fn on_draw(&mut self, gl: &gl::Gl) {
        write_log("arcore_jni::lib::on_draw");

        unsafe {
            //            gl.clear_color(0.9, 0.9, 0.9, 1.0);
            //            gl.clear(gl::DEPTH_BUFFER_BIT | gl::COLOR_BUFFER_BIT);
            //
            //            gl.enable(gl::CULL_FACE);
            //            gl.enable(gl::DEPTH_TEST);
            //            gl.enable(gl::BLEND);
            //            gl.blend_func(gl::SRC_ALPHA, gl::ONE_MINUS_SRC_ALPHA);


            write_log(&format!("arcore_jni::lib::on_draw background_renderer_ : {:?}", self.clone().background_renderer_));

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

            ArSession_setCameraTextureName(self.ar_session, texture_id);

            let mut ar_status_update: ArStatus = ArSession_update(self.ar_session, self.ar_frame);
            if ar_status_update != 0 {
                write_log(&format!("arcore_jni::lib::::on_draw ArSession_resume error, ar_status_update = {}", unsafe { ar_status_update }));
            }

            let mut out_camera: *mut ArCamera = ::std::ptr::null_mut();
            ArFrame_acquireCamera(self.ar_session, self.ar_frame, &mut out_camera);

            ArCamera_getViewMatrix(self.ar_session, out_camera, self.view_mat4x4.as_mut_ptr());

            ArCamera_getProjectionMatrix(self.ar_session,
                                         out_camera as *const ArCamera,
                                         0.1,
                                         100.0,
                                         self.proj_mat4x4.as_mut_ptr());

            write_log(&format!("arcore_jni::lib::on_draw ArCamera_getViewMatrix : {:?}", self.view_mat4x4));
            write_log(&format!("arcore_jni::lib::on_draw ArCamera_getProjectionMatrix : {:?}", self.proj_mat4x4));
            write_log(&format!("arcore_jni::lib::on_draw ArPose_getMatrix : {:?}", self.mode_mat4x4));


            let mut v = glm::mat4(self.view_mat4x4[0], self.view_mat4x4[1], self.view_mat4x4[2], self.view_mat4x4[3],
                                  self.view_mat4x4[4], self.view_mat4x4[5], self.view_mat4x4[6], self.view_mat4x4[7],
                                  self.view_mat4x4[8], self.view_mat4x4[9], self.view_mat4x4[10], self.view_mat4x4[11],
                                  self.view_mat4x4[12], self.view_mat4x4[13], self.view_mat4x4[14], self.view_mat4x4[15]);

            let mut p = glm::mat4(self.proj_mat4x4[0], self.proj_mat4x4[1], self.proj_mat4x4[2], self.proj_mat4x4[3],
                                  self.proj_mat4x4[4], self.proj_mat4x4[5], self.proj_mat4x4[6], self.proj_mat4x4[7],
                                  self.proj_mat4x4[8], self.proj_mat4x4[9], self.proj_mat4x4[10], self.proj_mat4x4[11],
                                  self.proj_mat4x4[12], self.proj_mat4x4[13], self.proj_mat4x4[14], self.proj_mat4x4[15]);

            let mut m = glm::mat4(self.mode_mat4x4[0], self.mode_mat4x4[1], self.mode_mat4x4[2], self.mode_mat4x4[3],
                                  self.mode_mat4x4[4], self.mode_mat4x4[5], self.mode_mat4x4[6], self.mode_mat4x4[7],
                                  self.mode_mat4x4[8], self.mode_mat4x4[9], self.mode_mat4x4[10], self.mode_mat4x4[11],
                                  self.mode_mat4x4[12], self.mode_mat4x4[13], self.mode_mat4x4[14], self.mode_mat4x4[15]);

            let pv = p * v;

            let pv_array_vec4 = pv.as_array();

            let mut pv_array: Vec<f32> = Vec::new();

            for i in 0..pv_array_vec4.len() {
                for j in 0..4 {
                    pv_array.push(pv_array_vec4[i][j]);
                }
            }
            self.pv = pv_array;

//            write_log(&format!("arcore_jni::lib::on_draw mvp_array : {:?}", mvp_array));


            let mut camera_tracking_state: ArTrackingState = 0;
            ArCamera_getTrackingState(self.ar_session, out_camera as *const ArCamera, &mut camera_tracking_state as *mut ArTrackingState);
            ArCamera_release(out_camera);

            bgr.draw(gl, self.ar_session as *const ArSession, self.ar_frame as *const ArFrame);
            self.background_renderer_ = Some(bgr);

            write_log(&format!("arcore_jni::lib::on_draw camera_tracking_state : {:?}", camera_tracking_state));

            if camera_tracking_state != AR_TRACKING_STATE_TRACKING as i32 {
                return;
            }

            //            // Get light estimation value.
            //            let mut ar_light_estimate: *mut ArLightEstimate = ::std::ptr::null_mut();
            //            let mut ar_light_estimate_state: ArLightEstimateState = AR_LIGHT_ESTIMATE_STATE_NOT_VALID as i32;
            //            ArLightEstimate_create(self.ar_session as *const ArSession, &mut ar_light_estimate);
            //
            //            ArFrame_getLightEstimate(self.ar_session as *const ArSession, self.ar_frame as *const ArFrame, ar_light_estimate);
            //            ArLightEstimate_getState(self.ar_session as *const ArSession,
            //                                     ar_light_estimate as *const ArLightEstimate,
            //                                     &mut ar_light_estimate_state as *mut ArLightEstimateState);
            //
            //            let mut color_correction = [1.0, 1.0, 1.0, 1.0];
            //            if ar_light_estimate_state == AR_LIGHT_ESTIMATE_STATE_VALID as i32 {
            //                ArLightEstimate_getColorCorrection(self.ar_session as *const ArSession,
            //                                                   ar_light_estimate as *const ArLightEstimate,
            //                                                   color_correction.as_mut_ptr());
            //            }
            //            ArLightEstimate_destroy(ar_light_estimate);
            //            ar_light_estimate = ::std::ptr::null_mut();


            //            // Render Andy objects.
            //            for obj_iter in self.tracked_obj_set_.clone() {
            //                let mut tracking_state: ArTrackingState = AR_TRACKING_STATE_STOPPED as i32;
            //                ArAnchor_getTrackingState(self.ar_session as *const ArSession,
            //                                          obj_iter as *const ArAnchor,
            //                                          &mut tracking_state as *mut ArTrackingState);
            //                if tracking_state == AR_TRACKING_STATE_TRACKING as i32 {
            //                    util::get_transform_matrix_from_anchor(self.ar_session, obj_iter, self.mode_mat4x4.as_mut_ptr());
            //                }
            //            }


            // Update and render planes.
            let mut plane_list: *mut ArTrackableList = ::std::ptr::null_mut();
            ArTrackableList_create(self.ar_session as *const ArSession, &mut plane_list);
            if plane_list == ::std::ptr::null_mut() {
                write_log("arcore_jni::lib::on_draw plane_list is null");
            }

            let plane_tracked_type: ArTrackableType = AR_TRACKABLE_PLANE as i32;
            ArSession_getAllTrackables(self.ar_session as *const ArSession, plane_tracked_type, plane_list);

            let mut plane_list_size = 0;
            ArTrackableList_getSize(self.ar_session as *const ArSession,
                                    plane_list as *const ArTrackableList,
                                    &mut plane_list_size as *mut i32);
            self.plane_count_ = plane_list_size;

            write_log(&format!("arcore_jni::lib::on_draw plane_list_size : {:?}", plane_list_size));

            for i in 0..plane_list_size {
                let mut ar_trackable: *mut ArTrackable = ::std::ptr::null_mut();
                ArTrackableList_acquireItem(self.ar_session as *const ArSession,
                                            plane_list as *const ArTrackableList,
                                            i,
                                            &mut ar_trackable);
                let ar_plane: *mut ArPlane = ::std::mem::transmute::<*mut ArTrackable, *mut ArPlane>(ar_trackable);
                let mut out_tracking_state: ArTrackingState = 0;
                ArTrackable_getTrackingState(self.ar_session as *const ArSession,
                                             ar_trackable as *const ArTrackable,
                                             &mut out_tracking_state as *mut ArTrackingState);
                let mut subsume_plane: *mut ArPlane = ::std::ptr::null_mut();
                ArPlane_acquireSubsumedBy(self.ar_session as *const ArSession,
                                          ar_plane as *const ArPlane,
                                          &mut subsume_plane);
                if subsume_plane != ::std::ptr::null_mut() {
                    ArTrackable_release(::std::mem::transmute::<*mut ArPlane, *mut ArTrackable>(subsume_plane));
                    continue
                }

                if out_tracking_state != AR_TRACKING_STATE_TRACKING as i32 {
                    continue
                }

                //                let plane_tracking_state: ArTrackingState;
                //                ArTrackable_getTrackingState(self.ar_session as *const ArSession,
                //                                             ::std::mem::transmute::<*mut ArPlane, *mut ArTrackable>(ar_plane) as *const ArTrackable,
                //                                             &plane_tracking_state as *mut ArTrackingState);
                //                if plane_tracking_state == AR_TRACKING_STATE_TRACKING as i32 {
                //                    let iter = match self.plane_color_map_.get(ar_plane) {
                //                        Some(x) => x;
                //                        None => ::glm::vec3(0.0, 0.0, 0.0);
                //                    };
                //                }

                let mut plr =
                    match self.clone().plane_renderer_ {
                        Some(p) => {
                            p
                        }
                        None => {
                            let plr = ::plane_renderer::PlaneRenderer::initializel_content(gl);
                            self.plane_renderer_ = Some(plr.clone());
                            plr
                        }
                    };
                plr.draw(gl, p, v, self.ar_session, ar_plane, ::glm::vec3(0.0, 0.0, 0.0));
            }


            ArTrackableList_destroy(plane_list);
            plane_list = ::std::ptr::null_mut();


            // Update and render point cloud.
            let mut ar_point_cloud: *mut ArPointCloud = ::std::ptr::null_mut();
            let point_cloud_status = ArFrame_acquirePointCloud(self.ar_session as *const ArSession,
                                                               self.ar_frame as *const ArFrame,
                                                               &mut ar_point_cloud);
            if point_cloud_status == AR_SUCCESS as i32 {
                let mut pcr =
                    match self.clone().point_cloud_renderer_ {
                        Some(p) => {
                            p
                        }
                        None => {
                            let pcr = ::point_cloud_renderer::PointCloudRenderer::initializel_content(gl);
                            self.point_cloud_renderer_ = Some(pcr.clone());
                            pcr
                        }
                    };
                pcr.draw(gl, p * v, self.ar_session, ar_point_cloud);
                ArPointCloud_release(ar_point_cloud);
            }
        }
    }

    pub fn on_touched(&mut self, x: f32, y: f32) {
        write_log("arcore_jni::ArCore::on_touched");

        unsafe {
            if self.ar_session != ::std::ptr::null_mut() && self.ar_frame != ::std::ptr::null_mut() {
                let mut hit_result_list: *mut ArHitResultList = ::std::ptr::null_mut();
                ArHitResultList_create(self.ar_session as *const ArSession, &mut hit_result_list);
                ArFrame_hitTest(self.ar_session as *const ArSession, self.ar_frame as *const ArFrame,
                                x, y,
                                hit_result_list);

                write_log(&format!("arcore_jni::lib::on_touched x = {}, y = {}", x, y));

                let mut hit_result_list_size = 0;
                ArHitResultList_getSize(self.ar_session as *const ArSession,
                                        hit_result_list as *const ArHitResultList,
                                        &mut hit_result_list_size as *mut i32);


                write_log(&format!("arcore_jni::lib::on_touched hit_result_list_size = {}", hit_result_list_size));

                let mut ar_hit_result: *mut ArHitResult = ::std::ptr::null_mut();
                for i in 0..hit_result_list_size {
                    let mut ar_hit: *mut ArHitResult = ::std::ptr::null_mut();
                    ArHitResult_create(self.ar_session as *const ArSession, &mut ar_hit);

                    ArHitResultList_getItem(self.ar_session as *const ArSession,
                                            hit_result_list as *const ArHitResultList,
                                            i,
                                            ar_hit);

                    if ar_hit == ::std::ptr::null_mut() {
                        write_log(&format!("arcore_jni::lib::on_touched ArHitResultList_getItem error"));
                        return;
                    }

                    let mut ar_trackable: *mut ArTrackable = ::std::ptr::null_mut();
                    ArHitResult_acquireTrackable(self.ar_session as *const ArSession,
                                                 ar_hit as *const ArHitResult,
                                                 &mut ar_trackable);

                    let mut ar_trackable_type: ArTrackableType = AR_TRACKABLE_NOT_VALID as i32;

                    ArTrackable_getType(self.ar_session as *const ArSession,
                                        ar_trackable as *const ArTrackable,
                                        &mut ar_trackable_type as *mut ArTrackableType);

                    if ar_trackable_type == AR_TRACKABLE_PLANE as i32 {
                        let mut ar_pose: *mut ArPose = ::std::ptr::null_mut();
                        ArPose_create(self.ar_session as *const ArSession, 0 as *const _,
                                      &mut ar_pose);

                        ArHitResult_getHitPose(self.ar_session as *const ArSession,
                                               ar_hit_result as *const ArHitResult,
                                               ar_pose);
                        let mut in_polygon = 0;
                        //                        ArPlane* ar_plane = ArAsPlane(ar_trackable);
                        let mut ar_plane: *mut ArPlane = ::std::mem::transmute::<*mut ArTrackable, *mut ArPlane>(ar_trackable);

                        ArPlane_isPoseInPolygon(self.ar_session as *const ArSession,
                                                ar_plane as *const ArPlane,
                                                ar_pose as *const ArPose,
                                                &mut in_polygon as *mut i32);

                        ArPose_destroy(ar_pose);

                        if in_polygon == 0 {
                            continue
                        }

                        ar_hit_result = ar_hit;

                        break
                    } else {
                        let mut ar_point: *mut ArPoint = ::std::mem::transmute::<*mut ArTrackable, *mut ArPoint>(ar_trackable);
                        let mut mode: ArPointOrientationMode = 0;

                        ArPoint_getOrientationMode(self.ar_session as *const ArSession,
                                                   ar_point as *const ArPoint,
                                                   &mut mode as *mut ArPointOrientationMode);

                        if mode == AR_POINT_ORIENTATION_ESTIMATED_SURFACE_NORMAL as i32 {
                            ar_hit_result = ar_hit;
                            break
                        }
                    }
                }

                if ar_hit_result != ::std::ptr::null_mut() {
                    let mut anchor: *mut ArAnchor = ::std::ptr::null_mut();
                    let ar_status: ArStatus = ArHitResult_acquireNewAnchor(self.ar_session,
                                                                           ar_hit_result,
                                                                           &mut anchor);
                    if ar_status != AR_SUCCESS as i32 {
                        write_log(&format!("arcore_jni::lib::on_touched ArHitResult_acquireNewAnchor error"));
                        return;
                    }

                    let mut tracking_state: ArTrackingState = AR_TRACKING_STATE_STOPPED as i32;
                    ArAnchor_getTrackingState(self.ar_session as *const ArSession,
                                              anchor as *const ArAnchor,
                                              &mut tracking_state as *mut ArTrackingState);

                    if tracking_state != AR_TRACKING_STATE_TRACKING as i32 {
                        ArAnchor_release(anchor);
                        return;
                    }

                    if self.tracked_obj_set_.len() >= kMaxNumberOfAndroidsToRender {
                        self.tracked_obj_set_.remove(0);
                    }

                    self.tracked_obj_set_.push(anchor);
                    ArHitResult_destroy(ar_hit_result);
                    ar_hit_result = ::std::ptr::null_mut();

                    ArHitResultList_destroy(hit_result_list);
                    hit_result_list = ::std::ptr::null_mut();
                }
            }
        }
    }

    pub fn on_finish(&self) {
        write_log("arcore_jni::lib::on_finish");

        unsafe {
            if self.ar_session != ::std::ptr::null_mut() {
                ArSession_destroy(self.ar_session);
                ArFrame_destroy(self.ar_frame);
            }
        }
    }
}