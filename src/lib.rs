#[cfg(target_os = "android")]
extern crate android_injected_glue;
extern crate glm;
extern crate lodepng;
extern crate nalgebra_glm;
extern crate ndk_sys;
extern crate rgb;
extern crate sparkle;

mod ffi_arcore {
    include!(concat!(env!("OUT_DIR"), "/arcore_bindings.rs"));
}

mod augmented_face;
mod augmented_image;
mod jni_interface;
pub mod log;
mod renderer_background;
mod renderer_plane;
mod renderer_point_cloud;
mod util;

use self::ffi_arcore::*;

use android_injected_glue::ffi::JNIEnv;
use android_injected_glue::ffi::jobject;
use std::collections::HashMap;
use sparkle::gl;

const K_MAX_NUMBER_OF_ANDROIDS_TO_RENDER: usize = 20;


pub fn init_arcore() -> ArCore {
    log::i("arcore::lib::init_arcore");

    let (env, context) = jni_interface::init_jni();
    ArCore::new(env, context)
}

#[derive(Clone, Debug)]
pub struct ColoredAnchor {
    anchor: *mut ArAnchor,
    color: [f32; 4],
}

#[derive(Clone)]
pub struct ArCore {
    ar_session: *mut ArSession,
    ar_frame: *mut ArFrame,
    loop_num: i32,

    width_: i32,
    height_: i32,
    display_rotation_: i32,

    renderer_background_: Option<renderer_background::BackgroundRenderer>,
    renderer_point_cloud_: Option<renderer_point_cloud::PointCloudRenderer>,
    renderer_plane_: Option<renderer_plane::PlaneRenderer>,

    plane_obj_map_: HashMap<i32, ColoredAnchor>,
    point_obj_map_: HashMap<i32, ColoredAnchor>,
    image_obj_map_: HashMap<i32, ColoredAnchor>,
    faces_obj_map_: HashMap<i32, ColoredAnchor>,
    number_to_render: usize,

    view_mat4x4: [f32; 16],
    proj_mat4x4: [f32; 16],
}

impl ArCore {
    pub fn new(env: *mut JNIEnv, context: jobject) -> ArCore {
        log::i("arcore::lib::new");

        unsafe {

            // Create ArSession
            let mut out_session_pointer: *mut ArSession = ::std::ptr::null_mut();
            let mut ar_status_create: ArStatus = ArSession_create(env as *mut ::std::os::raw::c_void, context as *mut ::std::os::raw::c_void, &mut out_session_pointer);
            if ar_status_create != 0 {
                log::e(&format!("arcore::lib::new ArSession_create error, ar_status_create = {}", ar_status_create));
            }

            // Create ArConfig
            let mut out_config: *mut ArConfig = ::std::ptr::null_mut();
            ArConfig_create(out_session_pointer as *const ArSession, &mut out_config);

            // Check ArSession_checkSupported
            let mut ar_status_check: ArStatus = ArSession_checkSupported(out_session_pointer as *const ArSession, out_config);
            if ar_status_check != 0 {
                log::e(&format!("arcore::lib::new ArSession_checkSupported error, ar_status_check = {}", ar_status_check));
            }

            // Create Augmented Image Database
            let mut ar_augmented_image_database: *mut ArAugmentedImageDatabase = ::augmented_image::init_augmented_image_database(out_session_pointer as *const ArSession);
            ArConfig_setAugmentedImageDatabase(out_session_pointer as *const ArSession, out_config, ar_augmented_image_database);
            ArAugmentedImageDatabase_destroy(ar_augmented_image_database);

            // Check ArSession_configure
            let mut ar_status_configure: ArStatus = ArSession_configure(out_session_pointer, out_config);
            if ar_status_configure != 0 {
                log::e(&format!("arcore::lib::new ArSession_configure error, ar_status_configure = {}", ar_status_configure));
            }
            ArConfig_destroy(out_config);

            // Create ArFrame
            let mut out_frame: *mut ArFrame = ::std::ptr::null_mut();
            ArFrame_create(out_session_pointer as *const ArSession, &mut out_frame);

            ArSession_setDisplayGeometry(out_session_pointer, 0, 1, 1);

            let mut ar_status_resume: ArStatus = ArSession_resume(out_session_pointer);
            if ar_status_resume != 0 {
                log::e(&format!("arcore::lib::new ArSession_resume error, ar_status_resume = {}", ar_status_resume));
            }


            ArCore {
                ar_session: out_session_pointer,
                ar_frame: out_frame,
                loop_num: 0,

                width_: 1,
                height_: 1,
                display_rotation_: 0,

                renderer_background_: None,
                renderer_plane_: None,
                renderer_point_cloud_: None,

                plane_obj_map_: HashMap::new(),
                point_obj_map_: HashMap::new(),
                image_obj_map_: HashMap::new(),
                faces_obj_map_: HashMap::new(),
                number_to_render: 0,

                view_mat4x4: [0.0; 16],
                proj_mat4x4: [0.0; 16],
            }
        }
    }

    pub fn get_proj_matrix(&self) -> [f32; 16] {
        log::print_matrix("arcore::lib::::get_proj_matrix", &self.proj_mat4x4);
        self.proj_mat4x4
    }

    pub fn get_view_matrix(&self) -> [f32; 16] {
        log::print_matrix("arcore::lib::::get_view_matrix", &self.view_mat4x4);
        self.view_mat4x4
    }

    // 1: plane, 2: point, 3: images, 4: faces
    pub fn get_mode_matrix(&self, track_type: i32, index: i32) -> [f32; 16] {
        log::i(&format!("arcore::lib::get_mode_matrix track_type = {}, index = {}", track_type, index));
        match track_type {
            1 => get_matrix_by_anchor_and_index(self.ar_session, &self.plane_obj_map_, index),
            2 => get_matrix_by_anchor_and_index(self.ar_session, &self.point_obj_map_, index),
            3 => get_matrix_by_anchor_and_index(self.ar_session, &self.image_obj_map_, index),
            4 => get_matrix_by_anchor_and_index(self.ar_session, &self.faces_obj_map_, index),
            _ => get_matrix_by_anchor_and_index(self.ar_session, &self.plane_obj_map_, index),
        }
    }

    pub fn get_view_mode_matrix(&self, track_type: i32, index: i32) -> [f32; 16] {
        log::i(&format!("arcore::lib::get_view_mode_matrix track_type = {}, index = {}", track_type, index));
        let mode_mat4x4 = self.get_mode_matrix(track_type, index);
        let vm = util::get_mat4_from_array(self.get_view_matrix()) * util::get_mat4_from_array(mode_mat4x4);
        log::print_matrix("arcore::lib::get_view_mode_matrix", &util::get_array_from_mat4(vm));
        util::get_array_from_mat4(vm)
    }

    pub fn get_light_estimation(&mut self) -> [f32; 4] {
        log::i("arcore::lib::light_estimation");
        // Get light estimation value.
        unsafe {
            let mut ar_light_estimate: *mut ArLightEstimate = ::std::ptr::null_mut();
            ArLightEstimate_create(self.ar_session as *const ArSession, &mut ar_light_estimate);
            ArFrame_getLightEstimate(self.ar_session as *const ArSession, self.ar_frame as *const ArFrame, ar_light_estimate);

            let mut ar_light_estimate_state: ArLightEstimateState = AR_LIGHT_ESTIMATE_STATE_NOT_VALID as i32;
            ArLightEstimate_getState(self.ar_session as *const ArSession,
                                     ar_light_estimate as *const ArLightEstimate,
                                     &mut ar_light_estimate_state as *mut ArLightEstimateState);

            let mut color_correction = [1.0, 1.0, 1.0, 1.0];
            if ar_light_estimate_state == AR_LIGHT_ESTIMATE_STATE_VALID as i32 {
                ArLightEstimate_getColorCorrection(self.ar_session as *const ArSession,
                                                   ar_light_estimate as *const ArLightEstimate,
                                                   color_correction.as_mut_ptr());
            }
            ArLightEstimate_destroy(ar_light_estimate);
            ar_light_estimate = ::std::ptr::null_mut();

            color_correction
        }
    }

    pub fn on_display_changed(&mut self, gl: &gl::Gl, display_rotation: i32, width: i32, height: i32) {
        log::i(&format!("arcore::lib::on_display_changed display_rotation = {}, width = {}, height = {}", display_rotation, width, height));

        self.init_renderers(gl);

        self.display_rotation_ = display_rotation;
        self.width_ = width;
        self.height_ = height;
        if self.ar_session != ::std::ptr::null_mut() {
            unsafe { ArSession_setDisplayGeometry(self.ar_session, display_rotation, width, height) };
        }
    }

    pub fn on_draw(&mut self, gl: &gl::Gl) {
        log::i("arcore::lib::on_draw");

        unsafe {
            let mut texture_id = self.clone().renderer_background_.unwrap().get_texture_id();
            ArSession_setCameraTextureName(self.ar_session, texture_id);

            let mut ar_status_update: ArStatus = ArSession_update(self.ar_session, self.ar_frame);

            let mut out_camera: *mut ArCamera = ::std::ptr::null_mut();
            ArFrame_acquireCamera(self.ar_session as *const ArSession, self.ar_frame as *const ArFrame, &mut out_camera);

            let mut camera_tracking_state: ArTrackingState = 0;
            ArCamera_getTrackingState(self.ar_session as *const ArSession, out_camera as *const ArCamera, &mut camera_tracking_state as *mut ArTrackingState);
            ArCamera_release(out_camera);

            if ar_status_update != 0 {
                log::e(&format!("arcore::lib::::on_draw ArSession_resume error, ar_status_update = {}", ar_status_update));
            } else if camera_tracking_state != AR_TRACKING_STATE_TRACKING as i32 {
                log::e(&format!("arcore::lib::::on_draw ArCamera_getTrackingState error, camera_tracking_state = {}", camera_tracking_state));
            } else {
                self.update_proj_view_matrix();

                let p = util::get_mat4_from_array(self.proj_mat4x4);
                let v = util::get_mat4_from_array(self.view_mat4x4);

                self.render_background(gl);

                self.render_point_cloud(gl, p * v);

//                self.render_planes(gl);

                ::augmented_image::track_images(self.ar_session as *const ArSession, self.ar_frame as *const ArFrame, &mut self.image_obj_map_);

                ::augmented_face::track_faces(self.ar_session as *const ArSession, self.ar_frame as *const ArFrame, &mut self.faces_obj_map_);

                if self.loop_num % 10 == 0 {
                    log::d(&format!("arcore::lib::::on_draw loop_num = {}", &self.loop_num));
                }
            }
        }
        self.loop_num += 1;
    }

    pub fn on_touched(&mut self, x: f32, y: f32) -> i32 {
        log::i("arcore::lib::on_touched");

        // Update loop if not hited, in onDraw
        unsafe {
            if self.ar_session != ::std::ptr::null_mut() && self.ar_frame != ::std::ptr::null_mut() {
                let mut hit_result_list: *mut ArHitResultList = ::std::ptr::null_mut();
                ArHitResultList_create(self.ar_session as *const ArSession, &mut hit_result_list);
                ArFrame_hitTest(self.ar_session as *const ArSession, self.ar_frame as *const ArFrame,
                                x, y,
                                hit_result_list);

                log::d(&format!("arcore::lib::on_touched x = {}, y = {}", x, y));

                let mut hit_result_list_size = 0;
                ArHitResultList_getSize(self.ar_session as *const ArSession,
                                        hit_result_list as *const ArHitResultList,
                                        &mut hit_result_list_size as *mut i32);

                log::d(&format!("arcore::lib::on_touched hit_result_list_size = {}", hit_result_list_size));
                let mut trackable_type: ArTrackableType = AR_TRACKABLE_NOT_VALID as i32;


                for i in 0..hit_result_list_size {
                    let mut ar_hit_result: *mut ArHitResult = ::std::ptr::null_mut();
                    ArHitResult_create(self.ar_session as *const ArSession, &mut ar_hit_result);
                    ArHitResultList_getItem(self.ar_session as *const ArSession,
                                            hit_result_list as *const ArHitResultList,
                                            i,
                                            ar_hit_result);
                    if ar_hit_result == ::std::ptr::null_mut() {
                        log::e(&format!("arcore::lib::on_touched ArHitResultList_getItem error"));
                        return -1;
                    }


                    let mut ar_trackable: *mut ArTrackable = ::std::ptr::null_mut();
                    ArHitResult_acquireTrackable(self.ar_session as *const ArSession,
                                                 ar_hit_result as *const ArHitResult,
                                                 &mut ar_trackable);

                    let mut ar_trackable_type: ArTrackableType = AR_TRACKABLE_NOT_VALID as i32;
                    ArTrackable_getType(self.ar_session as *const ArSession,
                                        ar_trackable as *const ArTrackable,
                                        &mut ar_trackable_type as *mut ArTrackableType);

                    if ar_trackable_type == AR_TRACKABLE_PLANE as i32 {
                        let mut ar_pose: *mut ArPose = ::std::ptr::null_mut();
                        ArPose_create(self.ar_session as *const ArSession, 0 as *const _, &mut ar_pose);

                        ArHitResult_getHitPose(self.ar_session as *const ArSession,
                                               ar_hit_result as *const ArHitResult,
                                               ar_pose);

                        let mut in_polygon = 0;
                        let mut ar_plane: *mut ArPlane = ::std::mem::transmute::<*mut ArTrackable, *mut ArPlane>(ar_trackable);

                        ArPlane_isPoseInPolygon(self.ar_session as *const ArSession,
                                                ar_plane as *const ArPlane,
                                                ar_pose as *const ArPose,
                                                &mut in_polygon as *mut i32);

                        ArPose_destroy(ar_pose);

                        if in_polygon != 0 {
                            log::e(&format!("arcore::lib::on_touched in_polygon = {}", in_polygon));
                            return -1;
                        }

                        let mut anchor: *mut ArAnchor = ::std::ptr::null_mut();
                        let ar_status: ArStatus = ArHitResult_acquireNewAnchor(self.ar_session,
                                                                               ar_hit_result,
                                                                               &mut anchor);
                        if ar_status != AR_SUCCESS as i32 {
                            log::e(&format!("arcore::lib::on_touched ArHitResult_acquireNewAnchor error"));
                            return -1;
                        }

                        let mut tracking_state: ArTrackingState = AR_TRACKING_STATE_STOPPED as i32;
                        ArAnchor_getTrackingState(self.ar_session as *const ArSession,
                                                  anchor as *const ArAnchor,
                                                  &mut tracking_state as *mut ArTrackingState);

                        if tracking_state != AR_TRACKING_STATE_TRACKING as i32 {
                            log::e(&format!("arcore::lib::on_touched tracking_state = {}", tracking_state));
                            ArAnchor_release(anchor);
                            return -1;
                        }


                        let mut color = [0.0, 0.0, 0.0, 0.0];
                        match trackable_type as u32 {
                            AR_TRACKABLE_POINT => {
                                color[0] = 66.0;
                                color[1] = 133.0;
                                color[2] = 244.0;
                                color[3] = 255.0;
                            }
                            AR_TRACKABLE_PLANE => {
                                color[0] = 139.0;
                                color[1] = 195.0;
                                color[2] = 74.0;
                                color[3] = 255.0;
                            }
                            _ => {
                                color[0] = 0.0;
                                color[1] = 0.0;
                                color[2] = 0.0;
                                color[3] = 0.0;
                            }
                        }

                        let colored_anchor = ColoredAnchor { anchor, color };

                        log::d(&format!("arcore::lib::on_touched i : {}, colored_anchor: {:?}", i, &colored_anchor));

                        self.plane_obj_map_.insert(i, colored_anchor);

                        return i;
                    } else if AR_TRACKABLE_POINT as i32 == ar_trackable_type {
//                        let mut ar_point: *mut ArPoint = ::std::mem::transmute::<*mut ArTrackable, *mut ArPoint>(ar_trackable);
//                        let mut mode: ArPointOrientationMode = 0;
//
//                        ArPoint_getOrientationMode(self.ar_session as *const ArSession,
//                                                   ar_point as *const ArPoint,
//                                                   &mut mode as *mut ArPointOrientationMode);
//
//                        if mode == AR_POINT_ORIENTATION_ESTIMATED_SURFACE_NORMAL as i32 {
//                            ar_hit_result = ar_hit;
//                            trackable_type = ar_trackable_type;
//                            break;
//                        }
                        return -1;
                    }
                }

                ArHitResultList_destroy(hit_result_list);
                hit_result_list = ::std::ptr::null_mut();

                return -1;
            }
            -1
        }
    }

    pub fn on_pause(&self) {
        log::i("arcore::lib::on_pause");

        unsafe {
            if self.ar_session != ::std::ptr::null_mut() {
                ArSession_pause(self.ar_session);
            }
        }
    }

    pub fn on_finish(&self) {
        log::i("arcore::lib::on_finish");

        unsafe {
            if self.ar_session != ::std::ptr::null_mut() {
                ArSession_destroy(self.ar_session);
                ArFrame_destroy(self.ar_frame);
            }
        }
    }


    fn init_renderers(&mut self, gl: &gl::Gl) {
        log::i("arcore::lib::init_renderers");

        let bgr = ::renderer_background::BackgroundRenderer::new(gl);
        self.renderer_background_ = Some(bgr);

        let plr = ::renderer_plane::PlaneRenderer::new(gl);
        self.renderer_plane_ = Some(plr);

        let pcr = ::renderer_point_cloud::PointCloudRenderer::new(gl);
        self.renderer_point_cloud_ = Some(pcr);
    }

    fn render_background(&mut self, gl: &gl::Gl) {
        self.clone().renderer_background_.unwrap().draw(gl, self.ar_session as *const ArSession, self.ar_frame as *const ArFrame);
    }

    fn render_point_cloud(&mut self, gl: &gl::Gl, mvp_matrix: ::glm::Mat4) {
        // Update and render point cloud.
        unsafe {
            let mut ar_point_cloud: *mut ArPointCloud = ::std::ptr::null_mut();
            let point_cloud_status =
                ArFrame_acquirePointCloud(self.ar_session as *const ArSession,
                                          self.ar_frame as *const ArFrame,
                                          &mut ar_point_cloud);

            if point_cloud_status == AR_SUCCESS as i32 {
                self.clone().renderer_point_cloud_.unwrap().draw(gl, mvp_matrix, self.ar_session, ar_point_cloud);
                ArPointCloud_release(ar_point_cloud);
            }
        }
    }

    fn render_planes(&mut self, gl: &gl::Gl) {
        // Update loop, in onDraw
        unsafe {

            // Update and render planes.
            let mut plane_list: *mut ArTrackableList = ::std::ptr::null_mut();
            ArTrackableList_create(self.ar_session as *const ArSession, &mut plane_list);
            if plane_list == ::std::ptr::null_mut() {
                log::e("arcore::lib::render_planes plane_list is null");
            }

            let plane_tracked_type: ArTrackableType = AR_TRACKABLE_PLANE as i32;
            ArSession_getAllTrackables(self.ar_session as *const ArSession, plane_tracked_type, plane_list);

            let mut plane_list_size = 0;
            ArTrackableList_getSize(self.ar_session as *const ArSession,
                                    plane_list as *const ArTrackableList,
                                    &mut plane_list_size as *mut i32);

            log::i(&format!("arcore::lib::render_planes plane_list_size : {:?}", plane_list_size));

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
                    continue;
                }

                if out_tracking_state != AR_TRACKING_STATE_TRACKING as i32 {
                    continue;
                }

                let p = util::get_mat4_from_array(self.proj_mat4x4);
                let v = util::get_mat4_from_array(self.view_mat4x4);

                self.clone().renderer_plane_.unwrap().draw(gl, p, v, self.ar_session, ar_plane, ::glm::vec3(255.0, 255.0, 255.0));
            }

            ArTrackableList_destroy(plane_list);
            plane_list = ::std::ptr::null_mut();
        }
    }

    fn update_proj_view_matrix(&mut self) {
        unsafe {
            let mut out_camera: *mut ArCamera = ::std::ptr::null_mut();
            ArFrame_acquireCamera(self.ar_session as *const ArSession,
                                  self.ar_frame as *const ArFrame,
                                  &mut out_camera);

            ArCamera_getProjectionMatrix(self.ar_session as *const ArSession,
                                         out_camera as *const ArCamera,
                                         0.1,
                                         100.0,
                                         self.proj_mat4x4.as_mut_ptr());
            ArCamera_getViewMatrix(self.ar_session as *const ArSession,
                                   out_camera as *const ArCamera,
                                   self.view_mat4x4.as_mut_ptr());

            ArCamera_release(out_camera);
        }
    }
}

fn get_matrix_by_anchor_and_index(session: *mut ArSession, obj_map: &HashMap<i32, ColoredAnchor>, index: i32) -> [f32; 16] {
    let mut mode_mat4x4 = [0.0; 16];
    match obj_map.get(&index) {
        Some(colored_anchor) => {
            let mut tracking_state: ArTrackingState = AR_TRACKING_STATE_STOPPED as i32;
            unsafe { ArAnchor_getTrackingState(session as *const ArSession, colored_anchor.anchor, &mut tracking_state as *mut ArTrackingState) };
            if tracking_state == AR_TRACKING_STATE_TRACKING as i32 {
                get_transform_matrix_from_anchor(session, colored_anchor.anchor, mode_mat4x4.as_mut_ptr())
            }
        }
        None => {}
    }
    log::print_matrix("arcore::lib::::get_mode_matrix", &mode_mat4x4);
    mode_mat4x4
}

fn get_transform_matrix_from_anchor(session: *mut ArSession, anchor: *mut ArAnchor, out_model_mat: *mut f32) {
    unsafe {
        if out_model_mat == ::std::ptr::null_mut() {
            return;
        }
        let mut out_pose: *mut ArPose = 0 as *mut _;
        ArPose_create(session as *const ArSession, 0 as *const _, &mut out_pose);
        ArAnchor_getPose(session as *const ArSession, anchor as *const ArAnchor, out_pose as *mut ArPose);
        ArPose_getMatrix(session as *const ArSession, out_pose as *const ArPose, out_model_mat);
    }
}

fn get_ar_image_from_camera(session: *const ArSession, frame: *const ArFrame) -> *mut ArImage {
    unsafe {
        let mut ar_image: *mut ArImage = ::std::ptr::null_mut();
        let ar_status: ArStatus = ArFrame_acquireCameraImage(session as *mut ArSession, frame as *mut ArFrame, &mut ar_image);
        if ar_status != AR_SUCCESS as i32 {
            log::d(&format!("arcore::lib::get_ar_image_from_camera ar_status : {:?}", &ar_status));
        }
        ar_image
    }
}

fn get_a_image(ar_image: *const ArImage) -> *const ::ndk_sys::AImage {
    unsafe {
        let mut out_ndk_image: *const AImage = ::std::ptr::null_mut();
        ArImage_getNdkImage(ar_image, &mut out_ndk_image as *mut *const AImage);
        let a_image = ::std::mem::transmute::<*const AImage, *mut ::ndk_sys::AImage>(out_ndk_image); // transfor arcore AImage to NDK AImage
        a_image
    }
}

// java.lang.UnsatisfiedLinkError: dlopen failed: cannot locate symbol "ArImage_getFormat"
//fn get_ar_image_format(session: *const ArSession, ar_image: *const ArImage) -> i32 {
//    unsafe {
//        let mut height = 0;
//        ArImage_getFormat(session, ar_image, &mut height as *mut i32);
//        height
//    }
//}

fn get_ar_image_height(session: *const ArSession, ar_image: *const ArImage) -> i32 {
    unsafe {
        let mut height = 0;
        ArImage_getHeight(session, ar_image, &mut height as *mut i32);
        height
    }
}

fn get_ar_image_width(session: *const ArSession, ar_image: *const ArImage) -> i32 {
    unsafe {
        let mut width = 0;
        ArImage_getWidth(session, ar_image, &mut width as *mut i32);
        width
    }
}

fn get_ar_image_timestamp(session: *const ArSession, ar_image: *const ArImage) -> i64 {
    unsafe {
        let mut timestamp = 0;
        ArImage_getTimestamp(session, ar_image, &mut timestamp as *mut i64);
        timestamp
    }
}