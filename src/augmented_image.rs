use std::collections::HashMap;

use ffi_arcore::*;
use jni_interface;
use log;
use util;


pub fn init_augmented_image_database(ar_session: *const ArSession) -> *mut ArAugmentedImageDatabase {
    let mut ar_augmented_image_database = create_augmented_image_database(ar_session);

    let k_sample_image_name: &str = "default.jpg";
    add_asset_image_to_ar_database(ar_session, ar_augmented_image_database, k_sample_image_name);
    ar_augmented_image_database
}

pub fn create_augmented_image_database(ar_session: *const ArSession) -> *mut ArAugmentedImageDatabase {
    unsafe {
        let mut ar_augmented_image_database: *mut ArAugmentedImageDatabase = ::std::ptr::null_mut();
        ArAugmentedImageDatabase_create(ar_session, &mut ar_augmented_image_database);
        ar_augmented_image_database
    }
}

pub fn add_asset_image_to_ar_database(ar_session: *const ArSession, database: *mut ArAugmentedImageDatabase, image_name: &str) {
    unsafe {
        // load_image_from_asset
        let mut width = 0;
        let mut height = 0;
        let mut stride = 0;
        let mut index = 0;

        let mut image_pixel_buffer: *mut u8 = ::std::ptr::null_mut();

        let load_image_result = jni_interface::load_image_from_assets(
            image_name,
            &mut width,
            &mut height,
            &mut stride,
            &mut image_pixel_buffer,
        );

        if !load_image_result {
            log::e(&format!("arcore::augmented_image load image failed: {}", &load_image_result));
        } else {
            log::d(&format!("arcore::augmented_image load image width = {}, height = {}, stride = {}, image_pixel_buffer = {:?}", &width, &height, &stride, &image_pixel_buffer));

            // convert_rgba_to_grayscale
            let mut grayscale_buffer: *mut u8 = ::std::ptr::null_mut();
            util::convert_rgba_to_grayscale(image_pixel_buffer, width, height, stride, &mut grayscale_buffer);

            // add image to ArAugmentedImageDatabase
            let grayscale_stride = stride / 4;
            log::i(&format!("arcore::augmented_image grayscale_stride : {:?}", &grayscale_stride));
            let ar_status: ArStatus = ArAugmentedImageDatabase_addImage(
                ar_session,
                database,
                image_name.as_ptr(),
                grayscale_buffer as *const u8,
                width as i32,
                height as i32,
                grayscale_stride as i32,
                &mut index as *mut i32);

            if ar_status != AR_SUCCESS as i32 {
                log::e(&format!("arcore::augmented_image ArAugmentedImageDatabase_addImage failed: {}", &ar_status));
            }
        }
    }
}

pub fn track_images(session: *const ArSession, frame: *const ArFrame, map: &mut HashMap<i32, ::ColoredAnchor>) {
    log::i("arcore::augmented_image::track_images");

// Update loop, in onDraw
    unsafe {
        let mut updated_image_list: *mut ArTrackableList = ::std::ptr::null_mut();
        ArTrackableList_create(session, &mut updated_image_list);
        if updated_image_list == ::std::ptr::null_mut() {
            log::e("arcore::augmented_image::track_images updated_image_list is null");
        }
        ArFrame_getUpdatedTrackables(session, frame, AR_TRACKABLE_AUGMENTED_IMAGE as i32, updated_image_list);

        let mut image_list_size = 0;
        ArTrackableList_getSize(session,
                                updated_image_list as *const ArTrackableList,
                                &mut image_list_size as *mut i32);

        log::d(&format!("arcore::augmented_image::track_images image_list_size : {:?}", image_list_size));

        for i in 0..image_list_size {
            let mut ar_trackable: *mut ArTrackable = ::std::ptr::null_mut();
            ArTrackableList_acquireItem(session,
                                        updated_image_list as *const ArTrackableList,
                                        i,
                                        &mut ar_trackable);

            log::d(&format!("arcore::augmented_image::track_images ar_trackable : {:?}", &ar_trackable));

            let mut image: *mut ArAugmentedImage = ::std::mem::transmute::<*mut ArTrackable, *mut ArAugmentedImage>(ar_trackable);

            let mut image_tracking_state: ArTrackingState = 2;
            ArTrackable_getTrackingState(session,
                                         ar_trackable as *const ArTrackable,
                                         &mut image_tracking_state as *mut ArTrackingState);

            log::d(&format!("arcore::augmented_image::track_images image_tracking_state : {:?}", &image_tracking_state));

            if image_tracking_state == AR_TRACKING_STATE_TRACKING as i32 {
                let mut ar_pose: *mut ArPose = ::std::ptr::null_mut();
                ArPose_create(session, 0 as *const _, &mut ar_pose);
                ArAugmentedImage_getCenterPose(session, image as *const ArAugmentedImage, ar_pose);
                log::d(&format!("arcore::augmented_image::track_images ar_pose : {:?}", &ar_pose));

                let mut image_anchor: *mut ArAnchor = ::std::ptr::null_mut();
                let mut ar_status: ArStatus = ArTrackable_acquireNewAnchor(session as *mut ArSession,
                                                                           ar_trackable as *mut ArTrackable,
                                                                           ar_pose as *mut ArPose,
                                                                           &mut image_anchor as *mut *mut ArAnchor);
                if ar_status != AR_SUCCESS as i32 {
                    log::e(&format!("arcore::augmented_image::track_images ArTrackable_acquireNewAnchor ar_status: {}", &ar_status));
                    return;
                }

                let mut color = [0.0, 0.0, 0.0, 0.0];

                color[0] = 139.0;
                color[1] = 195.0;
                color[2] = 74.0;
                color[3] = 255.0;

                let colored_anchor = ::ColoredAnchor { anchor: image_anchor, color };

                log::d(&format!("arcore::augmented_image::track_images i : {}, colored_anchor : {:?}", i, &colored_anchor));

                map.insert(i, colored_anchor);
            }
        }

        ArTrackableList_destroy(updated_image_list);
        updated_image_list = ::std::ptr::null_mut();
    }
}
