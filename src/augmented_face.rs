use std::collections::HashMap;

use ffi_arcore::*;
use jni_interface;
use log;
use util;

pub fn track_faces(session: *const ArSession, frame: *const ArFrame, map: &mut HashMap<i32, ::ColoredAnchor>) {
    log::i("arcore::augmented_face::track_faces");

    // Update loop, in onDraw
    unsafe {
        let mut ar_trackable_list: *mut ArTrackableList = ::std::ptr::null_mut();
        ArTrackableList_create(session, &mut ar_trackable_list);
        if ar_trackable_list == ::std::ptr::null_mut() {
            log::e("arcore::augmented_face::track_faces ar_trackable_list is null");
        }
        ArFrame_getUpdatedTrackables(session, frame, AR_TRACKABLE_FACE as i32, ar_trackable_list);

        let mut ar_trackable_list_size = 0;
        ArTrackableList_getSize(session,
                                ar_trackable_list as *const ArTrackableList,
                                &mut ar_trackable_list_size as *mut i32);

        log::d(&format!("arcore::augmented_face::track_faces ar_trackable_list_size : {:?}", ar_trackable_list_size));

        for i in 0..ar_trackable_list_size {
            let mut ar_trackable: *mut ArTrackable = ::std::ptr::null_mut();
            ArTrackableList_acquireItem(session,
                                        ar_trackable_list as *const ArTrackableList,
                                        i,
                                        &mut ar_trackable);

            log::d(&format!("arcore::augmented_face::track_faces ar_trackable : {:?}", &ar_trackable));

            let mut face: *mut ArAugmentedFace = ::std::mem::transmute::<*mut ArTrackable, *mut ArAugmentedFace>(ar_trackable);

            let mut tracking_state: ArTrackingState = 2;
            ArTrackable_getTrackingState(session,
                                         ar_trackable as *const ArTrackable,
                                         &mut tracking_state as *mut ArTrackingState);

            log::d(&format!("arcore::augmented_face::track_faces tracking_state : {:?}", &tracking_state));

            if tracking_state == AR_TRACKING_STATE_TRACKING as i32 {

                let mut ar_pose: *mut ArPose = ::std::ptr::null_mut();
                ArPose_create(session, 0 as *const _, &mut ar_pose);
                ArAugmentedFace_getCenterPose(session, face as *const ArAugmentedFace, ar_pose);
                log::d(&format!("arcore::augmented_face::track_faces ar_pose : {:?}", &ar_pose));

                let mut ar_anchor: *mut ArAnchor = ::std::ptr::null_mut();
                let mut ar_status: ArStatus = ArTrackable_acquireNewAnchor(session as *mut ArSession,
                                                                           ar_trackable as *mut ArTrackable,
                                                                           ar_pose as *mut ArPose,
                                                                           &mut ar_anchor as *mut *mut ArAnchor);
                if ar_status != AR_SUCCESS as i32 {
                    log::e(&format!("arcore::augmented_face::track_faces ArTrackable_acquireNewAnchor ar_status: {}", &ar_status));
                    return;
                }

                let mut color = [0.0, 0.0, 0.0, 0.0];

                color[0] = 139.0;
                color[1] = 195.0;
                color[2] = 74.0;
                color[3] = 255.0;

                let colored_anchor = ::ColoredAnchor { anchor: ar_anchor, color };

                log::d(&format!("arcore::augmented_face::track_faces i : {}, colored_anchor : {:?}", i, &colored_anchor));

                map.insert(i, colored_anchor);
            }
        }

        ArTrackableList_destroy(ar_trackable_list);
        ar_trackable_list = ::std::ptr::null_mut();
    }
}
