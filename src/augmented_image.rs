use ffi_arcore::*;

use android_injected_glue::write_log;

use jni_interface;
use util;

pub fn create_augmented_image_database(ar_session: *const ArSession) -> *mut ArAugmentedImageDatabase {
    unsafe {
        let mut ar_augmented_image_database: *mut ArAugmentedImageDatabase = ::std::ptr::null_mut();
        ArAugmentedImageDatabase_create(ar_session, &mut ar_augmented_image_database);


        // load_image_from_asset
        let mut width = 0;
        let mut height = 0;
        let mut stride = 0;
        let mut index = 0;

        let mut image_pixel_buffer: *mut u8 = ::std::ptr::null_mut();
        let k_sample_image_name: &str = "default.jpg";
        let load_image_result = jni_interface::load_image_from_assets(
            k_sample_image_name,
            &mut width,
            &mut height,
            &mut stride,
            &mut image_pixel_buffer,
        );

        if !load_image_result {
            write_log(&format!("arcore_jni::lib::create_augmented_image_database load image failed: {}", &load_image_result));
        } else {
            write_log(&format!("arcore_jni::lib::create_augmented_image_database load image width = {}, height = {}, stride = {}, image_pixel_buffer = {:?}", &width, &height, &stride, &image_pixel_buffer));
        }


        // convert_rgba_to_grayscale
        let mut grayscale_buffer: *mut u8 = ::std::ptr::null_mut();
        util::convert_rgba_to_grayscale(image_pixel_buffer, width, height, stride, &mut grayscale_buffer);

        write_log(&format!("arcore_jni::lib::create_augmented_image_database ar_session : {:?}", &ar_session));
        write_log(&format!("arcore_jni::lib::create_augmented_image_database k_sample_image_name : {:?}", &k_sample_image_name));
        write_log(&format!("arcore_jni::lib::create_augmented_image_database grayscale_buffer : {:?}", &grayscale_buffer));
        write_log(&format!("arcore_jni::lib::create_augmented_image_database width : {:?}", &width));
        write_log(&format!("arcore_jni::lib::create_augmented_image_database height : {:?}", &height));
        write_log(&format!("arcore_jni::lib::create_augmented_image_database index : {:?}", &index));

        // add image to ArAugmentedImageDatabase
        let grayscale_stride = stride / 4;
        write_log(&format!("arcore_jni::lib::create_augmented_image_database grayscale_stride : {:?}", &grayscale_stride));
        let ar_status: ArStatus = ArAugmentedImageDatabase_addImage(
            ar_session,
            ar_augmented_image_database,
            k_sample_image_name.as_ptr(),
            grayscale_buffer as *const u8,
            width as i32,
            height as i32,
            grayscale_stride as i32,
            &mut index as *mut i32);

        if ar_status != AR_SUCCESS as i32 {
            write_log(&format!("arcore_jni::lib::create_augmented_image_database ArAugmentedImageDatabase_addImage failed: {}", &ar_status));
        }
        ar_augmented_image_database
    }
}
