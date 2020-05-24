use std::env;
use std::fs::File;
use std::io::Read;

use rgb::*;
use sparkle::gl::*;

use log;

pub fn load_shader(gl: &Gl, shader_type: GLenum, shader_source: &[u8]) -> GLuint {
    unsafe {
        let mut shader = gl.create_shader(shader_type);

        log::d(&format!("arcore::util::load_shader : shader = {}", shader));

        if shader == 0 {
            return shader;
        }
        gl.shader_source(shader, &[shader_source]);
        gl.compile_shader(shader);

        let mut compiled = [0];
        gl.get_shader_iv(shader, COMPILE_STATUS, &mut compiled);

        log::d(&format!("arcore::util::load_shader : compiled = {}", compiled[0]));

        if compiled[0] == 0 {
            let mut info_len = [0];
            gl.get_shader_iv(shader, INFO_LOG_LENGTH, &mut info_len);

            log::d(&format!("arcore::util::load_shader : info_len = {}", info_len[0]));

            if info_len[0] == 0 {
                return shader;
            }

            gl.delete_shader(shader);
            shader = 0;
        }

        shader
    }
}

pub fn create_program(gl: &Gl, vertex_source: &[u8], fragment_source: &[u8]) -> GLuint {
    unsafe {
        let vertex_shader = load_shader(gl, VERTEX_SHADER, vertex_source);

        log::d(&format!("arcore::util::create_program : vertex_shader = {}", vertex_shader));

        if vertex_shader == 0 {
            return 0;
        }

        let fragment_shader = load_shader(gl, FRAGMENT_SHADER, fragment_source);

        log::d(&format!("arcore::util::create_program : fragment_shader = {}", fragment_shader));

        if fragment_shader == 0 {
            return 0;
        }

        let mut program = gl.create_program();

        log::d(&format!("arcore::util::create_program : program = {}", program));

        if program != 0 {
            gl.attach_shader(program, vertex_shader);
            gl.attach_shader(program, fragment_shader);
            gl.link_program(program);

            let mut link_status = [0];
            gl.get_program_iv(program, LINK_STATUS, &mut link_status);

            log::d(&format!("arcore::util::create_program : link_status = {}", link_status[0]));

            if link_status[0] == 0 {
                gl.delete_program(program);
                program = 0;
            }
        }

        program
    }
}

pub fn read_image_from_resource(file: &str) {
    let mut path = env::current_exe().unwrap();
    log::d(&format!("arcore::util::read_image_from_resource path0 =  {:?}", &path));
    path = path.canonicalize().unwrap();
    log::d(&format!("arcore::util::read_image_from_resource path1 =  {:?}", &path));
    while path.pop() {
        path.push("resources");
        if path.is_dir() {
            break;
        }
        path.pop();
    }
    path.push(file);

    let mut state = ::lodepng::State::new();
//    state.remember_unknown_chunks(true);

    log::d(&format!("arcore::util::read_image_from_resource path =  {:?}", path));
    match state.decode_file(&path) {
        Ok(image) =>
            match image {
                ::lodepng::Image::RGBA(bitmap) => {
                    let image_u8: &[u8] = bitmap.buffer.as_bytes();
                    log::d(&format!("arcore::util::read_image_from_resource width = {} height = x {}", &bitmap.width, &bitmap.height));
                    log::d(&format!("arcore::util::read_image_from_resource path =  {:?}", image_u8));
                }
                other => {
                    log::e(&format!("arcore::util::read_image_from_resource  Could not load file , other: {:?}", other));
                }
            },
        Err(e) => {
            log::e(&format!("arcore::util::read_image_from_resource  Could not load file , because: {}", e));
        }
    }
}

// Read file from resourse
pub fn read_file_from_resource(file: &str) -> Vec<u8> {
    let mut path = env::current_exe().unwrap();
    path = path.canonicalize().unwrap();
    while path.pop() {
        path.push("resources");
        if path.is_dir() {
            break;
        }
        path.pop();
    }
    path.push(file);
    let mut buffer = vec![];
    File::open(path)
        .expect(&format!("Can't find file: {}", file))
        .read_to_end(&mut buffer)
        .expect("Can't read file");
    buffer
}


pub fn convert_rgba_to_grayscale(
    image_pixel_buffer: *mut u8,
    width: u32,
    height: u32,
    stride: u32,
    out_grayscale_buffer: *mut *mut u8,
) {
    log::d("arcore::util::convert_rgba_to_grayscale");

    let grayscale_stride = stride / 4;  // Only support RGBA_8888 format
    let grayscale_buffer_len = grayscale_stride * height;
    let mut grayscale_buffer: *mut u8 = Vec::with_capacity(grayscale_buffer_len as usize).as_mut_ptr();
    log::e(&format!("arcore::util::convert_rgba_to_grayscale grayscale_stride : {:?}", &grayscale_stride));
    for h in 0..height {
        for w in 0..width {
            let pixel = unsafe { image_pixel_buffer.offset((w * 4 + h * stride) as isize) };

            let r = unsafe { *pixel.offset(0 as isize) as f32};
            let g = unsafe { *pixel.offset(1 as isize) as f32};
            let b = unsafe { *pixel.offset(2 as isize) as f32};

//            if w == 100 && h % 10 == 0 {
//                log::d(&format!("arcore::util::convert_rgba_to_grayscale r : {:?}", &r));
//                log::d(&format!("arcore::util::convert_rgba_to_grayscale g : {:?}", &g));
//                log::d(&format!("arcore::util::convert_rgba_to_grayscale b : {:?}", &b));
//                log::d(&format!("arcore::util::convert_rgba_to_grayscale buffer i = {}, v = {:?}", w + h * grayscale_stride, (0.213 * r + 0.715 * g + 0.072 * b) as u8));
//            }
            unsafe { *grayscale_buffer.offset((w + h * grayscale_stride) as isize) = (0.213 * r + 0.715 * g + 0.072 * b) as u8 };

        }
    }
    log::e(&format!("arcore::util::convert_rgba_to_grayscale grayscale_buffer : {:?}", unsafe { *grayscale_buffer.offset(521520 as isize) } ));
    unsafe { *out_grayscale_buffer = grayscale_buffer as *mut u8 };
}

pub fn from_slice(bytes: &[f32]) -> [f32; 16] {
    let mut array = [0.0; 16];
    let bytes = &bytes[..array.len()]; // panics if not enough data
    array.copy_from_slice(bytes);
    array
}

pub fn get_array_from_mat4(mat: glm::Mat4) -> [f32; 16]{
    let mat_array_vec4 = mat.as_array();
    let mut mat_array: Vec<f32> = Vec::new();
    for i in 0..mat_array_vec4.len() {
        for j in 0..4 {
            mat_array.push(mat_array_vec4[i][j]);
        }
    }
    from_slice(mat_array.as_slice())
}

pub fn get_mat4_from_array(array: [f32; 16]) -> glm::Mat4 {
    glm::mat4(array[0], array[1], array[2], array[3],
              array[4], array[5], array[6], array[7],
              array[8], array[9], array[10], array[11],
              array[12], array[13], array[14], array[15])
}