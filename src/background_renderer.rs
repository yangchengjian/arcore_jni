use std::ffi::CStr;
use std::mem;
use std::ptr;
use std::rc::Rc;

use android_injected_glue::write_log;
use ffi_arcore::*;
use gleam::gl;
use util;

pub const K_NUM_VERTICES: i32 = 4;

const VS_SRC: &'static [u8] = b"
    attribute vec4 vertex;
    attribute vec2 textureCoords;
    varying vec2 v_textureCoords;
    void main() {
    v_textureCoords = textureCoords;
        gl_Position = vertex;
    }
\0";

const FS_SRC: &'static [u8] = b"
    #extension GL_OES_EGL_image_external : require
    precision mediump float;
    uniform samplerExternalOES texture;
    varying vec2 v_textureCoords;
    void main() {
      gl_FragColor = texture2D(texture, v_textureCoords);
    }
\0";

const K_VERTICES: [f32; 12] = [
    -1.0, -1.0, 0.0,
    1.0, -1.0, 0.0,
    -1.0, 1.0, 0.0,
    1.0, 1.0, 0.0,
];

// UVs of the quad vertices (S, T)
const K_UVS: [f32; 8] = [
    0.0, 1.0,
    1.0, 1.0,
    0.0, 0.0,
    1.0, 0.0,
];

#[derive(Clone, Debug)]
pub struct BackgroundRenderer {
    //    pub gl: Rc<gl::Gl>,
    shader_program_: gl::types::GLuint,
    texture_id_: gl::types::GLuint,
    attribute_vertices_: gl::types::GLuint,
    attribute_uvs_: gl::types::GLuint,
    uniform_texture_: gl::types::GLuint,
    transformed_uvs_: [f32; 8],
    uvs_initialized_: bool,
}

impl BackgroundRenderer {
    pub fn initializel_content(gl: &gl::Gl) -> BackgroundRenderer {
        write_log("servo_jni::BackgroundRenderer::initializel_content");
        unsafe {
            let shader_program = util::create_program(gl, VS_SRC, FS_SRC);

            if shader_program == 0 {
                write_log("Could not create program.");
            }

            let texture_id = gl.gen_textures(1)[0];
            gl.bind_texture(gl::TEXTURE_EXTERNAL_OES, texture_id);
            gl.tex_parameter_i(gl::TEXTURE_EXTERNAL_OES, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl.tex_parameter_i(gl::TEXTURE_EXTERNAL_OES, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);


            let uniform_texture = gl.get_uniform_location(shader_program, "texture") as u32;
            let attribute_vertices = gl.get_attrib_location(shader_program, "vertex") as u32;
            let attribute_uvs = gl.get_attrib_location(shader_program, "textureCoords") as u32;

            let transformed_uvs: [f32; 8] = [0.0; 8];

            write_log(&format!("shader_program : {}", shader_program));
            write_log(&format!("texture_id : {:?}", texture_id));
            write_log(&format!("uniform_texture : {}", uniform_texture));
            write_log(&format!("attribute_vertices : {}", attribute_vertices));
            write_log(&format!("attribute_uvs : {}", attribute_uvs));

            BackgroundRenderer {
                shader_program_: shader_program,
                texture_id_: texture_id,
                attribute_vertices_: attribute_vertices,
                attribute_uvs_: attribute_uvs,
                uniform_texture_: uniform_texture,
                transformed_uvs_: transformed_uvs,
                uvs_initialized_: false,
            }
        }
    }

    pub fn draw(&mut self, gl: &gl::Gl, session: *const ArSession, frame: *const ArFrame) {
        write_log("servo_jni::BackgroundRenderer::draw");
        unsafe {
            write_log(&format!("self.shader_program_ : {}", self.shader_program_));
            write_log(&format!("self.texture_id_ : {:?}", self.texture_id_));
            write_log(&format!("self.uniform_texture_ : {}", self.uniform_texture_));
            write_log(&format!("self.attribute_vertices_ : {}", self.attribute_vertices_));
            write_log(&format!("self.attribute_uvs_ : {}", self.attribute_uvs_));
            write_log(&format!("self.uvs_initialized : {}", self.uvs_initialized_));

            let mut x = 0;
            let geometry_changed: *mut i32 = &mut x;
            ArFrame_getDisplayGeometryChanged(session, frame, geometry_changed);

            write_log(&format!("geometry_changed : {}", *geometry_changed));

            if (*geometry_changed != 0 || !self.uvs_initialized_) {
                ArFrame_transformDisplayUvCoords(session, frame, K_NUM_VERTICES * 2, &K_UVS as *const f32, self.transformed_uvs_.as_mut_ptr());
                self.uvs_initialized_ = true;
                write_log(&format!("self.uvs_initialized_ : {}", self.uvs_initialized_));
            }

            gl.use_program(self.shader_program_);
            gl.depth_mask(false);

            gl.uniform_1i(self.uniform_texture_ as i32, 0);
            gl.active_texture(gl::TEXTURE0);
            gl.bind_texture(gl::TEXTURE_EXTERNAL_OES, self.texture_id_);

//            gl.tex_image_2d(gl::TEXTURE_2D, 0, gl::RGB as gl::GLint, 2, 2, 0, gl::RGB, gl::UNSIGNED_BYTE, Some(&PIXELS));

            gl.enable_vertex_attrib_array(self.attribute_vertices_);
            gl.vertex_attrib_pointer_ptr(self.attribute_vertices_, 3, false, 0, K_VERTICES.as_ptr() as *const gl::GLvoid);

            gl.enable_vertex_attrib_array(self.attribute_uvs_);
            gl.vertex_attrib_pointer_ptr(self.attribute_uvs_, 2, false, 0, self.transformed_uvs_.as_ptr() as *const gl::GLvoid);


            //            let vbo = gl.gen_buffers(2);
            //            let vao = gl.gen_vertex_arrays(2);
            //
            //            gl.bind_buffer(gl::ARRAY_BUFFER, vbo[0]);
            //            gl::buffer_data(gl, gl::ARRAY_BUFFER, &K_VERTICES, gl::STATIC_DRAW);
            //
            //            gl.bind_vertex_array(vao[0]);
            //            gl.enable_vertex_attrib_array(self.attribute_vertices_);
            //            gl.vertex_attrib_pointer(self.attribute_vertices_, 3, gl::FLOAT, false, 0, 0);
            //
            //            gl.bind_buffer(gl::ARRAY_BUFFER, vbo[1]);
            //            gl::buffer_data(&*gl, gl::ARRAY_BUFFER, &self.transformed_uvs_, gl::STATIC_DRAW);
            //
            //            gl.bind_vertex_array(vao[1]);
            //            gl.enable_vertex_attrib_array(self.attribute_uvs_);
            //            gl.vertex_attrib_pointer(self.attribute_uvs_, 2, gl::FLOAT, false, 0, (3 * ::std::mem::size_of::<f64>()) as gl::GLuint);

            gl.draw_arrays(gl::TRIANGLE_STRIP, 0, 4);

            gl.use_program(0);
            gl.depth_mask(true);
        }
    }

    pub fn get_texture_id(&self) -> gl::types::GLuint {
        write_log("servo_jni::BackgroundRenderer::get_texture_id");
        self.texture_id_.clone()
    }
}