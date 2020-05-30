use std::ffi::CStr;
use std::mem;
use std::ptr;
use std::rc::Rc;

use ffi_arcore::*;
use sparkle::gl;

use log;
use util;

pub const TEXTURE_EXTERNAL_OES: gl::ffi::types::GLenum = 0x8D65;

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
    shader_program_: gl::types::GLuint,
    texture_id_: gl::types::GLuint,
    attribute_vertices_: gl::types::GLuint,
    attribute_uvs_: gl::types::GLuint,
    uniform_texture_: gl::types::GLuint,
    transformed_uvs_: [f32; 8],
    uvs_initialized_: bool,
}

impl BackgroundRenderer {
    pub fn new(gl: &gl::Gl) -> BackgroundRenderer {
        log::i("arcore::background_renderer::new");
        unsafe {
            let shader_program = util::create_program(gl, VS_SRC, FS_SRC);

            if shader_program == 0 {
                log::e("arcore::background_renderer::new Could not create program.");
            }

            let texture_id = gl.gen_textures(1)[0];
            gl.bind_texture(TEXTURE_EXTERNAL_OES, texture_id);
            gl.tex_parameter_i(TEXTURE_EXTERNAL_OES, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl.tex_parameter_i(TEXTURE_EXTERNAL_OES, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);


            let uniform_texture = gl.get_uniform_location(shader_program, "texture") as u32;
            let attribute_vertices = gl.get_attrib_location(shader_program, "vertex") as u32;
            let attribute_uvs = gl.get_attrib_location(shader_program, "textureCoords") as u32;

            let transformed_uvs: [f32; 8] = [0.0; 8];

//            log::d(&format!("arcore::background_renderer::new shader_program : {}", shader_program));
//            log::d(&format!("arcore::background_renderer::new texture_id : {:?}", texture_id));
//            log::d(&format!("arcore::background_renderer::new uniform_texture : {}", uniform_texture));
//            log::d(&format!("arcore::background_renderer::new attribute_vertices : {}", attribute_vertices));
//            log::d(&format!("arcore::background_renderer::new attribute_uvs : {}", attribute_uvs));

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
        log::d("arcore::background_renderer::draw");
        unsafe {
            let mut x = 0;
            let geometry_changed: *mut i32 = &mut x;
            ArFrame_getDisplayGeometryChanged(session, frame, geometry_changed);

            log::d(&format!("arcore::background_renderer::draw geometry_changed : {}", *geometry_changed));

            if *geometry_changed != 0 || !self.uvs_initialized_ {
                ArFrame_transformDisplayUvCoords(session, frame, K_NUM_VERTICES * 2, &K_UVS as *const f32, self.transformed_uvs_.as_mut_ptr());
                self.uvs_initialized_ = true;
                log::d(&format!("arcore::background_renderer::draw self.uvs_initialized_ : {}", self.uvs_initialized_));
            }

            gl.use_program(self.shader_program_);
            gl.depth_mask(false);

            gl.uniform_1i(self.uniform_texture_ as i32, 1);
            gl.active_texture(gl::TEXTURE1);
            gl.bind_texture(TEXTURE_EXTERNAL_OES, self.texture_id_);


            let vbo = gl.gen_buffers(2);

            gl.bind_buffer(gl::ARRAY_BUFFER, vbo[0]);
            gl::buffer_data(gl, gl::ARRAY_BUFFER, &K_VERTICES, gl::STATIC_DRAW);

            gl.enable_vertex_attrib_array(self.attribute_vertices_);
            gl.vertex_attrib_pointer(self.attribute_vertices_, 3, gl::FLOAT, false, 0, 0);

            gl.bind_buffer(gl::ARRAY_BUFFER, vbo[1]);
            gl::buffer_data(&*gl, gl::ARRAY_BUFFER, &self.transformed_uvs_, gl::STATIC_DRAW);

            gl.enable_vertex_attrib_array(self.attribute_uvs_);
            gl.vertex_attrib_pointer(self.attribute_uvs_, 2, gl::FLOAT, false, 0, 0);

            gl.draw_arrays(gl::TRIANGLE_STRIP, 0, 4);

            gl.use_program(0);
            gl.depth_mask(true);
        }
    }

    pub fn get_texture_id(&self) -> gl::types::GLuint {
        self.texture_id_.clone()
    }
}