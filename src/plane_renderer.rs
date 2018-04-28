use std::error::Error;
use std::io::BufReader;
use std::fs::File;
use std::path::Path;

use android_injected_glue::write_log;
use ffi_arcore::*;
use gleam::gl;
use util;

const VS_SRC: &'static [u8] = b"
    precision highp float;
    precision highp int;
    attribute vec3 vertex;
    varying vec2 v_textureCoords;
    varying float v_alpha;

    uniform mat4 mvp;
    uniform mat4 texture_mat;
    void main() {
      gl_Position = mvp * vec4(vertex.x, 0.0, vertex.y, 1.0);
      // Vertex Z value is used as the alpha in this shader.
      v_alpha = vertex.z;
      v_textureCoords = (texture_mat * vec4(vertex.x, 0.0, vertex.y, 1.0)).xz;
    }
\0";

const FS_SRC: &'static [u8] = b"
    precision highp float;
    precision highp int;
    uniform sampler2D texture;
    uniform vec3 color;
    varying vec2 v_textureCoords;
    varying float v_alpha;
    void main() {
      float r = texture2D(texture, v_textureCoords).r;
      gl_FragColor = vec4(color.xyz, r * v_alpha);
    }
\0";

#[derive(Clone, Debug)]
pub struct PlaneRenderer {
    //    std::vector<glm::vec3> vertices_;
    //    std::vector<GLushort> triangles_;
    //    glm::mat4 model_mat_ = glm::mat4(1.0f);
    //
    //    GLuint vertex_buffers_[2];
    //    GLuint texture_id_;
    //
    //    GLuint shader_program_;
    //    GLuint attri_vertices_;
    //    GLuint uniform_mvp_mat_;
    //    GLuint uniform_texture_;
    //    GLuint uniform_texture_mat_;
    //    GLuint uniform_color_;
    vertices_: Vec<::glm::Vec3>,
    triangles_: Vec<gl::types::GLushort>,
    model_mat_: [f32; 16],

    vertex_buffers_: [gl::types::GLuint; 2],
    texture_id_: gl::types::GLuint,

    shader_program_: gl::types::GLuint,
    attri_vertices_: gl::types::GLuint,
    uniform_mvp_mat_: gl::types::GLuint,
    uniform_texture_: gl::types::GLuint,
    uniform_texture_mat_: gl::types::GLuint,
    uniform_color_: gl::types::GLuint,
}

impl PlaneRenderer {
    pub fn initializel_content(gl: &gl::Gl) -> PlaneRenderer {
        write_log("arcore_jni::plane_renderer::initializel_content");

        let shader_program = util::create_program(gl, VS_SRC, FS_SRC);

        if shader_program == 0 {
            write_log("arcore_jni::plane_renderer::initializel_content Could not create program.");
        }

        let uniform_mvp_mat_ = gl.get_uniform_location(shader_program, "mvp") as u32;
        let uniform_texture_ = gl.get_uniform_location(shader_program, "texture") as u32;
        let uniform_texture_mat_ = gl.get_uniform_location(shader_program, "texture_mat") as u32;
        let uniform_color_ = gl.get_uniform_location(shader_program, "color") as u32;

        let attri_vertices_ = gl.get_attrib_location(shader_program, "vertex") as u32;

        let texture_id = gl.gen_textures(1)[0];
        gl.bind_texture(gl::TEXTURE_2D, texture_id);
        gl.tex_parameter_i(gl::TEXTURE_2D, gl::TEXTURE_WRAP_S, gl::REPEAT as i32);
        gl.tex_parameter_i(gl::TEXTURE_2D, gl::TEXTURE_WRAP_T, gl::REPEAT as i32);

        gl.tex_parameter_i(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::LINEAR_MIPMAP_LINEAR as i32);
        gl.tex_parameter_i(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);

        //        let path = Path::new("assests/models/trigrid.png");
        //        let mut file = match File::open(&path) {
        //            Err(why) => panic!("couldn't open : {}", why.description()),
        //            Ok(file) => file,
        //        };
        //        let mut image = Vec::new();
        //        file.read(&mut image);


        gl.tex_image_2d(gl::TEXTURE_2D, 0, gl::RGB as gl::GLint, 2 as gl::GLsizei,
                        2 as gl::GLsizei, 0, gl::RGB, gl::UNSIGNED_BYTE, Some(&[0, 1]));

        gl.generate_mipmap(gl::TEXTURE_2D);

        gl.bind_texture(gl::TEXTURE_2D, 0);

        //        util::CheckGlError("plane_renderer::InitializeGlContent()");

        PlaneRenderer {
            vertices_: Vec::new(),
            triangles_: Vec::new(),
            model_mat_: [0.0; 16],

            vertex_buffers_: [0; 2],
            texture_id_: texture_id,

            shader_program_: shader_program,
            attri_vertices_: attri_vertices_,
            uniform_mvp_mat_: uniform_mvp_mat_,
            uniform_texture_: uniform_texture_,
            uniform_texture_mat_: uniform_texture_mat_,
            uniform_color_: uniform_color_,
        }
    }

    pub fn draw(&mut self, gl: &gl::Gl, projection_mat: ::glm::Mat4, view_mat: ::glm::Mat4, session: *const ArSession, plane: *const ArPlane, color: ::glm::Vec3) {
        write_log("arcore_jni::plane_renderer::on_draw");

        if self.shader_program_ == 0 {
            write_log("arcore_jni::plane_renderer::on_draw shader_program is null.");
            return;
        }

        self.update_for_plane(session, plane);

        gl.use_program(self.shader_program_);
        gl.depth_mask(false);

        gl.uniform_1i(self.uniform_texture_ as i32, 0);
        gl.active_texture(gl::TEXTURE0);

        gl.bind_texture(gl::TEXTURE_2D, self.texture_id_);

        // Compose final mvp matrix for this plane renderer.
        let mut model_mat = ::glm::mat4(self.model_mat_[0], self.model_mat_[1], self.model_mat_[2], self.model_mat_[3],
                                        self.model_mat_[4], self.model_mat_[5], self.model_mat_[6], self.model_mat_[7],
                                        self.model_mat_[8], self.model_mat_[9], self.model_mat_[10], self.model_mat_[11],
                                        self.model_mat_[12], self.model_mat_[13], self.model_mat_[14], self.model_mat_[15]);

        let mvp = projection_mat * view_mat * model_mat;
        let mvp_array_vec4 = mvp.as_array();

        let mut mvp_array: Vec<f32> = Vec::new();

        for i in 0..mvp_array_vec4.len() {
            for j in 0..4 {
                mvp_array.push(mvp_array_vec4[i][j]);
            }
        }

        gl.uniform_matrix_4fv(self.uniform_mvp_mat_ as i32, false, &mvp_array);

        gl.uniform_matrix_4fv(self.uniform_texture_mat_ as i32, false, &self.model_mat_);

        gl.uniform_3f(self.uniform_color_ as i32, color.x, color.y, color.z);

        let vbo = gl.gen_buffers(1);

        gl.bind_buffer(gl::ARRAY_BUFFER, vbo[0]);
        gl::buffer_data(gl, gl::ARRAY_BUFFER, &self.vertices_, gl::STATIC_DRAW);

        gl.enable_vertex_attrib_array(self.attri_vertices_);
        gl.vertex_attrib_pointer(self.attri_vertices_, 3, gl::FLOAT, false, 0, 0);

        gl.draw_elements(gl::TRIANGLES, self.triangles_.len() as i32, gl::UNSIGNED_SHORT, 0);

        gl.use_program(0);
        gl.depth_mask(true);
    }

    pub fn update_for_plane(&mut self, session: *const ArSession, plane: *const ArPlane) {
        write_log("arcore_jni::PlaneRenderer::update_for_plane");

        unsafe {
            self.vertices_.clear();
            self.triangles_.clear();

            let mut polygon_length: i32 = 0;
            ArPlane_getPolygonSize(session, plane, &mut polygon_length as *mut i32);
            if polygon_length == 0 {
                write_log("arcore_jni::PlaneRenderer::update_for_plane no valid plane polygon is found.");
            }

            let mut vertices_size: usize = (polygon_length / 2) as usize;

            write_log(&format!("arcore_jni::PlaneRenderer::update_for_plane polygon_length = {}, vertices_size =  {}", polygon_length, vertices_size));

            let mut raw_vertices_f32: Vec<f32> = Vec::with_capacity(polygon_length as usize);
            raw_vertices_f32.set_len(polygon_length as usize);

            ArPlane_getPolygon(session, plane, raw_vertices_f32.as_mut_ptr() as *mut f32);

            let mut raw_vertices: Vec<::glm::Vec2> = Vec::new();
            for i in 0..vertices_size {
                raw_vertices.push(::glm::vec2(raw_vertices_f32[2 * i], raw_vertices_f32[2 * i + 1]))
            }

            for i in 0..vertices_size {
                self.vertices_.push(::glm::vec3(raw_vertices_f32[2 * i], raw_vertices_f32[2 * i + 1], 0.0))
            }

            let mut pose: *mut ArPose = ::std::ptr::null_mut();
            ArPose_create(session, 0 as *const f32, &mut pose);
            ArPlane_getCenterPose(session, plane, pose);
            ArPose_getMatrix(session, pose as *const ArPose,
                             self.model_mat_.as_mut_ptr());


            // Get plane center in XZ axis.
            let plane_center: ::glm::Vec2 = ::glm::vec2(self.model_mat_[12], self.model_mat_[14]);
            // Feather distance 0.2 meters.
            let kFeatherLength = 0.2;
            // Feather scale over the distance between plane center and vertices.
            let kFeatherScale = 0.2;

            // Fill vertex 0 to 3, with alpha set to 1.
            for i in 0..vertices_size {
                let v: ::glm::Vec2 = raw_vertices[i];


                // Vector from plane center to current point.
                let d: ::glm::Vec2 = v - plane_center;
                let kf = kFeatherLength / ::glm::length(d);
                let mut kfinal = 1.0;
                if kf < kFeatherScale {
                    kfinal = kf;
                } else {
                    kfinal = kFeatherScale;
                }
                let scale = 1.0 - kfinal;
                let result_v: ::glm::Vec2 = d * scale + plane_center;
                self.vertices_.push(::glm::vec3(result_v.x, result_v.y, 1.0));
            }

            write_log("arcore_jni::PlaneRenderer::update_for_plane E");

            let vertices_length = self.vertices_.len();
            let half_vertices_length = vertices_length / 2;
            // Generate triangle (4, 5, 6) and (4, 6, 7).
            for i in half_vertices_length..vertices_length - 1 {
                self.triangles_.push(half_vertices_length as u16);
                self.triangles_.push(i as u16);
                self.triangles_.push((i + 1) as u16);
            }
            // Generate triangle (0, 1, 4), (4, 1, 5), (5, 1, 2), (5, 2, 6),
            // (6, 2, 3), (6, 3, 7), (7, 3, 0), (7, 0, 4)
            for i in 0..half_vertices_length {
                self.triangles_.push(i as u16);
                self.triangles_.push(((i + 1) % half_vertices_length) as u16);
                self.triangles_.push((i + half_vertices_length) as u16);

                self.triangles_.push((i + half_vertices_length) as u16);
                self.triangles_.push(((i + 1) % half_vertices_length) as u16);
                self.triangles_.push(((i + half_vertices_length + 1) % half_vertices_length + half_vertices_length) as u16);
            }
            write_log(&format!("arcore_jni::PlaneRenderer::update_for_plane self.vertices_ : {:?}", self.vertices_));
            write_log(&format!("arcore_jni::PlaneRenderer::update_for_plane self.triangles_ : {:?}", self.triangles_));
        }
    }
}