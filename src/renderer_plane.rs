use std::error::Error;
use std::io::BufReader;
use std::fs::File;
use std::path::Path;

use ffi_arcore::*;
use rgb::*;
use sparkle::gl;

use jni_interface;
use log;
use util;

const VS_SRC: &'static [u8] = b"
precision highp float;
precision highp int;
attribute vec3 vertex;
varying vec2 v_textureCoords;
varying float v_alpha;

uniform mat4 mvp;
uniform mat4 model_mat;
uniform vec3 normal;

void main() {
  // Vertex Z value is used as the alpha in this shader.
  v_alpha = vertex.z;

  vec4 local_pos = vec4(vertex.x, 0.0, vertex.y, 1.0);
  gl_Position = mvp * local_pos;
  vec4 world_pos = model_mat * local_pos;

  // Construct two vectors that are orthogonal to the normal.
  // This arbitrary choice is not co-linear with either horizontal
  // or vertical plane normals.
  const vec3 arbitrary = vec3(1.0, 1.0, 0.0);
  vec3 vec_u = normalize(cross(normal, arbitrary));
  vec3 vec_v = normalize(cross(normal, vec_u));

  // Project vertices in world frame onto vec_u and vec_v.
  v_textureCoords = vec2(
  dot(world_pos.xyz, vec_u), dot(world_pos.xyz, vec_v));
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
    vertices_: Vec<::glm::Vec3>,
    triangles_: Vec<gl::types::GLushort>,
    model_mat_: [f32; 16],
    normal_vec_: ::glm::Vec3,

    texture_id_: gl::types::GLuint,

    shader_program_: gl::types::GLuint,
    attri_vertices_: gl::types::GLuint,
    uniform_mvp_mat_: gl::types::GLuint,
    uniform_texture_: gl::types::GLuint,
    uniform_model_mat_: gl::types::GLuint,
    uniform_normal_vec_: gl::types::GLuint,
    uniform_color_: gl::types::GLuint,
}

impl PlaneRenderer {
    pub fn new(gl: &gl::Gl) -> PlaneRenderer {
        log::i("arcore::renderer_plane::new");

        let shader_program = util::create_program(gl, VS_SRC, FS_SRC);

        if shader_program == 0 {
            log::e("arcore::renderer_plane::new Could not create program.");
        }

        let uniform_mvp_mat_ = gl.get_uniform_location(shader_program, "mvp") as u32;
        let uniform_texture_ = gl.get_uniform_location(shader_program, "texture") as u32;
        let uniform_model_mat_ = gl.get_uniform_location(shader_program, "model_mat") as u32;
        let uniform_normal_vec_ = gl.get_uniform_location(shader_program, "normal") as u32;
        let uniform_color_ = gl.get_uniform_location(shader_program, "color") as u32;
        let attri_vertices_ = gl.get_attrib_location(shader_program, "vertex") as u32;

        let texture_id = gl.gen_textures(1)[0];
        gl.bind_texture(gl::TEXTURE_2D, texture_id);
        gl.tex_parameter_i(gl::TEXTURE_2D, gl::TEXTURE_WRAP_S, gl::REPEAT as i32);
        gl.tex_parameter_i(gl::TEXTURE_2D, gl::TEXTURE_WRAP_T, gl::REPEAT as i32);

        gl.tex_parameter_i(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::LINEAR_MIPMAP_LINEAR as i32);
        gl.tex_parameter_i(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);

        jni_interface::load_png_from_assets(gl::TEXTURE_2D as i32, "trigrid.png");

        gl.generate_mipmap(gl::TEXTURE_2D);

        gl.bind_texture(gl::TEXTURE_2D, 0);

        PlaneRenderer {
            vertices_: Vec::new(),
            triangles_: Vec::new(),
            model_mat_: [0.0; 16],
            normal_vec_: ::glm::vec3(0.0, 0.0, 0.0),

            texture_id_: texture_id,

            shader_program_: shader_program,
            attri_vertices_: attri_vertices_,
            uniform_mvp_mat_: uniform_mvp_mat_,
            uniform_texture_: uniform_texture_,
            uniform_model_mat_: uniform_model_mat_,
            uniform_normal_vec_: uniform_normal_vec_,
            uniform_color_: uniform_color_,
        }
    }

    pub fn draw(&mut self, gl: &gl::Gl, projection_mat: ::glm::Mat4, view_mat: ::glm::Mat4, session: *const ArSession, plane: *const ArPlane, color: ::glm::Vec3) {
        log::d("arcore::renderer_plane::draw");

        if self.shader_program_ == 0 {
            log::e("arcore::renderer_plane::draw shader_program is null.");
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

        let mvp_array = util::get_array_from_mat4(projection_mat * view_mat * model_mat);

        gl.uniform_matrix_4fv(self.uniform_mvp_mat_ as i32, false, &mvp_array);
        gl.uniform_matrix_4fv(self.uniform_model_mat_ as i32, false, &self.model_mat_);

        gl.uniform_3f(self.uniform_normal_vec_ as i32, self.normal_vec_.x, self.normal_vec_.y, self.normal_vec_.z);
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
        log::d("arcore::renderer_plane::update_for_plane");

        unsafe {
            self.vertices_.clear();
            self.triangles_.clear();

            let mut polygon_length: i32 = 0;
            ArPlane_getPolygonSize(session, plane, &mut polygon_length as *mut i32);
            if polygon_length == 0 {
                log::e("arcore::renderer_plane::update_for_plane no valid plane polygon is found.");
            }

            let mut vertices_size: usize = (polygon_length / 2) as usize;

            log::d(&format!("arcore::renderer_plane::update_for_plane polygon_length = {}, vertices_size =  {}", polygon_length, vertices_size));

            let mut raw_vertices: Vec<glm::Vec2> = Vec::with_capacity(vertices_size);
            raw_vertices.set_len(vertices_size);

            ArPlane_getPolygon(session, plane, raw_vertices.as_mut_ptr() as *mut f32);

            for i in 0..vertices_size {
                self.vertices_.push(::glm::vec3(raw_vertices[i].x, raw_vertices[i].y, 0.0))
            }

            let mut pose: *mut ArPose = ::std::ptr::null_mut();
            ArPose_create(session, 0 as *const f32, &mut pose);
            ArPlane_getCenterPose(session, plane, pose);
            ArPose_getMatrix(session, pose as *const ArPose,
                             self.model_mat_.as_mut_ptr());

//            normal_vec_ = get_plane_normal(session, pose);

            // Feather distance 0.2 meters.
            let k_feather_length = 0.2;
            // Feather scale over the distance between plane center and vertices.
            let k_feather_scale = 0.2;

            // Fill vertex 0 to 3, with alpha set to 1.
            for i in 0..vertices_size {
                // Vector from plane center to current point.
                let v: ::glm::Vec2 = raw_vertices[i];
                let kf = k_feather_length / ::glm::length(v);
                let mut kfinal = 1.0;
                if kf < k_feather_scale {
                    kfinal = kf;
                } else {
                    kfinal = k_feather_scale;
                }
                let scale = 1.0 - kfinal;
                let result_v: ::glm::Vec2 = v * scale;
                self.vertices_.push(::glm::vec3(result_v.x, result_v.y, 1.0));
            }

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
//            log::d(&format!("arcore::PlaneRenderer::update_for_plane self.vertices_ : {:?}", self.vertices_));
//            log::d(&format!("arcore::PlaneRenderer::update_for_plane self.triangles_ : {:?}", self.triangles_));
        }
    }
}


//pub fn get_plane_normal(session: *const ArSession, pose: *const ArPose) -> glm::Vec3 {
//
//    let mut plane_pose_raw = [0.0; 7];
//    unsafe {
//        ArPose_getPoseRaw(session, pose, plane_pose_raw.as_mut_ptr());
//    }
//
//    let plane_quaternion: nalgebra_glm::Quat = nalgebra_glm::quat(plane_pose_raw[3], plane_pose_raw[0],
//                                                                  plane_pose_raw[1], plane_pose_raw[2]);
//
//    nalgebra_glm::quat_rotate_vec3(&plane_quaternion, &glm::vec3(0.0, 1.0, 0.0))
//}