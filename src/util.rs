use std::ptr;
use std::rc::Rc;

use gleam::gl::*;
use android_injected_glue::write_log;
use ffi_arcore::*;

pub fn load_shader(gl: &Gl, shader_type: GLenum, shader_source: &[u8]) -> GLuint {
    //    GLuint shader = glCreateShader(shader_type);
    //    if (!shader) {
    //    return shader;
    //    }
    //
    //    glShaderSource(shader, 1, &shader_source, nullptr);
    //    glCompileShader(shader);
    //    GLint compiled = 0;
    //    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    //
    //    if (!compiled) {
    //    GLint info_len = 0;
    //
    //    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &info_len);
    //    if (!info_len) {
    //    return shader;
    //    }
    //
    //    char* buf = reinterpret_cast<char*>(malloc(info_len));
    //    if (!buf) {
    //    return shader;
    //    }
    //
    //    glGetShaderInfoLog(shader, info_len, nullptr, buf);
    //    LOGE("hello_ar::util::Could not compile shader %d:\n%s\n", shader_type,
    //    buf);
    //    free(buf);
    //    glDeleteShader(shader);
    //    shader = 0;
    //    }
    //
    //    return shader;
    //    }
    unsafe {
        let mut shader = gl.create_shader(shader_type);

        write_log(&format!("arcore_jni::util::load_shader : shader = {}", shader));

        if shader == 0 {
            return shader;
        }
        gl.shader_source(shader, &[shader_source]);
        gl.compile_shader(shader);

        let mut compiled = 0;
        compiled = gl.get_shader_iv(shader, COMPILE_STATUS);

        write_log(&format!("arcore_jni::util::load_shader : compiled = {}", compiled));

        if compiled == 0 {
            let mut info_len = 0;
            info_len = gl.get_shader_iv(shader, INFO_LOG_LENGTH);

            write_log(&format!("arcore_jni::util::load_shader : info_len = {}", info_len));

            if info_len == 0 {
                return shader;
            }

            gl.delete_shader(shader);
            shader = 0;
        }

        shader
    }
}

pub fn create_program(gl: &Gl, vertex_source: &[u8], fragment_source: &[u8]) -> GLuint {
    //    GLuint vertexShader = LoadShader(GL_VERTEX_SHADER, vertex_source);
    //    if (!vertexShader) {
    //        return 0;
    //    }
    //
    //    GLuint fragment_shader = LoadShader(GL_FRAGMENT_SHADER, fragment_source);
    //    if (!fragment_shader) {
    //        return 0;
    //    }
    //
    //    GLuint program = glCreateProgram();
    //    if (program) {
    //        glAttachShader(program, vertexShader);
    //        CheckGlError("hello_ar::util::glAttachShader");
    //        glAttachShader(program, fragment_shader);
    //        CheckGlError("hello_ar::util::glAttachShader");
    //        glLinkProgram(program);
    //        GLint link_status = GL_FALSE;
    //        glGetProgramiv(program, GL_LINK_STATUS, &link_status);
    //        if (link_status != GL_TRUE) {
    //            GLint buf_length = 0;
    //            glGetProgramiv(program, GL_INFO_LOG_LENGTH, &buf_length);
    //            if (buf_length) {
    //                char* buf = reinterpret_cast<char*>(malloc(buf_length));
    //                if (buf) {
    //                    glGetProgramInfoLog(program, buf_length, nullptr, buf);
    //                    LOGE("hello_ar::util::Could not link program:\n%s\n", buf);
    //                    free(buf);
    //                }
    //            }
    //            glDeleteProgram(program);
    //            program = 0;
    //        }
    //    }
    //    return program;
    //}
    unsafe {
        let vertex_shader = load_shader(gl, VERTEX_SHADER, vertex_source);

        write_log(&format!("arcore_jni::util::create_program : vertex_shader = {}", vertex_shader));

        if vertex_shader == 0 {
            return 0;
        }

        let fragment_shader = load_shader(gl, FRAGMENT_SHADER, fragment_source);

        write_log(&format!("arcore_jni::util::create_program : fragment_shader = {}", fragment_shader));

        if fragment_shader == 0 {
            return 0;
        }

        let mut program = gl.create_program();

        write_log(&format!("arcore_jni::util::create_program : program = {}", program));

        if program != 0 {
            gl.attach_shader(program, vertex_shader);
            gl.attach_shader(program, fragment_shader);
            gl.link_program(program);

            let mut link_status = 0;
            link_status = gl.get_program_iv(program, LINK_STATUS);

            write_log(&format!("arcore_jni::util::create_program : link_status = {}", link_status));

            if link_status == 0 {
                gl.delete_program(program);
                program = 0;
            }
        }

        //        let vs = gl.create_shader(VERTEX_SHADER);
        //
        //        gl.shader_source(vs, &[vertex_source]);
        //        gl.compile_shader(vs);
        //
        //        let fs = gl.create_shader(FRAGMENT_SHADER);
        //        gl.shader_source(fs, &[fragment_source]);
        //        gl.compile_shader(fs);
        //
        //        let program = gl.create_program();
        //        gl.attach_shader(program, vs);
        //        gl.attach_shader(program, fs);
        //        gl.link_program(program);

        program
    }
}


pub fn get_transform_matrix_from_anchor(session: *mut ArSession, anchor: *mut ArAnchor, out_model_mat: *mut f32) {
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