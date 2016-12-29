
//
// Minimal OpenGL demo of Implicit Skinning
// Johan Berdat and Laurent Valette
//

// Demo parameters
#define WIDTH 1280
#define HEIGHT 768
#define MULTISAMPLING 16
#define FPS 0
#define DEBUG 1

// GLEW (http://glew.sourceforge.net/)
#define GLEW_STATIC
#include <GL/glew.h>

// GLM (http://glm.g-truc.net/)
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/dual_quaternion.hpp>

// GLFW 3 (http://www.glfw.org/)
#include <GLFW/glfw3.h>

// OpenMesh (http://www.openmesh.org/)
#define OM_STATIC_BUILD
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

// Standard C++ Library
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#if FPS
#include <chrono>
#include <thread>
#endif

// Skinning code
#include "skinning.h"

// Debug log callback
#if DEBUG
static void APIENTRY debugCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, GLchar const * message, void const * userParam) {
  std::cout << message << std::endl;
}
#endif

// Program entry point
int main(int argc, char** argv) {
  srand(time(NULL));

  // Load mesh
  Implicit<HRBF> implicit;
  if (!OpenMesh::IO::read_mesh(implicit.mesh, "../matlab/mesh/neutral.obj")) {
    std::cerr << "failed to open mesh" << std::endl;
    return -1;
  }
  implicit.mesh.update_normals();
  
  // Load skeleton
  std::ifstream skeleton("../matlab/mesh/neutral.skl");
  std::vector<int> parents;
  std::vector<glm::vec3> centers;
  std::vector<std::vector<glm::mat4>> poses;
  poses.push_back(std::vector<glm::mat4>(18));
  while (skeleton) {
    std::string token;
    skeleton >> token;
    if (token == "w") {
      glm::vec4 w;
      glm::ivec4 i;
      skeleton >> i.x;
      skeleton.ignore();
      skeleton >> w.x >> i.y;
      skeleton.ignore();
      skeleton >> w.y >> i.z;
      skeleton.ignore();
      skeleton >> w.z >> i.w;
      skeleton.ignore();
      skeleton >> w.w;
      implicit.weights.push_back(w);
      implicit.weights_bones.push_back(i - glm::ivec4(1, 1, 1, 1));
    } else if (token == "b") {
      int p;
      glm::vec3 c;
      skeleton >> p >> p >> c.x >> c.y >> c.z;
      parents.push_back(p - 1);
      centers.push_back(c);
      implicit.sdfs.push_back(HRBF());
    } else if (token == "h") {
      glm::vec3 p;
      glm::vec4 c;
      skeleton >> p.x >> p.y >> p.z >> c.x >> c.y >> c.z >> c.w;
      implicit.sdfs[centers.size() - 1].centers.push_back(p);
      implicit.sdfs[centers.size() - 1].coefficients.push_back(c);
    } else if (token == "p") {
      poses.push_back(std::vector<glm::mat4>());
    } else if (token == "t") {
      glm::mat4 t;
      for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
          skeleton >> t[j][i];
      poses[poses.size() - 1].push_back(t);
    }
  }
  implicit.initialize();
  
  // Create OpenGL context
  glfwInit();
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
#if DEBUG
  glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, 1);
#endif
#if MULTISAMPLING
  glfwWindowHint(GLFW_SAMPLES, MULTISAMPLING);
#endif
  GLFWwindow * window = glfwCreateWindow(WIDTH, HEIGHT, "Demo", NULL, NULL);
  glfwMakeContextCurrent(window);
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  glewExperimental = GL_TRUE;
  glewInit();
  
  // Enable debug callback
#if DEBUG
  glEnable(GL_DEBUG_OUTPUT);
  glDebugMessageCallback(debugCallback, nullptr);
  glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, GL_TRUE);
  glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_NOTIFICATION, 0, nullptr, GL_FALSE);
  glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
#endif
  
  // Create vertex buffer
  GLuint vertices;
  glGenBuffers(1, &vertices);
  glBindBuffer(GL_ARRAY_BUFFER, vertices);
  glBufferData(GL_ARRAY_BUFFER, implicit.count * sizeof(float) * 3, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, (void *)0);
  
  // Create normal buffer
  GLuint normals;
  glGenBuffers(1, &normals);
  glBindBuffer(GL_ARRAY_BUFFER, normals);
  glBufferData(GL_ARRAY_BUFFER, implicit.count * sizeof(float) * 3, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, false, 0, (void *)0);
  
  // Create index buffer
  GLuint indices;
  glGenBuffers(1, &indices);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices);
  std::vector<GLuint> indices_data;
  indices_data.reserve(implicit.mesh.n_faces() * 3);
  for (Mesh::FaceIter f = implicit.mesh.faces_begin(); f != implicit.mesh.faces_end(); ++f)
    for (Mesh::FaceVertexIter v = implicit.mesh.fv_iter(*f); v.is_valid(); ++v)
      indices_data.push_back(v->idx());
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_data.size() * sizeof(GLuint), &indices_data[0], GL_STATIC_DRAW);
  
  // Create shader
  char const * vsrc = "\
    #version 330 core                                                           \n\
    layout(location = 0) in vec3 position;                                      \n\
    layout(location = 1) in vec3 normal;                                        \n\
    out vec3 nor;                                                               \n\
    uniform mat4 projection, world, object;                                     \n\
    void main() {                                                               \n\
      gl_Position = projection * world * object * vec4(position, 1.0);          \n\
      nor = (object * vec4(normal, 0.0)).xyz;                                   \n\
    }                                                                           \n\
  ";
  char const * fsrc = "\
    #version 330                                                                \n\
    in vec3 nor;                                                                \n\
    out vec4 color;                                                             \n\
    void main() {                                                               \n\
      float d = dot(normalize(nor), normalize(vec3(-1.0, -1.0, -1.0)));         \n\
      float x = d * 0.5 + 0.5;                                                  \n\
      color = vec4(x, x, x * 0.9 + 0.1, 1.0);                                   \n\
    }                                                                           \n\
  ";
  GLuint program = glCreateProgram();
  GLuint vs = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vs, 1, &vsrc, NULL);
  glCompileShader(vs);
  glAttachShader(program, vs);
  glDeleteShader(vs);
  GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fs, 1, &fsrc, NULL);
  glCompileShader(fs);
  glAttachShader(program, fs);
  glDeleteShader(fs);
  glLinkProgram(program);
  
  // Configure rendering
#if MULTISAMPLING
  glEnable(GL_MULTISAMPLE);
#endif
  glEnable(GL_DEPTH_TEST);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glUseProgram(program);
  glm::mat4 projection = glm::perspectiveFov(1.3f, (float)width, (float)height, 0.1f, 100.0f);
  glUniformMatrix4fv(glGetUniformLocation(program, "projection"), 1, false, glm::value_ptr(projection));
  glm::mat4 object = glm::transpose(glm::mat4(
    // TODO compute this in C++ (mean, std, pca)
     0.5584,  0.1422,  0.2803, -2.2876,
    -0.2289,  0.5755,  0.1641, -0.9378,
    -0.2154, -0.2431,  0.5523, -6.4991,
     0     ,  0     ,  0     ,  1
  ));
  glUniformMatrix4fv(glGetUniformLocation(program, "object"), 1, false, glm::value_ptr(object));
  float camh = 0, camd = 8;
  glm::mat4 world = glm::lookAt(glm::vec3(camd * std::cos(camh), camd * std::sin(camh), 1.0f), glm::vec3(0.0f, 0.0f, -4.0f), glm::vec3(0.0f, 0.0f, 1.0f));
  glUniformMatrix4fv(glGetUniformLocation(program, "world"), 1, false, glm::value_ptr(world));
  
  // Main loop
  int previous = 0;
  int next = 1;
  float alpha = 0;
  int mode = 1;
  double mousex = 0, mousey = 0, mousepx, mousepy;
  glfwSetTime(0.01);
  double now = 0.0, before = 0.0;
  double lastfps = 0.0, averagedelta = 0.0;
  do {
  
    // Update time and handle FPS capper
    before = now;
    now = glfwGetTime();
    double delta = now - before;
#if FPS
    if (delta < 1.0 / FPS) {
      std::this_thread::sleep_for(std::chrono::duration<double>(1.0 / FPS - delta));
      now = glfwGetTime();
      delta = now - before;
      glfwPollEvents();
    }
#endif
    averagedelta = delta * 0.05 + averagedelta * 0.95;
    if (now - lastfps > 1.0) {
      lastfps = now;
      std::cout << (1.0 / averagedelta) << " FPS" << std::endl;
    }
    
    // Update mouse
    mousepx = mousex;
    mousepy = mousey;
    glfwGetCursorPos(window, &mousex, &mousey);
    
    // Update camera
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) != GLFW_RELEASE) {
      camh += 4.0 * (mousex - mousepx) / WIDTH;
      camd += 8.0 * (mousey - mousepy) / HEIGHT;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT) != GLFW_RELEASE)
      camh -= delta;
    if (glfwGetKey(window, GLFW_KEY_RIGHT) != GLFW_RELEASE)
      camh += delta;
    if (glfwGetKey(window, GLFW_KEY_UP) != GLFW_RELEASE)
      camd -= delta * 4;
    if (glfwGetKey(window, GLFW_KEY_DOWN) != GLFW_RELEASE)
      camd += delta * 4;
    world = glm::lookAt(glm::vec3(camd * std::cos(camh), camd * std::sin(camh), 1.0f), glm::vec3(0.0f, 0.0f, -4.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    glUniformMatrix4fv(glGetUniformLocation(program, "world"), 1, false, glm::value_ptr(world));
    
    // Update mode
    if (glfwGetKey(window, GLFW_KEY_1) != GLFW_RELEASE)
      mode = 1;
    if (glfwGetKey(window, GLFW_KEY_2) != GLFW_RELEASE)
      mode = 2;
    if (glfwGetKey(window, GLFW_KEY_3) != GLFW_RELEASE)
      mode = 3;
  
    // Update transforms
    alpha += delta * 0.3;
    if (alpha >= 1) {
      alpha -= 1;
      previous = next;
      while (next == previous)
        next = rand() % poses.size();
    }
    for (unsigned i = 0; i < 18; ++i) {
      glm::quat q1(glm::quat_cast(poses[previous][i]));
      glm::quat q2(glm::quat_cast(poses[next][i]));
      glm::quat q(glm::slerp(q1, q2, alpha));
      glm::mat4 r(glm::mat3_cast(q));
      if (parents[i] >= 0) {
        glm::mat4 t, it;
        glm::vec3 c = centers[i];
        t[3].x = c.x;
        t[3].y = c.y;
        t[3].z = c.z;
        it[3].x = -c.x;
        it[3].y = -c.y;
        it[3].z = -c.z;
        implicit.transforms[i] = implicit.transforms[parents[i]] * t * r * it;
      }
      else
        implicit.transforms[i] = r; // TODO add translation
    }
    
    // Skin
    switch (mode) {
    case 1:
      implicit.skin_linear();
      break;
    case 2:
      implicit.skin_dualquat();
      break;
    case 3:
      implicit.skin_implicit();
    }
    
    // Clear screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Draw mesh
    glBindBuffer(GL_ARRAY_BUFFER, vertices);
    glBufferSubData(GL_ARRAY_BUFFER, 0, implicit.count * sizeof(float) * 3, implicit.vertices);
    glBindBuffer(GL_ARRAY_BUFFER, normals);
    glBufferSubData(GL_ARRAY_BUFFER, 0, implicit.count * sizeof(float) * 3, implicit.normals);
    glDrawElements(GL_TRIANGLES, indices_data.size(), GL_UNSIGNED_INT, (void*)0);
    
    // Swap
    glfwSwapBuffers(window);
    glfwPollEvents();
  } while (!glfwWindowShouldClose(window));
  
  // Release resources
  glDeleteProgram(program);
  glDeleteBuffers(1, &vertices);
  glDeleteBuffers(1, &normals);
  glDeleteBuffers(1, &indices);
  glfwDestroyWindow(window);
  glfwTerminate();
  return 0;
}
