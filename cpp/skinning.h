
//
// Minimal implementation of Implicit Skinning
// Johan Berdat an Laurent Valette
//

#ifndef SKINNING_H
#define SKINNING_H

// GLM (http://glm.g-truc.net/)
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtx/dual_quaternion.hpp>

// OpenMesh (http://www.openmesh.org/)
#define OM_STATIC_BUILD
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

// Standard C++ Library
#include <cmath>
#include <vector>


// Hermite Radial Basis Function using triharmonic kernel
struct HRBF {
  
  // Precomputed values
  std::vector<glm::vec3> centers;
  std::vector<glm::vec4> coefficients;
  
  // Evaluate signed distance function
  float apply(glm::vec3 const & point, glm::vec3 & gradient) const {
    
    // For each center
    float distance = 0;
    gradient = glm::vec3();
    for (unsigned i = 0; i < centers.size(); ++i) {
    
      // Get coefficients
      float alpha = coefficients[i].w;
      glm::vec3 beta(coefficients[i]);
      
      // Compute difference
      glm::vec3 delta = point - centers[i];
      float norm = glm::length(delta);
      
      // Add contribution
      distance += norm * (alpha * norm * norm + 3.0f * glm::dot(beta, delta));
      gradient += (alpha * norm + glm::dot(beta, delta) / norm) * delta + norm * beta;
      
    }
    gradient *= 3.0f;
    return distance;
    
  }
  
};


// Single-precision OpenMesh mesh with per-vertex normals
struct Traits : public OpenMesh::DefaultTraits {
  VertexAttributes(OpenMesh::Attributes::Normal);
  FaceAttributes(OpenMesh::Attributes::Normal);
};
typedef OpenMesh::TriMesh_ArrayKernelT<Traits> Mesh;


// Implicit skinning data structure
template<typename SDF>
struct Implicit {
  
  // Inputs at initialization
  Mesh mesh;
  std::vector<glm::vec4> weights;
  std::vector<glm::ivec4> weights_bones;
  std::vector<SDF> sdfs;
  float threshold = 2.0f;
  
  // Inputs before skinning
  std::vector<glm::mat4> transforms;
  
  // Internal data
  unsigned count = 0;
  unsigned bones = 0;
  glm::vec3 * vertices;
  glm::vec3 * normals;
  std::vector<glm::vec3> original_vertices;
  std::vector<glm::vec3> original_normals;
  std::vector<float> original_offsets;
  std::vector<std::vector<unsigned> > neighbors;
  std::vector<std::vector<float> > neighbors_weights;
  std::vector<float> offsets;
  std::vector<glm::vec3> previous_gradients;
  std::vector<float> betas;
  std::vector<float> averages;
  std::vector<glm::vec3> centroids;
  
  // Setup internal data
  void initialize() {
  
    // Get mesh infos
    count = mesh.n_vertices();
    bones = sdfs.size();
    vertices = (glm::vec3*)mesh.points();
    normals = (glm::vec3*)&mesh.normal(*mesh.vertices_begin());
    
    // Copy original vertices
    original_vertices.assign(vertices, vertices + count);
    original_normals.assign(normals, normals + count);
    
    // Set transforms to identity
    transforms.assign(bones, glm::mat4());
    
    // Compute distances at rest pose
    original_offsets.resize(count);
    glm::vec3 dummy;
    for (unsigned i = 0; i < count; ++i)
      original_offsets[i] = distance(i, dummy);
    
    // Get one-ring neighbors
    neighbors.resize(count);
    for (unsigned i = 0; i < count; ++i) {
      Mesh::VertexHandle v(i);
      for (Mesh::VertexVertexIter n = mesh.vv_iter(v); n.is_valid(); ++n)
        neighbors[i].push_back(n->idx());
    }
    
    // Compute mean value coordinates
    neighbors_weights.resize(count);
    for (unsigned i = 0; i < count; ++i) {
      
      // Project neighbors on tangent plane
      std::vector<glm::vec3> projections;
      for (unsigned j = 0; j < neighbors[i].size(); ++j) {
        glm::vec3 delta = vertices[neighbors[i][j]] - vertices[i];
        projections.push_back(delta - normals[i] * glm::dot(normals[i], delta));
      }
      
      // Compute angles
      std::vector<float> angles;
      for (unsigned j = 0; j < projections.size(); ++j) {
        float cosine = glm::dot(glm::normalize(projections[j]), glm::normalize(projections[(j + 1) % projections.size()]));
        angles.push_back(std::acos(cosine));
      }
      
      // Compute barycentric coordinates
      float sum = 0;
      for (unsigned j = 0; j < projections.size(); ++j) {
        float length = glm::length(projections[j]);
        float tan1 = std::tan(angles[(j + projections.size() - 1) % projections.size()] * 0.5f);
        float tan2 = std::tan(angles[j] * 0.5f);
        float weight = (tan1 + tan2) / length;
        neighbors_weights[i].push_back(weight);
        sum += weight;
      }
      
      // Normalize weights
      for (unsigned j = 0; j < projections.size(); ++j)
        neighbors_weights[i][j] /= sum;
    }
    
    // Allocate buffers
    offsets.resize(count);
    betas.resize(count);
    previous_gradients.resize(count);
    averages.resize(count);
    centroids.resize(count);
    
  }

  // Apply linear transform blending to update vertices and normals
  void skin_linear() {
  
    // For each vertex
    for (unsigned i = 0; i < count; ++i) {
      
      // Compute weighted transform
      glm::mat4 transform =
        weights[i].x * transforms[weights_bones[i].x] +
        weights[i].y * transforms[weights_bones[i].y] +
        weights[i].z * transforms[weights_bones[i].z] +
        weights[i].w * transforms[weights_bones[i].w];
      
      // Transform vertex and normal
      vertices[i] = glm::vec3(transform * glm::vec4(original_vertices[i], 1));
      normals[i] = glm::vec3(transform * glm::vec4(original_normals[i], 0));
      
    }
  }
  
  // Apply dual quaternion blending to update vertices and normals
  void skin_dualquat() {
    
    // Convert rotations to dual quaternions
    std::vector<glm::dualquat> duals(bones);
    for (unsigned i = 0; i < bones; ++i) {
      glm::quat q(glm::quat_cast(transforms[i]));
      glm::vec3 t(transforms[i][3]);
      glm::dualquat d(q, t);
      duals[i] = d;
    }
    
    // For each vertex
    for (unsigned i = 0; i < count; ++i) {
      
      // Compute weighted transform
      glm::dualquat d = glm::normalize(
        weights[i].x * duals[weights_bones[i].x] +
        weights[i].y * duals[weights_bones[i].y] +
        weights[i].z * duals[weights_bones[i].z] +
        weights[i].w * duals[weights_bones[i].w]
      );
      
      // Transform vertex and normal
      vertices[i] = d * original_vertices[i];
      normals[i] = d.real * original_normals[i];
      
    }
  }
  
  // Apply implicit skinning to update vertices and normals
  void skin_implicit() {
    
    // Initial guess using dual quaternion skinning
    skin_dualquat();
    
    // Surface tracking
    betas.assign(count, 0.0f);
    // TODO tweak iteration count
    for (unsigned n = 0; n < 10; ++n) {
    
      // Vertex projection
      for (unsigned i = 0; i < count; ++i) {
      
        // Ignore if converged
        if (betas[i] > 0.0f)
          continue;
      
        // Compute distance
        glm::vec3 gradient;
        offsets[i] = distance(i, gradient);
        
        // Check if vertex has caused self-intersection (angle between gradients is larger than 55 degrees)
        if (n > 0 && glm::dot(gradient, previous_gradients[i]) < 0.5736f) {
          betas[i] = 1.0f;
          continue;
        }
        previous_gradients[i] = gradient;
        
        // Check if vertex has reached its isovalue
        float delta = original_offsets[i] - offsets[i];
        if (std::abs(delta) < 0.0001f) {
          betas[i] = 0.001f;
          continue;
        }
        
        // Move vertex along gradient
        vertices[i] += 0.35f * delta * gradient;
        offsets[i] += 0.35f * delta; // TODO correct?
      }
      
      // Tangential relaxation
      for (unsigned i = 0; i < count; ++i) {
        
        // Ignore if converged
        if (betas[i] > 0.0f)
          continue;
        
        // Compute relaxation strength
        float mu = std::abs(offsets[i] - original_offsets[i]) - 1.0f;
        mu = std::max(0.0f, 1.0f - mu * mu * mu * mu);
        
        // Compute expected barycenter
        // TODO project on tangential plane?
        glm::vec3 barycenter;
        for (unsigned j = 0; j < neighbors[i].size(); ++j)
          barycenter += neighbors_weights[i][j] * vertices[neighbors[i][j]];
        
        // Keep calm and relax
        vertices[i] = (1.0f - mu) * vertices[i] + mu * barycenter;
      }
    }
    
    // Diffuse coefficients
    // TODO tweak that
    for (unsigned n = 0; n < 3; ++n) {
      for (unsigned i = 0; i < count; ++i) {
        float sum = betas[i];
        for (unsigned j = 0; j < neighbors[i].size(); ++j)
          sum += betas[neighbors[i][j]];
        averages[i] = sum / (neighbors[i].size() + 1);
      }
      std::swap(betas, averages);
    }
    
    // Compute centroids
    for (unsigned i = 0; i < count; ++i) {
      glm::vec3 centroid;
      for (unsigned j = 0; j < neighbors[i].size(); ++j)
        centroid += vertices[neighbors[i][j]];
      centroids[i] = centroid / (float)neighbors[i].size();
    }
    
    // Laplacian smoothing
    for (unsigned i = 0; i < count; ++i)
      vertices[i] = (1.0f - betas[i]) * vertices[i] + betas[i] * centroids[i];
    
  }
  
  // Distance reparametrization (0 outside, 1 inside)
  float reparam(float d) const {
    if (d <= -threshold)
      return 1.0f;
    if (d >= threshold)
      return 0.0f;
    d /= threshold;
    return
      (-3.0f / 16.0f) * d * d * d * d * d +
      (5.0f / 8.0f) * d * d * d +
      (-15.0f / 16.0f) * d +
      0.5f;
  }
  
  // Compute signed distance function for specified bone and point
  float distance(unsigned bone, glm::vec3 const & point, glm::vec3 & gradient) const {
    // TODO precompute matrix inverses
    glm::mat4 inverse = glm::inverse(transforms[bone]);
    
    // Transform point into local bone space
    float result = sdfs[bone].apply(glm::vec3(inverse * glm::vec4(point, 1.0f)), gradient);
    
    // Transform back and reparametrize
    gradient = -glm::vec3(transforms[bone] * glm::vec4(gradient, 0.0f));
    return reparam(result);
  }
  
  // Compute signed distance function for specified vertex
  float distance(unsigned index, glm::vec3 & gradient) const {
    glm::vec3 g;
    float d, result = -1.0f / 0.0f;
    
    // First bone
    if (weights[index].x > 0.0f) {
      d = distance(weights_bones[index].x, vertices[index], g);
      if (d > result) {
        result = d;
        gradient = g;
      }
    }
    
    // Second bone
    if (weights[index].y > 0.0f) {
      d = distance(weights_bones[index].y, vertices[index], g);
      if (d > result) {
        result = d;
        gradient = g;
      }
    }
    
    // Third bone
    if (weights[index].z > 0.0f) {
      d = distance(weights_bones[index].z, vertices[index], g);
      if (d > result) {
        result = d;
        gradient = g;
      }
    }
    
    // Fourth bone
    if (weights[index].w > 0.0f) {
      d = distance(weights_bones[index].w, vertices[index], g);
      if (d > result) {
        result = d;
        gradient = g;
      }
    }
    
    return result;
  }
  
};

#endif
