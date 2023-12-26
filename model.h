#pragma once
#include "math.h"
#include "shape.h"
#include "material.h"
#include "light.h"
#include "textures.h"

#include <tiny_gltf.h>
#include <vector>
#include <cstring>
#include <iostream>
#include <memory>
using std::make_shared;
using std::shared_ptr;

class Model {
public:
	Model() = default;

	Model(std::vector<shared_ptr<Primitive>>& _world, std::string const& filename) : world(_world) {
		tinygltf::Model gltfModel;
		tinygltf::TinyGLTF loader;
		std::string err;
		std::string warn;

		bool ret = loader.LoadBinaryFromFile(&gltfModel, &err, &warn, filename);
		if (!warn.empty()) std::cout << "Warn: " << warn << std::endl;
		if (!err.empty()) std::cerr << "Err: " << err << std::endl;
		if (!ret) std::cerr << "Failed to load glTF: " << filename << std::endl;

		model = std::move(gltfModel);
		std::cerr << "Loaded glTF file: " << filename << std::endl;

		LoadTextures();
		LoadMaterials();
		if (model.defaultScene >= 0) {
			tinygltf::Scene& scene = model.scenes[model.defaultScene];
			for (size_t i = 0; i < scene.nodes.size(); ++i)
				ProcessNode(model.nodes[scene.nodes[i]], Transform());
		}
	}

	void LoadTextures() {
		images.resize(model.textures.size());
		for (size_t i = 0; i < model.textures.size(); ++i) {
			const auto& texture = model.textures[i];
			int imageIndex = texture.source;
			assert(imageIndex > -1);
			images[i] = make_shared<tinygltf::Image>(model.images[imageIndex]);
		}
	}

	UVMapping GetUVMapping(tinygltf::ExtensionMap extensionMap) {
		const auto& extension = extensionMap.find("KHR_texture_transform");
		if (extension == extensionMap.end()) return UVMapping();
		double su = 1, sv = 1, du = 0, dv = 0;
		if (extension->second.Has("offset")) {
			const auto& offset = extension->second.Get("offset");
			du = offset.Get(0).Get<double>();
			dv = offset.Get(1).Get<double>();
		}
		if (extension->second.Has("scale")) {
			const auto& scale = extension->second.Get("scale");
			su = scale.Get(0).Get<double>();
			sv = scale.Get(1).Get<double>();
		}
		return UVMapping(su, sv, du, dv);
	}

	void LoadMaterials() {
		materials.resize(model.materials.size());
		for (size_t i = 0; i < model.materials.size(); ++i) {
			const auto& material = model.materials[i];

			const auto& pbrTexture = material.pbrMetallicRoughness;
			shared_ptr<SpectrumTexture> baseColor;
			if (pbrTexture.baseColorTexture.index != -1) {
				const auto& textureIndex = pbrTexture.baseColorTexture.index;
				const auto& extension = pbrTexture.baseColorTexture.extensions;
				baseColor = make_shared<ImageTexture>(pbrTexture.baseColorTexture.texCoord, GetUVMapping(extension), images[textureIndex]);
			}
			else {
				const auto& baseColorFactor = pbrTexture.baseColorFactor;
				baseColor = make_shared<SpectrumConstantTexture>(vec3(baseColorFactor[0], baseColorFactor[1], baseColorFactor[2]));
			}
			auto baseMaterial = make_shared<DiffuseMaterial>(baseColor);

			shared_ptr<SpectrumTexture> metallicRoughness;
			if (pbrTexture.metallicRoughnessTexture.index != -1) {
				const auto& textureIndex = pbrTexture.metallicRoughnessTexture.index;
				const auto& extension = pbrTexture.metallicRoughnessTexture.extensions;
				metallicRoughness = make_shared<ImageTexture>(pbrTexture.metallicRoughnessTexture.texCoord, GetUVMapping(extension), images[textureIndex]);
			}
			else {
				const auto& metallicFactor = pbrTexture.metallicFactor;
				const auto& roughnessFactor = pbrTexture.roughnessFactor;
				metallicRoughness = make_shared<SpectrumConstantTexture>(vec3(0, roughnessFactor, metallicFactor));
			}
			auto layerMaterial = make_shared<DielectricMaterial>(metallicRoughness);

			materials[i] = make_shared<MixMaterial>(baseMaterial, layerMaterial);

			if (material.normalTexture.index > -1) {
				const auto& textureIndex = material.normalTexture.index;
				const auto& extension = material.normalTexture.extensions;
				materials[i]->SetNormalMap(make_shared<ImageTexture>(material.normalTexture.texCoord, GetUVMapping(extension), images[textureIndex]));
			}
		}
	}

	Transform GetNodeLocalTransform(tinygltf::Node const& node) {
		Transform translate = node.translation.empty() ? Transform() : Transform::Translate(vec3(node.translation[0], node.translation[1], node.translation[2]));
		Transform rotation = node.rotation.empty() ? Transform() : Transform::Rotate(node.rotation[0], node.rotation[1], node.rotation[2], node.rotation[3]);
		Transform scale = node.scale.empty() ? Transform() : Transform::Scale(node.scale[0], node.scale[1], node.scale[2]);
		return translate * rotation * scale;
	}

	void ProcessNode(tinygltf::Node& node, const Transform& transform) {
		Transform currentTransform = transform * GetNodeLocalTransform(node);
		if (node.mesh >= 0)
			ProcessMesh(model.meshes[node.mesh], currentTransform);
		for (int i = 0; i < node.children.size(); ++i)
			ProcessNode(model.nodes[node.children[i]], currentTransform);
	}

	size_t GetComponentSize(int componentType) {
		switch (componentType) {
		case TINYGLTF_COMPONENT_TYPE_BYTE:           return sizeof(int8_t);
		case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:  return sizeof(uint8_t);
		case TINYGLTF_COMPONENT_TYPE_SHORT:          return sizeof(int16_t);
		case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT: return sizeof(uint16_t);
		case TINYGLTF_COMPONENT_TYPE_INT:            return sizeof(int32_t);
		case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:   return sizeof(uint32_t);
		case TINYGLTF_COMPONENT_TYPE_FLOAT:          return sizeof(float);
		default:                                     return 0;
		}
	}

	size_t GetTypeCount(int type) {
		switch (type) {
		case TINYGLTF_TYPE_SCALAR: return 1;
		case TINYGLTF_TYPE_VEC2:   return 2;
		case TINYGLTF_TYPE_VEC3:   return 3;
		case TINYGLTF_TYPE_VEC4:   return 4;
		case TINYGLTF_TYPE_MAT2:   return 4;
		case TINYGLTF_TYPE_MAT3:   return 9;
		case TINYGLTF_TYPE_MAT4:   return 16;
		default:                   return 0;
		}
	}

	void ProcessMesh(const tinygltf::Mesh& mesh, const Transform& transform) {
		for (const auto& primitive : mesh.primitives) {
			bool uvExists = primitive.attributes.find("TEXCOORD_0") != primitive.attributes.end();
			bool UVExists = primitive.attributes.find("TEXCOORD_1") != primitive.attributes.end();
			bool normalExists = primitive.attributes.find("NORMAL") != primitive.attributes.end();

			const tinygltf::Accessor& indexAccessor = model.accessors[primitive.indices];
			const tinygltf::BufferView& bufferView = model.bufferViews[indexAccessor.bufferView];
			const tinygltf::Buffer& buffer = model.buffers[bufferView.buffer];

			std::vector<int> vertexIndices(indexAccessor.count);
			const void* dataPtr = &(buffer.data[bufferView.byteOffset + indexAccessor.byteOffset]);
			for (size_t i = 0; i < indexAccessor.count; i++) {
				unsigned int a;
				switch (indexAccessor.componentType) {
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
					a = static_cast<unsigned int>(reinterpret_cast<const uint8_t*>(dataPtr)[i]);
					break;
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
					a = static_cast<unsigned int>(reinterpret_cast<const uint16_t*>(dataPtr)[i]);
					break;
				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
					a = reinterpret_cast<const uint32_t*>(dataPtr)[i];
					break;
				default:
					std::cerr << "Unknown component type." << std::endl;
					return;
				}
				vertexIndices[i] = a;
			}

			const int nVertices = model.accessors[primitive.attributes.find("POSITION")->second].count;
			std::vector<Vertex> vertices(nVertices);

			for (auto& attribute : primitive.attributes) {
				const std::string& attrName = attribute.first;
				int accessorIndex = attribute.second;
				const tinygltf::Accessor& accessor = model.accessors[accessorIndex];
				const tinygltf::BufferView& bufferView = model.bufferViews[accessor.bufferView];
				const tinygltf::Buffer& buffer = model.buffers[bufferView.buffer];

				size_t componentSize = GetComponentSize(accessor.componentType);
				size_t typeCount = GetTypeCount(accessor.type);
				size_t stride = bufferView.byteStride ? bufferView.byteStride : componentSize * typeCount;

				const unsigned char* basePtr = &(buffer.data[bufferView.byteOffset + accessor.byteOffset]);
				for (size_t i = 0; i < accessor.count; i++) {
					const float* dataPtr = reinterpret_cast<const float*>(basePtr + i * stride);
					if (attrName == "POSITION") vertices[i].p = vec3(dataPtr[0], dataPtr[1], dataPtr[2]);
					else if (attrName == "NORMAL") vertices[i].n = vec3(dataPtr[0], dataPtr[1], dataPtr[2]);
					else if (attrName == "TEXCOORD_0") vertices[i].uv = vec2(dataPtr[0], dataPtr[1]);
					else if (attrName == "TEXCOORD_1") vertices[i].UV = vec2(dataPtr[0], dataPtr[1]);
					else std::cerr << "Unknown attribute: " << attrName << std::endl;
				}
			}

			TriangleMesh triangleMesh(transform, vertexIndices, vertices, uvExists, UVExists, normalExists);
			meshes.push_back(triangleMesh);

			shared_ptr<Material> material;
			if (primitive.material == -1) {
				shared_ptr<SpectrumTexture> baseColor = make_shared<SpectrumConstantTexture>(vec3(1.0));
				material = make_shared<DiffuseMaterial>(baseColor);
			}
			else material = materials[primitive.material];
			for (int i = 0; i < meshes.back().nTriangles; ++i)
				world.emplace_back(make_shared<SimplePrimitive>(make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, i), material));
		}
	}

private:
	tinygltf::Model model;
	std::vector<shared_ptr<Primitive>>& world;
	std::vector<shared_ptr<tinygltf::Image>> images;
	std::vector<shared_ptr<Material>> materials;
};
