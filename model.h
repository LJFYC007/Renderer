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

	Model(std::vector<shared_ptr<Primitive>>& _world, std::string const& filename) : world(_world)
	{
		tinygltf::Model gltfModel;
		tinygltf::TinyGLTF loader;
		std::string err;
		std::string warn;

		bool ret = loader.LoadASCIIFromFile(&gltfModel, &err, &warn, filename);
		if (!warn.empty()) std::cout << "Warn: " << warn << std::endl;
		if (!err.empty()) std::cerr << "Err: " << err << std::endl;
		if (!ret) std::cerr << "Failed to load glTF: " << filename << std::endl;

		model = std::move(gltfModel);
		std::cerr << "Loaded glTF file: " << filename << std::endl;

		if (model.defaultScene >= 0) {
			tinygltf::Scene& scene = model.scenes[model.defaultScene];
			for (size_t i = 0; i < scene.nodes.size(); ++i)
				ProcessNode(model.nodes[scene.nodes[i]], Transform());
		}
	}

	Transform GetNodeLocalTransform(tinygltf::Node const& node) {
		Transform translate = node.translation.empty() ? Transform() : Transform::Translate(vec3(node.translation[0], node.translation[1], node.translation[2]));
		Transform rotation = node.rotation.empty() ? Transform() : Transform::Rotate(node.rotation[0], node.rotation[1], node.rotation[2], node.rotation[3]);
		Transform scale = node.scale.empty() ? Transform() : Transform::Scale(node.scale[0], node.scale[1], node.scale[2]);
		return translate * rotation * scale;
	}

	void ProcessNode(tinygltf::Node& node, const Transform& transform)
	{
		Transform currentTransform = transform * GetNodeLocalTransform(node);
		if (node.mesh >= 0) 
			ProcessMesh(model.meshes[node.mesh], currentTransform);
		for (int i = 0; i < node.children.size(); ++i)
			ProcessNode(model.nodes[node.children[i]], currentTransform);
	}

	void ProcessMesh(const tinygltf::Mesh& mesh, const Transform& transform) {
		for (const auto& primitive : mesh.primitives) {
			bool uvExists = primitive.attributes.find("TEXCOORD_0") != primitive.attributes.end();
			bool normalExists = primitive.attributes.find("NORMAL") != primitive.attributes.end();

			const tinygltf::Accessor& indexAccessor = model.accessors[primitive.indices];
			const tinygltf::BufferView& bufferView = model.bufferViews[indexAccessor.bufferView];
			const tinygltf::Buffer& buffer = model.buffers[bufferView.buffer];

			std::vector<int> vertexIndices(indexAccessor.count);
			const void* dataPtr = &(buffer.data[bufferView.byteOffset + indexAccessor.byteOffset]);
			for (size_t i = 0; i < indexAccessor.count; i ++) {
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
				const float* dataPtr = reinterpret_cast<const float*>(&(buffer.data[bufferView.byteOffset + accessor.byteOffset]));
				for (size_t i = 0; i < accessor.count; i++) {
					if (attrName == "POSITION") vertices[i].p = vec3(dataPtr[i * accessor.type], dataPtr[i * accessor.type + 1], dataPtr[i * accessor.type + 2]);
					else if (attrName == "NORMAL") vertices[i].n = vec3(dataPtr[i * accessor.type], dataPtr[i * accessor.type + 1], dataPtr[i * accessor.type + 2]);
					else if (attrName == "TEXCOORD_0") vertices[i].uv = vec2(dataPtr[i * accessor.type], dataPtr[i * accessor.type + 1]);
					else if (attrName == "TEXCOORD_1");
					else std::cerr << "Unknown attribute: " << attrName << std::endl;
				}
			}

			TriangleMesh triangleMesh(transform, vertexIndices, vertices, uvExists, normalExists);
			meshes.push_back(triangleMesh);		
			shared_ptr<Material> material = make_shared<DiffuseMaterial>(make_shared<SpectrumConstantTexture>(vec3(1)));
			for (int i = 0; i < meshes.back().nTriangles; ++i)
				world.emplace_back(make_shared<SimplePrimitive>(make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, i), material));
		}
	}

private:
	tinygltf::Model model;
	std::vector<shared_ptr<Primitive>>& world;
};


/*
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <assimp/DefaultLogger.hpp>
#include <vector>
#include <string>
#include <memory>
#include <iostream>

using std::shared_ptr;
using std::make_shared;
using std::vector;

class MyErrorStream : public Assimp::LogStream {
public:
	void write(const char* message) override {
		std::cerr << message << std::endl;
	}
};

class Model
{
public:
	Model(std::vector<shared_ptr<Primitive>>& _world, std::string const& path): world(_world)
	{
		/*
		Assimp::DefaultLogger::create("", Assimp::Logger::VERBOSE);
		Assimp::Logger* pLogger = Assimp::DefaultLogger::get();
		MyErrorStream myErrStream;
		pLogger->attachStream(&myErrStream, Assimp::Logger::Err | Assimp::Logger::Warn);

		Assimp::Importer importer;
		const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_CalcTangentSpace);
		if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
		{
			assert(-1);
			return;
		}
		processNode(scene->mRootNode, scene, Transform::RotateY(pi));
	}

private:
	vector<shared_ptr<ImageTexture>> texturesLoaded;
	std::vector<shared_ptr<Primitive>>& world;

	void processNode(const aiNode* node, const aiScene* scene, const Transform& parentTransform = Transform())
	{
		Transform currentTransform = parentTransform * Transform(node->mTransformation);
		for (unsigned int i = 0; i < node->mNumMeshes; ++i)
		{
			aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
			processMesh(mesh, scene, currentTransform);
		}
		for (unsigned int i = 0; i < node->mNumChildren; ++i)
			processNode(node->mChildren[i], scene, currentTransform);
	}

	aiTextureType getTextureType(const std::string& type) {
		if (type == "DIFFUSE")
			return aiTextureType_DIFFUSE;
		else if (type == "NORMALS")
			return aiTextureType_NORMALS;
		else if (type == "HEIGHT")
			return aiTextureType_HEIGHT;
		else if (type == "ROUGHNESS")
			return aiTextureType_DIFFUSE_ROUGHNESS;
	}

	shared_ptr<SpectrumTexture> loadTexture(const std::string& type, const aiMaterial* mat) {
		aiString str;
		if (mat->GetTexture(getTextureType(type), 0, &str) != AI_SUCCESS)
		{
			aiColor3D color(0.0, 0.0, 0.0);
			if (type == "DIFFUSE") mat->Get(AI_MATKEY_COLOR_DIFFUSE, color);
			else if (type == "NORMALS") mat->Get(AI_MATKEY_COLOR_DIFFUSE, color);
			else if (type == "ROUGHNESS") mat->Get(AI_MATKEY_COLOR_DIFFUSE_ROUGHNESS, color);
			return make_shared<SpectrumConstantTexture>(UVMapping(), vec3(color.r, color.g, color.b), 1);
		}
		std::string texturePath(str.C_Str());
		texturePath = "resources/" + texturePath;

		for (auto& loadedTexture : texturesLoaded) {
			if (std::strcmp(texturePath.c_str(), loadedTexture->GetPath().c_str()) == 0) {
				return loadedTexture;
			}
		}
		UVMapping mapping;
		auto texture = make_shared<ImageTexture>(mapping, texturePath, 1);
		texturesLoaded.emplace_back(texture);
		return texture;
	}

	void processMesh(const aiMesh* mesh, const aiScene* scene, const Transform& transform)
	{
		aiMaterial* mat = scene->mMaterials[mesh->mMaterialIndex];
		shared_ptr<SpectrumTexture> texture = loadTexture("DIFFUSE", mat);
		shared_ptr<SpectrumTexture> roughnessTexture = loadTexture("ROUGHNESS", mat);
		shared_ptr<SpectrumTexture> normalTexture = loadTexture("NORMALS", mat);

		shared_ptr<ConductorMaterial> material = make_shared<ConductorMaterial>(texture, roughnessTexture);
		material->SetNormalMap(normalTexture);
		/*
		shared_ptr<ImageTexture> heightTexture = loadTexture("HEIGHT", mat);
		material->SetBumpMap(heightTexture);

		vector<Vertex> vertices;
		vector<int> indices;
		for (unsigned int i = 0; i < mesh->mNumVertices; ++i)
		{
			Vertex vertex;

			vec3 vec;
			vec.a[0] = mesh->mVertices[i].x;
			vec.a[1] = mesh->mVertices[i].y;
			vec.a[2] = mesh->mVertices[i].z;
			vertex.p = vec;

			vec.a[0] = mesh->mNormals[i].x;
			vec.a[1] = mesh->mNormals[i].y;
			vec.a[2] = mesh->mNormals[i].z;
			vertex.n = vec;

			if (mesh->mTextureCoords[0])
			{
				vec2 Vec;
				Vec.a[0] = mesh->mTextureCoords[0][i].x;
				Vec.a[1] = mesh->mTextureCoords[0][i].y;
				vertex.uv = Vec;
			}
			else vertex.uv = vec2(0.0f);

			vertices.emplace_back(vertex);
		}

		for (unsigned int i = 0; i < mesh->mNumFaces; ++i)
		{
			aiFace face = mesh->mFaces[i];
			for (unsigned int j = 0; j < face.mNumIndices; ++j)
				indices.emplace_back(face.mIndices[j]);
		}

		meshes.emplace_back(transform, indices, vertices, true, true);
		for (int i = 0; i < meshes.back().nTriangles; ++i) {
			world.emplace_back(make_shared<SimplePrimitive>(make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, i), material));
		}
	}
};
*/
