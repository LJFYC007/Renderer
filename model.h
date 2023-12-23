#pragma once
#include "math.h"
#include "shape.h"
#include "material.h"
#include "light.h"
#include "textures.h"

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
		*/

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
		if (type == "DIFFUSE") {
			return aiTextureType_DIFFUSE;
		}
		else if (type == "NORMALS") {
			return aiTextureType_NORMALS;
		}
		else if (type == "HEIGHT") {
			return aiTextureType_HEIGHT;
		}
	}

	shared_ptr<ImageTexture> loadTexture(const std::string& type, const aiMaterial* mat) {
		aiString str;
		if (mat->GetTexture(getTextureType(type), 0, &str) != AI_SUCCESS)
			return nullptr;
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
		shared_ptr<ImageTexture> texture = loadTexture("DIFFUSE", mat);
		shared_ptr<DiffuseMaterial> material = make_shared<DiffuseMaterial>(texture);
		shared_ptr<ImageTexture> normalTexture = loadTexture("NORMAL", mat);
		material->SetNormalMap(normalTexture);
		/*
		shared_ptr<ImageTexture> heightTexture = loadTexture("HEIGHT", mat);
		material->SetBumpMap(heightTexture);
		*/

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
