#pragma once
#include "math.h"
#include "shape.h"
#include "material.h"
#include "light.h"
#include "textures.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <vector>
#include <string>
#include <memory>

using std::shared_ptr;
using std::make_shared;
using std::vector;

class Model
{
public:
	Model(std::vector<shared_ptr<Primitive>>& World, std::string const& path) {
		Assimp::Importer importer;
		const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_CalcTangentSpace);
		if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
		{
			assert(-1);
			return;
		}
		processNode(World, scene->mRootNode, scene);
	}

private:
	vector<shared_ptr<ImageTexture>> texturesLoaded;

	void processNode(std::vector<shared_ptr<Primitive>>& World, aiNode* node, const aiScene* scene)
	{
		for (unsigned int i = 0; i < node->mNumMeshes; ++i)
		{
			aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
			processMesh(World, mesh, scene);
		}
		for (unsigned int i = 0; i < node->mNumChildren; ++i)
			processNode(World, node->mChildren[i], scene);
	}

	aiTextureType getTextureType(const std::string& type) {
		if (type == "DIFFUSE") {
			return aiTextureType_DIFFUSE;
		}
		else if (type == "NORMALS") {
			return aiTextureType_NORMALS;
		}
	}

	shared_ptr<ImageTexture> loadTexture(std::string type, aiMaterial* mat) {
		aiString str;
		mat->GetTexture(getTextureType(type), 0, &str);
		std::string texturePath = "resources/" + std::string(str.C_Str());
		for (auto& loadedTexture : texturesLoaded) {
			if (std::strcmp(texturePath.c_str(), loadedTexture->GetPath().c_str()) == 0) {
				return loadedTexture;			
			}
		}
		UVMapping mapping;
		auto texture = make_shared<ImageTexture>(mapping, texturePath, 1);
		texturesLoaded.push_back(texture);
		return texture;
	}

	void processMesh(std::vector<shared_ptr<Primitive>>& World, aiMesh* mesh, const aiScene* scene)
	{
		vector<Vertex> vertices;
		vector<int> indices;
		shared_ptr<ImageTexture> texture, normalTexture;
		if (mesh->mMaterialIndex >= 0) {
			aiMaterial* mat = scene->mMaterials[mesh->mMaterialIndex];
			texture = loadTexture("DIFFUSE", mat);
			normalTexture = loadTexture("NORMALS", mat);
		}
		shared_ptr<DiffuseMaterial> material = make_shared<DiffuseMaterial>(texture);
		//material->SetNormalMap(normalTexture);

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

			vertices.push_back(vertex);
		}

		for (unsigned int i = 0; i < mesh->mNumFaces; ++i)
		{
			aiFace face = mesh->mFaces[i];
			for (unsigned int j = 0; j < face.mNumIndices; ++j)
				indices.push_back(face.mIndices[j]);
		}

		if (mesh->mMaterialIndex >= 0)
			aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];

		meshes.push_back(TriangleMesh(Transform::Translate(vec3(0, -100, 278)) * Transform::Scale(4, 4, 4), indices, vertices));
		for (int i = 0; i < meshes.back().nTriangles; ++i) {
			World.push_back(make_shared<SimplePrimitive>(make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, i), material));
		}
	}
};
