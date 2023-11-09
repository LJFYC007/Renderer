#pragma once
#include "math.h"
#include "shape.h"
#include "material.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <vector>
#include <string>
#include <memory>

using std::shared_ptr;
using std::make_shared;
using std::string;
using std::vector;

class Model
{
public:
	bool gammaCorrection;

	Model(string const& path, shared_ptr<material> mat) { loadModel(path, mat); }

private:
	void loadModel(string const& path, shared_ptr<material> mat)
	{
		Assimp::Importer importer;
		const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs);
		if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
		{
			exit(-1);
			return;
		}
		processNode(scene->mRootNode, scene, mat);
	}

	void processNode(aiNode* node, const aiScene* scene, shared_ptr<material> mat)
	{
		for (unsigned int i = 0; i < node->mNumMeshes; ++i)
		{
			aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
			meshes.push_back(processMesh(mesh, scene, mat));
		}
		for (unsigned int i = 0; i < node->mNumChildren; ++i)
			processNode(node->mChildren[i], scene, mat);
	}

	TriangleMesh processMesh(aiMesh* mesh, const aiScene* scene, shared_ptr<material> mat)
	{
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

		return TriangleMesh(SquareMatrix<4>(), indices, vertices, mat);
	}
};
