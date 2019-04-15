#pragma once

// GLEW
#include <GL/glew.h>

// SDL
#include <SDL.h>
#include <SDL_opengl.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>

#include <vector>


class CMyApp
{
public:
	CMyApp(void);
	~CMyApp(void);

	bool Init();
	void Clean();

	void Update();
	void Render();

	void KeyboardDown(SDL_KeyboardEvent&);
	void KeyboardUp(SDL_KeyboardEvent&);
	void MouseMove(SDL_MouseMotionEvent&);
	void MouseDown(SDL_MouseButtonEvent&);
	void MouseUp(SDL_MouseButtonEvent&);
	void MouseWheel(SDL_MouseWheelEvent&);
	void Resize(int, int);
protected:
	// segédfüggvények
	glm::vec3 GetUV(float u, float v);

	// shaderekhez szükséges változók
	GLuint m_programID; // shaderek programja

	// OpenGL-es dolgok
	GLuint m_vaoID[2]; // vertex array object erõforrás azonosító
	GLuint m_vboID[2]; // vertex buffer object erõforrás azonosító
	GLuint m_ibID;  // index buffer object erõforrás azonosító

	// transzformációs mátrixok
	glm::mat4 m_matWorld;
	glm::mat4 m_matView;
	glm::mat4 m_matProj;

	// mátrixok helye a shaderekben
	GLuint	m_loc_mvp; // a három mátrixunk szorzatát adjuk át a hatékonyság érdekében

	struct Vertex
	{
		glm::vec3 p;
		glm::vec3 c;
	};

	
	double look_at_from_x = 0;
	double look_at_from_y = 0;
	double look_at_from_z = 15;

	double deg_horiz = 90;
	double deg_vert = 0;

	std::vector<glm::vec3> octa_vert;
	int last_time = 0;
	int delta_time = 0;

	glm::vec3 compute_pos(const float& time);

	// NxM darab négyszöggel közelítjük a parametrikus felületünket => (N+1)x(M+1) pontban kell kiértékelni
	static const int N	= 20;
	static const int M	= 10;
};

