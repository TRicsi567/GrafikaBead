#include "MyApp.h"
#include "GLUtils.hpp"
#include <iostream>
#include <math.h>
#include <imgui/imgui.h>


CMyApp::CMyApp(void)
{
	m_vaoID[0] = 0;
	m_vaoID[1] = 0;

	m_vboID[0] = 0;
	m_vboID[1] = 0;

	m_ibID = 0;

	m_programID = 0;

	// m_camera.SetView(glm::vec3(5, 5, 5), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
}


CMyApp::~CMyApp(void)
{
}

//
// egy parametrikus fel�let (u,v) param�ter�rt�kekhez tartoz� pontj�nak
// kisz�m�t�s�t v�gz� f�ggv�ny
//
glm::vec3	CMyApp::GetUV(float u, float v)
{
	// Henger
	// R sugar�, z, tengely�, h magass�g�
	// x(u,v) = R * cos(2*pi*u)
	// y(u,v) = R *sinc(2*pi*u)
	// z(u,v) = h * v

	float r = 1;
	float h = 3;
	
	float x = r * cosf(2 * (float)M_PI* u);
	float y = r * sinf(2 * (float)M_PI * u);
	float z = h * v;
	
	/*
	float x = r * (1 - v) * cosf(2 * (float)M_PI*u);
	float y = r * (1 - v) * sinf(2 * (float)M_PI*u);
	float z = h * v;
	*/
	return glm::vec3(x, z, -y);
}

bool CMyApp::Init()
{
	//2 �lhossz�s�g� orig� k�z�ppont� okta�der 6 db cs�csa
	octa_vert.push_back(glm::vec3( 2,  0,  0));
	octa_vert.push_back(glm::vec3(-2,  0,  0));
	octa_vert.push_back(glm::vec3( 0,  0,  2));
	octa_vert.push_back(glm::vec3( 0,  0, -2));
	octa_vert.push_back(glm::vec3( 0,  2,  0));
	octa_vert.push_back(glm::vec3( 0, -2,  0));


	// t�rl�si sz�n legyen k�kes
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);

	glEnable(GL_CULL_FACE); // kapcsoljuk be a hatrafele nezo lapok eldobasat
	glEnable(GL_DEPTH_TEST); // m�lys�gi teszt bekapcsol�sa (takar�s)
	glCullFace(GL_BACK); // GL_BACK: a kamer�t�l "elfel�" n�z� lapok, GL_FRONT: a kamera fel� n�z� lapok

	//
	// geometria letrehozasa
	//

	std::vector<Vertex> circle_vert;

	circle_vert.push_back(Vertex{ glm::vec3(0, 0, 0), glm::vec3(0.8f, 0.8f, 0.8f) });
	float radius = 1.f;
	float deg_to_rad = (float)M_PI / 180.f;

	for (int angle = 1; angle < 361; ++angle)
	{
		const glm::vec3& center = circle_vert[0].p;
		Vertex foo;
		foo.p = glm::vec3((center.x + cosf(angle * deg_to_rad) * radius), (center.z), -(center.y + sinf(angle * deg_to_rad) * radius));
		foo.c = glm::vec3(0.6f, 0.2f, 0.5f);
		circle_vert.push_back(foo);
	}

	circle_vert.push_back(circle_vert.at(1));
	/*
	for (unsigned int i = 0; i < circle_vert.size(); ++i) {
		std::cout << circle_vert[i].c.x << " " << circle_vert[i].c.y << " " << circle_vert[i].c.z << std::endl;
	}
	*/

	
	// NxM darab n�gysz�ggel k�zel�tj�k a parametrikus fel�let�nket => (N+1)x(M+1) pontban kell ki�rt�kelni
	Vertex vert[(N+1)*(M+1)];
	for (int i = 0; i <= N; ++i)
	{
		for (int j = 0; j <= M; ++j)
		{
			float u = i / (float)N;
			float v = j / (float)M;

			vert[i + j * (N + 1)].p = GetUV(u, v);
			vert[i + j * (N + 1)].c = glm::normalize(vert[i + j * (N + 1)].p);
			// vert[i + j * (N + 1)].c = glm::vec3(0.6f, 0.2f, 0.5f);
		}
	}

	// indexpuffer adatai: NxM n�gysz�g = 2xNxM h�romsz�g = h�romsz�glista eset�n 3x2xNxM index
    GLushort indices[3*2*(N)*(M)];
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			// minden n�gysz�gre csin�ljunk kett� h�romsz�get, amelyek a k�vetkez� 
			// (i,j) indexekn�l sz�letett (u_i, v_i) param�ter�rt�kekhez tartoz�
			// pontokat k�tik �ssze:
			//
			//		(i,j+1)
			//		  o-----o(i+1,j+1)
			//		  |\    |			a = p(u_i, v_i)
			//		  | \   |			b = p(u_{i+1}, v_i)
			//		  |  \  |			c = p(u_i, v_{i+1})
			//		  |   \ |			d = p(u_{i+1}, v_{i+1})
			//		  |	   \|
			//	(i,j) o-----o(i+1, j)
			//
			// - az (i,j)-hez tart�z� 1D-s index a VBO-ban: i+j*(N+1)
			// - az (i,j)-hez tart�z� 1D-s index az IB-ben: i*6+j*6*(N+1) 
			//		(mert minden n�gysz�gh�z 2db h�romsz�g = 6 index tartozik)
			//
			indices[6 * i + j * 3 * 2 * (N)+0] = (i)+(j)*	(N + 1);
			indices[6 * i + j * 3 * 2 * (N)+1] = (i + 1) + (j)*	(N + 1);
			indices[6 * i + j * 3 * 2 * (N)+2] = (i)+(j + 1)*(N + 1);
			indices[6 * i + j * 3 * 2 * (N)+3] = (i + 1) + (j)*	(N + 1);
			indices[6 * i + j * 3 * 2 * (N)+4] = (i + 1) + (j + 1)*(N + 1);
			indices[6 * i + j * 3 * 2 * (N)+5] = (i)+(j + 1)*(N + 1);
		}
	}


	// 2 db VAO foglalasa
	glGenVertexArrays(2, m_vaoID);
	// a frissen gener�lt VAO beallitasa akt�vnak
	glBindVertexArray(m_vaoID[0]);
	
	// hozzunk l�tre egy �j VBO er�forr�s nevet
	glGenBuffers(2, m_vboID); 

	glBindBuffer(GL_ARRAY_BUFFER, m_vboID[0]); // tegy�k "akt�vv�" a l�trehozott VBO-t
	glBufferData(GL_ARRAY_BUFFER, sizeof(vert), vert, GL_STATIC_DRAW);

	glEnableVertexAttribArray(0); // ez lesz majd a poz�ci�
	glVertexAttribPointer(
		0,				// a VB-ben tal�lhat� adatok k�z�l a 0. "index�" attrib�tumait �ll�tjuk be
		3,				// komponens szam
		GL_FLOAT,		// adatok tipusa
		GL_FALSE,		// normalizalt legyen-e
		sizeof(Vertex),	// stride (0=egymas utan)
		0				// a 0. index� attrib�tum hol kezd�dik a sizeof(Vertex)-nyi ter�leten bel�l
	); 

	// a m�sodik attrib�tumhoz pedig a VBO-ban sizeof(Vertex) ugr�s ut�n sizeof(glm::vec3)-nyit menve �jabb 3 float adatot tal�lunk (sz�n)
	
	glEnableVertexAttribArray(1); // ez lesz majd a sz�n
	glVertexAttribPointer(
		1,
		3, 
		GL_FLOAT,
		GL_FALSE,
		sizeof(Vertex),
		(void*)(sizeof(glm::vec3)) );
	
	// index puffer l�trehoz�sa
	glGenBuffers(1, &m_ibID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
	

	glBindVertexArray(0); // felt�lt�tt�k a VAO-t, kapcsoljuk le
	glBindBuffer(GL_ARRAY_BUFFER, 0); // felt�lt�tt�k a VBO-t is, ezt is vegy�k le
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); // felt�lt�tt�k a VBO-t is, ezt is vegy�k le

	///////////////////////K�R///////////////////////////

	glBindVertexArray(m_vaoID[1]);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboID[1]);
	glBufferData(GL_ARRAY_BUFFER, circle_vert.size() * sizeof(Vertex), &circle_vert[0], GL_STATIC_DRAW);
	
	glEnableVertexAttribArray(0); // ez lesz majd a poz�ci�
	glVertexAttribPointer(
		0,				// a VB-ben tal�lhat� adatok k�z�l a 0. "index�" attrib�tumait �ll�tjuk be
		3,				// komponens szam
		GL_FLOAT,		// adatok tipusa
		GL_FALSE,		// normalizalt legyen-e
		sizeof(Vertex),	// stride (0=egymas utan)
		0				// a 0. index� attrib�tum hol kezd�dik a sizeof(Vertex)-nyi ter�leten bel�l
	);

	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),	(void*)(sizeof(glm::vec3)) );



	glBindVertexArray(0); // felt�lt�tt�k a VAO-t, kapcsoljuk le
	glBindBuffer(GL_ARRAY_BUFFER, 0); // felt�lt�tt�k a VBO-t is, ezt is vegy�k le
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); // felt�lt�tt�k a VBO-t is, ezt is vegy�k le


	//
	// shaderek bet�lt�se
	//
	GLuint vs_ID = loadShader(GL_VERTEX_SHADER,		"myVert.vert");
	GLuint fs_ID = loadShader(GL_FRAGMENT_SHADER,	"myFrag.frag");

	// a shadereket t�rol� program l�trehoz�sa
	m_programID = glCreateProgram();

	// adjuk hozz� a programhoz a shadereket
	glAttachShader(m_programID, vs_ID);
	glAttachShader(m_programID, fs_ID);

	// VAO-beli attrib�tumok hozz�rendel�se a shader v�ltoz�khoz
	// FONTOS: linkel�s el�tt kell ezt megtenni!
	glBindAttribLocation(	m_programID,	// shader azonos�t�ja, amib�l egy v�ltoz�hoz szeretn�nk hozz�rendel�st csin�lni
							0,				// a VAO-beli azonos�t� index
							"vs_in_pos");	// a shader-beli v�ltoz�n�v
	glBindAttribLocation( m_programID, 1, "vs_in_col");

	// illessz�k �ssze a shadereket (kimen�-bemen� v�ltoz�k �sszerendel�se stb.)
	glLinkProgram(m_programID);

	// linkeles ellenorzese
	GLint infoLogLength = 0, result = 0;

	glGetProgramiv(m_programID, GL_LINK_STATUS, &result);
	glGetProgramiv(m_programID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if ( GL_FALSE == result )
	{
		std::vector<char> ProgramErrorMessage( infoLogLength );
		glGetProgramInfoLog(m_programID, infoLogLength, NULL, &ProgramErrorMessage[0]);
		fprintf(stdout, "%s\n", &ProgramErrorMessage[0]);
		
		char* aSzoveg = new char[ProgramErrorMessage.size()];
		memcpy( aSzoveg, &ProgramErrorMessage[0], ProgramErrorMessage.size());

		std::cout << "[app.Init()] S�der Huba panasza: " << aSzoveg << std::endl;

		delete aSzoveg;
	}

	// mar nincs ezekre szukseg
	glDeleteShader( vs_ID );
	glDeleteShader( fs_ID );

	//
	// egy�b inicializ�l�s
	//

	// vet�t�si m�trix l�trehoz�sa
	m_matProj = glm::perspective( 45.0f, 640/480.0f, 1.0f, 1000.0f );

	// shader-beli transzform�ci�s m�trixok c�m�nek lek�rdez�se
	m_loc_mvp = glGetUniformLocation( m_programID, "MVP");

	last_time = SDL_GetTicks();

	return true;
}

void CMyApp::Clean()
{
	glDeleteBuffers(1, m_vboID);
	glDeleteBuffers(1, &m_ibID);
	glDeleteVertexArrays(1, m_vaoID);

	glDeleteProgram( m_programID );
}

void CMyApp::Update()
{
	// n�zeti transzform�ci� be�ll�t�sa
	float t = SDL_GetTicks()/1000.0f;
	m_matView = glm::lookAt(glm::vec3( look_at_from_x,  look_at_from_y,  look_at_from_z),		// honnan n�zz�k a sz�nteret
							glm::vec3( 0,  0,  0),		// a sz�nt�r melyik pontj�t n�zz�k
							glm::vec3( 0,  1,  0));		// felfel� mutat� ir�ny a vil�gban
}


void CMyApp::Render()
{
	// t�r�lj�k a frampuffert (GL_COLOR_BUFFER_BIT) �s a m�lys�gi Z puffert (GL_DEPTH_BUFFER_BIT)
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// shader bekapcsolasa
	glUseProgram( m_programID );

	// shader parameterek be�ll�t�sa
	/*

	GLM transzform�ci�s m�trixokra p�ld�k:
		glm::rotate<float>( sz�g, glm::vec3(tengely_x, tengely_y, tengely_z) ) <- tengely_{xyz} k�r�li elforgat�s
		glm::translate<float>( glm::vec3(eltol_x, eltol_y, eltol_z) ) <- eltol�s
		glm::scale<float>( glm::vec3(s_x, s_y, s_z) ) <- l�pt�kez�s

	*/
	m_matWorld = glm::mat4(1.0f);

	float now = SDL_GetTicks() * (360.f / 12000.f);


	glm::mat4 mvp = m_matProj * m_matView * m_matWorld;

	// majd k�ldj�k �t a megfelel� m�trixot!
	glUniformMatrix4fv( m_loc_mvp,// erre a helyre t�lts�nk �t adatot
						1,			// egy darab m�trixot
						GL_FALSE,	// NEM transzpon�lva
						&(mvp[0][0]) ); // innen olvasva a 16 x sizeof(float)-nyi adatot

	
	
	// kapcsoljuk be a VAO-t (a VBO j�n vele egy�tt)
	glBindVertexArray(m_vaoID[0]);

	// kirajzol�s henger
	for(unsigned int i = 0; i < octa_vert.size(); ++i)
	{
		mvp = m_matProj * m_matView * m_matWorld  * glm::translate<float>(compute_pos(now)) * glm::translate<float>(octa_vert[i]);
		glUniformMatrix4fv(m_loc_mvp, 1, GL_FALSE, &(mvp[0][0]));
		glDrawElements(GL_TRIANGLES, 3 * 2 * (N)*(M), GL_UNSIGNED_SHORT, 0);
	}
	/*
	
	glDrawElements(	GL_TRIANGLES,		// primit�v t�pus
					3*2*(N)*(M),		// hany csucspontot hasznalunk a kirajzolashoz
					GL_UNSIGNED_SHORT,	// indexek tipusa
					0);					// indexek cime
	*/

	glBindVertexArray(m_vaoID[1]);

	// kirajzol�s teteje
	for (unsigned int i = 0; i < octa_vert.size(); ++i)
	{
		mvp = m_matProj * m_matView * m_matWorld * glm::translate<float>(compute_pos(now)) * glm::translate<float>(octa_vert[i]) * glm::translate<float>(glm::vec3(0, 3, 0)) ;
		glUniformMatrix4fv(m_loc_mvp, 1, GL_FALSE, &(mvp[0][0]));
		glDrawArrays(GL_TRIANGLE_FAN, 0, 362);
	}

	// kirajzol�s alja
	for (unsigned int i = 0; i < octa_vert.size(); ++i)
	{
		mvp = m_matProj * m_matView * m_matWorld * glm::translate<float>(compute_pos(now)) * glm::translate<float>(octa_vert[i]) * glm::rotate<float>((float)M_PI, glm::vec3(1, 0, 0));
		glUniformMatrix4fv(m_loc_mvp, 1, GL_FALSE, &(mvp[0][0]));
		glDrawArrays(GL_TRIANGLE_FAN, 0, 362);
	}

	


	// mozgat�s az y = 0.02*x^2 parabola ment�n -10 �s 5 k�z�tt az XY s�kon oda vissza 12 sec alatt
	// �vhossz ~ 15.3
	/*
	int now = SDL_GetTicks();
	delta_time += (now - last_time);
	last_time = SDL_GetTicks();

	if (delta_time >= 1000)
	{
		std::cout << "eltelt 1 sec" << std::endl;
		last_time = now;
		delta_time = 0;
	}
	*/


	// VAO kikapcsolasa
	glBindVertexArray(0);

	// shader kikapcsolasa
	glUseProgram( 0 );
}

glm::vec3 CMyApp::compute_pos(const float& time)
{
	float period_time = 12.f;

	float scaler = 0.01f;
	
	float x = cos(time * ((float)M_PI/180)) * 7.5f - 2.5f;

	float y = 0.02f * pow(x, 2);
	
	return glm::vec3(x, y , 0);
}

void CMyApp::KeyboardDown(SDL_KeyboardEvent& key)
{
	double r_horiz = sqrt(pow(look_at_from_z, 2) + pow(look_at_from_x, 2));
	double r_vert = sqrt(pow(look_at_from_z, 2) + pow(look_at_from_y, 2));

	switch (key.keysym.sym)
	{
	case SDLK_w:
		look_at_from_x *= 0.9f;
		look_at_from_y *= 0.9f;
		look_at_from_z *= 0.9f;
		break;
	case SDLK_a:
		deg_horiz += 5;
		look_at_from_x = r_horiz * cos(deg_horiz * (M_PI / 180));
		look_at_from_z = r_horiz * sin(deg_horiz * (M_PI / 180));
		break;
	case SDLK_s:
		look_at_from_x *= 1.1f;
		look_at_from_y *= 1.1f;
		look_at_from_z *= 1.1f;
		break;
	case SDLK_d:
		deg_horiz -= 5;
		look_at_from_x = r_horiz * cos(deg_horiz * (M_PI / 180));
		look_at_from_z = r_horiz * sin(deg_horiz * (M_PI / 180));
		break;
	case SDLK_UP:
		deg_vert += 5;
		look_at_from_z = r_vert * cos(deg_vert * (M_PI / 180));
		look_at_from_y = r_vert * sin(deg_vert * (M_PI / 180));
		// ++look_at_from_y;
		break;
	case SDLK_DOWN:
		deg_vert -= 5;
		look_at_from_z = r_vert * cos(deg_vert * (M_PI / 180));
		look_at_from_y = r_vert * sin(deg_vert * (M_PI / 180));
		// --look_at_from_y;
		break;
	}

}

void CMyApp::KeyboardUp(SDL_KeyboardEvent& key)
{

}

void CMyApp::MouseMove(SDL_MouseMotionEvent& mouse)
{

}

void CMyApp::MouseDown(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseUp(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseWheel(SDL_MouseWheelEvent& wheel)
{
}

// a k�t param�terbe az �j ablakm�ret sz�less�ge (_w) �s magass�ga (_h) tal�lhat�
void CMyApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h);

	m_matProj = glm::perspective(  45.0f,		// 90 fokos nyilasszog
									_w/(float)_h,	// ablakmereteknek megfelelo nezeti arany
									0.01f,			// kozeli vagosik
									100.0f);		// tavoli vagosik
}