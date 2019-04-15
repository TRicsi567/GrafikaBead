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
// egy parametrikus felület (u,v) paraméterértékekhez tartozó pontjának
// kiszámítását végzõ függvény
//
glm::vec3	CMyApp::GetUV(float u, float v)
{
	// Henger
	// R sugarú, z, tengelyû, h magasságú
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
	//2 élhosszúságú origó középpontú oktaéder 6 db csúcsa
	octa_vert.push_back(glm::vec3( 2,  0,  0));
	octa_vert.push_back(glm::vec3(-2,  0,  0));
	octa_vert.push_back(glm::vec3( 0,  0,  2));
	octa_vert.push_back(glm::vec3( 0,  0, -2));
	octa_vert.push_back(glm::vec3( 0,  2,  0));
	octa_vert.push_back(glm::vec3( 0, -2,  0));


	// törlési szín legyen kékes
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);

	glEnable(GL_CULL_FACE); // kapcsoljuk be a hatrafele nezo lapok eldobasat
	glEnable(GL_DEPTH_TEST); // mélységi teszt bekapcsolása (takarás)
	glCullFace(GL_BACK); // GL_BACK: a kamerától "elfelé" nézõ lapok, GL_FRONT: a kamera felé nézõ lapok

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

	
	// NxM darab négyszöggel közelítjük a parametrikus felületünket => (N+1)x(M+1) pontban kell kiértékelni
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

	// indexpuffer adatai: NxM négyszög = 2xNxM háromszög = háromszöglista esetén 3x2xNxM index
    GLushort indices[3*2*(N)*(M)];
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			// minden négyszögre csináljunk kettõ háromszöget, amelyek a következõ 
			// (i,j) indexeknél született (u_i, v_i) paraméterértékekhez tartozó
			// pontokat kötik össze:
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
			// - az (i,j)-hez tartózó 1D-s index a VBO-ban: i+j*(N+1)
			// - az (i,j)-hez tartózó 1D-s index az IB-ben: i*6+j*6*(N+1) 
			//		(mert minden négyszöghöz 2db háromszög = 6 index tartozik)
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
	// a frissen generált VAO beallitasa aktívnak
	glBindVertexArray(m_vaoID[0]);
	
	// hozzunk létre egy új VBO erõforrás nevet
	glGenBuffers(2, m_vboID); 

	glBindBuffer(GL_ARRAY_BUFFER, m_vboID[0]); // tegyük "aktívvá" a létrehozott VBO-t
	glBufferData(GL_ARRAY_BUFFER, sizeof(vert), vert, GL_STATIC_DRAW);

	glEnableVertexAttribArray(0); // ez lesz majd a pozíció
	glVertexAttribPointer(
		0,				// a VB-ben található adatok közül a 0. "indexû" attribútumait állítjuk be
		3,				// komponens szam
		GL_FLOAT,		// adatok tipusa
		GL_FALSE,		// normalizalt legyen-e
		sizeof(Vertex),	// stride (0=egymas utan)
		0				// a 0. indexû attribútum hol kezdõdik a sizeof(Vertex)-nyi területen belül
	); 

	// a második attribútumhoz pedig a VBO-ban sizeof(Vertex) ugrás után sizeof(glm::vec3)-nyit menve újabb 3 float adatot találunk (szín)
	
	glEnableVertexAttribArray(1); // ez lesz majd a szín
	glVertexAttribPointer(
		1,
		3, 
		GL_FLOAT,
		GL_FALSE,
		sizeof(Vertex),
		(void*)(sizeof(glm::vec3)) );
	
	// index puffer létrehozása
	glGenBuffers(1, &m_ibID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
	

	glBindVertexArray(0); // feltöltüttük a VAO-t, kapcsoljuk le
	glBindBuffer(GL_ARRAY_BUFFER, 0); // feltöltöttük a VBO-t is, ezt is vegyük le
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); // feltöltöttük a VBO-t is, ezt is vegyük le

	///////////////////////KÖR///////////////////////////

	glBindVertexArray(m_vaoID[1]);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboID[1]);
	glBufferData(GL_ARRAY_BUFFER, circle_vert.size() * sizeof(Vertex), &circle_vert[0], GL_STATIC_DRAW);
	
	glEnableVertexAttribArray(0); // ez lesz majd a pozíció
	glVertexAttribPointer(
		0,				// a VB-ben található adatok közül a 0. "indexû" attribútumait állítjuk be
		3,				// komponens szam
		GL_FLOAT,		// adatok tipusa
		GL_FALSE,		// normalizalt legyen-e
		sizeof(Vertex),	// stride (0=egymas utan)
		0				// a 0. indexû attribútum hol kezdõdik a sizeof(Vertex)-nyi területen belül
	);

	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),	(void*)(sizeof(glm::vec3)) );



	glBindVertexArray(0); // feltöltüttük a VAO-t, kapcsoljuk le
	glBindBuffer(GL_ARRAY_BUFFER, 0); // feltöltöttük a VBO-t is, ezt is vegyük le
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); // feltöltöttük a VBO-t is, ezt is vegyük le


	//
	// shaderek betöltése
	//
	GLuint vs_ID = loadShader(GL_VERTEX_SHADER,		"myVert.vert");
	GLuint fs_ID = loadShader(GL_FRAGMENT_SHADER,	"myFrag.frag");

	// a shadereket tároló program létrehozása
	m_programID = glCreateProgram();

	// adjuk hozzá a programhoz a shadereket
	glAttachShader(m_programID, vs_ID);
	glAttachShader(m_programID, fs_ID);

	// VAO-beli attribútumok hozzárendelése a shader változókhoz
	// FONTOS: linkelés elõtt kell ezt megtenni!
	glBindAttribLocation(	m_programID,	// shader azonosítója, amibõl egy változóhoz szeretnénk hozzárendelést csinálni
							0,				// a VAO-beli azonosító index
							"vs_in_pos");	// a shader-beli változónév
	glBindAttribLocation( m_programID, 1, "vs_in_col");

	// illesszük össze a shadereket (kimenõ-bemenõ változók összerendelése stb.)
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

		std::cout << "[app.Init()] Sáder Huba panasza: " << aSzoveg << std::endl;

		delete aSzoveg;
	}

	// mar nincs ezekre szukseg
	glDeleteShader( vs_ID );
	glDeleteShader( fs_ID );

	//
	// egyéb inicializálás
	//

	// vetítési mátrix létrehozása
	m_matProj = glm::perspective( 45.0f, 640/480.0f, 1.0f, 1000.0f );

	// shader-beli transzformációs mátrixok címének lekérdezése
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
	// nézeti transzformáció beállítása
	float t = SDL_GetTicks()/1000.0f;
	m_matView = glm::lookAt(glm::vec3( look_at_from_x,  look_at_from_y,  look_at_from_z),		// honnan nézzük a színteret
							glm::vec3( 0,  0,  0),		// a színtér melyik pontját nézzük
							glm::vec3( 0,  1,  0));		// felfelé mutató irány a világban
}


void CMyApp::Render()
{
	// töröljük a frampuffert (GL_COLOR_BUFFER_BIT) és a mélységi Z puffert (GL_DEPTH_BUFFER_BIT)
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// shader bekapcsolasa
	glUseProgram( m_programID );

	// shader parameterek beállítása
	/*

	GLM transzformációs mátrixokra példák:
		glm::rotate<float>( szög, glm::vec3(tengely_x, tengely_y, tengely_z) ) <- tengely_{xyz} körüli elforgatás
		glm::translate<float>( glm::vec3(eltol_x, eltol_y, eltol_z) ) <- eltolás
		glm::scale<float>( glm::vec3(s_x, s_y, s_z) ) <- léptékezés

	*/
	m_matWorld = glm::mat4(1.0f);

	float now = SDL_GetTicks() * (360.f / 12000.f);


	glm::mat4 mvp = m_matProj * m_matView * m_matWorld;

	// majd küldjük át a megfelelõ mátrixot!
	glUniformMatrix4fv( m_loc_mvp,// erre a helyre töltsünk át adatot
						1,			// egy darab mátrixot
						GL_FALSE,	// NEM transzponálva
						&(mvp[0][0]) ); // innen olvasva a 16 x sizeof(float)-nyi adatot

	
	
	// kapcsoljuk be a VAO-t (a VBO jön vele együtt)
	glBindVertexArray(m_vaoID[0]);

	// kirajzolás henger
	for(unsigned int i = 0; i < octa_vert.size(); ++i)
	{
		mvp = m_matProj * m_matView * m_matWorld  * glm::translate<float>(compute_pos(now)) * glm::translate<float>(octa_vert[i]);
		glUniformMatrix4fv(m_loc_mvp, 1, GL_FALSE, &(mvp[0][0]));
		glDrawElements(GL_TRIANGLES, 3 * 2 * (N)*(M), GL_UNSIGNED_SHORT, 0);
	}
	/*
	
	glDrawElements(	GL_TRIANGLES,		// primitív típus
					3*2*(N)*(M),		// hany csucspontot hasznalunk a kirajzolashoz
					GL_UNSIGNED_SHORT,	// indexek tipusa
					0);					// indexek cime
	*/

	glBindVertexArray(m_vaoID[1]);

	// kirajzolás teteje
	for (unsigned int i = 0; i < octa_vert.size(); ++i)
	{
		mvp = m_matProj * m_matView * m_matWorld * glm::translate<float>(compute_pos(now)) * glm::translate<float>(octa_vert[i]) * glm::translate<float>(glm::vec3(0, 3, 0)) ;
		glUniformMatrix4fv(m_loc_mvp, 1, GL_FALSE, &(mvp[0][0]));
		glDrawArrays(GL_TRIANGLE_FAN, 0, 362);
	}

	// kirajzolás alja
	for (unsigned int i = 0; i < octa_vert.size(); ++i)
	{
		mvp = m_matProj * m_matView * m_matWorld * glm::translate<float>(compute_pos(now)) * glm::translate<float>(octa_vert[i]) * glm::rotate<float>((float)M_PI, glm::vec3(1, 0, 0));
		glUniformMatrix4fv(m_loc_mvp, 1, GL_FALSE, &(mvp[0][0]));
		glDrawArrays(GL_TRIANGLE_FAN, 0, 362);
	}

	


	// mozgatás az y = 0.02*x^2 parabola mentén -10 és 5 között az XY síkon oda vissza 12 sec alatt
	// Ívhossz ~ 15.3
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

// a két paraméterbe az új ablakméret szélessége (_w) és magassága (_h) található
void CMyApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h);

	m_matProj = glm::perspective(  45.0f,		// 90 fokos nyilasszog
									_w/(float)_h,	// ablakmereteknek megfelelo nezeti arany
									0.01f,			// kozeli vagosik
									100.0f);		// tavoli vagosik
}