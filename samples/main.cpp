// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#define _CRT_SECURE_NO_WARNINGS
//#define IMGUI_DISABLE_OBSOLETE_FUNCTIONS 1

#if defined(_WIN32)
#include <crtdbg.h>

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <windows.h>
#endif


#include "draw.h"
#include "sample.h"
#include "settings.h"

#include "solver2d/constants.h"
#include "solver2d/math.h"
#include "solver2d/timer.h"

#include <glad/glad.h>
// Keep glad.h before glfw3.h
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(_WIN32)
static int MyAllocHook(int allocType, void* userData, size_t size, int blockType, long requestNumber, const unsigned char* filename,
					   int lineNumber)
{
	// This hook can help find leaks
	if (size == 33660)
	{
		size += 0;
	}

	return 1;
}
#endif

GLFWwindow* g_mainWindow = nullptr;
static int32_t s_selection = 0;
static Sample* s_sample = nullptr;
static Settings s_settings;
static bool s_rightMouseDown = false;
static s2Vec2 s_clickPointWS = s2Vec2_zero;
static float s_windowScale = 1.0f;
static float s_framebufferScale = 1.0f;

void glfwErrorCallback(int error, const char* description)
{
	fprintf(stderr, "GLFW error occurred. Code: %d. Description: %s\n", error, description);
}

static inline int CompareSamples(const void* a, const void* b)
{
	SampleEntry* sa = (SampleEntry*)a;
	SampleEntry* sb = (SampleEntry*)b;

	int result = strcmp(sa->category, sb->category);
	if (result == 0)
	{
		result = strcmp(sa->name, sb->name);
	}

	return result;
}

static void SortTests()
{
	qsort(g_sampleEntries, g_sampleCount, sizeof(SampleEntry), CompareSamples);
}

static void RestartTest()
{
	delete s_sample;
	s_sample = g_sampleEntries[s_settings.m_sampleIndex].createFcn(s_settings);
}

static void CreateUI(GLFWwindow* window, const char* glslVersion)
{
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	bool success;
	success = ImGui_ImplGlfw_InitForOpenGL(window, false);
	if (success == false)
	{
		printf("ImGui_ImplGlfw_InitForOpenGL failed\n");
		assert(false);
	}

	success = ImGui_ImplOpenGL3_Init(glslVersion);
	if (success == false)
	{
		printf("ImGui_ImplOpenGL3_Init failed\n");
		assert(false);
	}

	// Search for font file
	const char* fontPath1 = "data/droid_sans.ttf";
	const char* fontPath2 = "../data/droid_sans.ttf";
	const char* fontPath = nullptr;
	FILE* file1 = fopen(fontPath1, "rb");
	FILE* file2 = fopen(fontPath2, "rb");
	if (file1)
	{
		fontPath = fontPath1;
		fclose(file1);
	}

	if (file2)
	{
		fontPath = fontPath2;
		fclose(file2);
	}

	if (fontPath)
	{
		ImFontConfig fontConfig;
		fontConfig.RasterizerMultiply = s_windowScale * s_framebufferScale;
		ImGui::GetIO().Fonts->AddFontFromFileTTF(fontPath, 14.0f, &fontConfig);
	}
}

static void DestroyUI()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
}

static void ResizeWindowCallback(GLFWwindow*, int width, int height)
{
	g_camera.m_width = int(width / s_windowScale);
	g_camera.m_height = int(height / s_windowScale);
	s_settings.m_windowWidth = int(width / s_windowScale);
	s_settings.m_windowHeight = int(height / s_windowScale);
}

static void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);
	if (ImGui::GetIO().WantCaptureKeyboard)
	{
		return;
	}

	if (action == GLFW_PRESS)
	{
		switch (key)
		{
			case GLFW_KEY_ESCAPE:
				// Quit
				glfwSetWindowShouldClose(g_mainWindow, GL_TRUE);
				break;

			case GLFW_KEY_LEFT:
				g_camera.m_center.x -= 0.5f;
				break;

			case GLFW_KEY_RIGHT:
				g_camera.m_center.x += 0.5f;
				break;

			case GLFW_KEY_DOWN:
				g_camera.m_center.y -= 0.5f;
				break;

			case GLFW_KEY_UP:
				g_camera.m_center.y += 0.5f;
				break;

			case GLFW_KEY_HOME:
				g_camera.ResetView();
				break;

			case GLFW_KEY_R:
				RestartTest();
				break;

			case GLFW_KEY_O:
				s_settings.m_singleStep = true;
				break;

			case GLFW_KEY_P:
				s_settings.m_pause = !s_settings.m_pause;
				break;

			case GLFW_KEY_LEFT_BRACKET:
				// Switch to previous test
				--s_selection;
				if (s_selection < 0)
				{
					s_selection = g_sampleCount - 1;
				}
				break;

			case GLFW_KEY_RIGHT_BRACKET:
				// Switch to next test
				++s_selection;
				if (s_selection == g_sampleCount)
				{
					s_selection = 0;
				}
				break;

			case GLFW_KEY_TAB:
				g_draw.m_showUI = !g_draw.m_showUI;

			default:
				if (s_sample)
				{
					s_sample->Keyboard(key);
				}
		}
	}
	else if (action == GLFW_RELEASE)
	{
		s_sample->KeyboardUp(key);
	}
	// else GLFW_REPEAT
}

static void CharCallback(GLFWwindow* window, unsigned int c)
{
	ImGui_ImplGlfw_CharCallback(window, c);
}

static void MouseButtonCallback(GLFWwindow* window, int32_t button, int32_t action, int32_t mods)
{
	ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);

	if (ImGui::GetIO().WantCaptureMouse)
	{
		return;
	}

	double xd, yd;
	glfwGetCursorPos(g_mainWindow, &xd, &yd);
	s2Vec2 ps = {float(xd) / s_windowScale, float(yd) / s_windowScale};

	// Use the mouse to move things around.
	if (button == GLFW_MOUSE_BUTTON_1)
	{
		s2Vec2 pw = g_camera.ConvertScreenToWorld(ps);
		if (action == GLFW_PRESS)
		{
			s_sample->MouseDown(pw, button, mods);
		}

		if (action == GLFW_RELEASE)
		{
			s_sample->MouseUp(pw, button);
		}
	}
	else if (button == GLFW_MOUSE_BUTTON_2)
	{
		if (action == GLFW_PRESS)
		{
			s_clickPointWS = g_camera.ConvertScreenToWorld(ps);
			s_rightMouseDown = true;
		}

		if (action == GLFW_RELEASE)
		{
			s_rightMouseDown = false;
		}
	}
}

static void MouseMotionCallback(GLFWwindow* window, double xd, double yd)
{
	s2Vec2 ps = {float(xd) / s_windowScale, float(yd) / s_windowScale};

	ImGui_ImplGlfw_CursorPosCallback(window, ps.x, ps.y);

	s2Vec2 pw = g_camera.ConvertScreenToWorld(ps);
	s_sample->MouseMove(pw);

	if (s_rightMouseDown)
	{
		s2Vec2 diff = s2Sub(pw, s_clickPointWS);
		g_camera.m_center.x -= diff.x;
		g_camera.m_center.y -= diff.y;
		s_clickPointWS = g_camera.ConvertScreenToWorld(ps);
	}
}

static void ScrollCallback(GLFWwindow* window, double dx, double dy)
{
	ImGui_ImplGlfw_ScrollCallback(window, dx, dy);
	if (ImGui::GetIO().WantCaptureMouse)
	{
		return;
	}

	if (dy > 0)
	{
		g_camera.m_zoom /= 1.1f;
	}
	else
	{
		g_camera.m_zoom *= 1.1f;
	}
}

static void UpdateUI()
{
	float menuWidth = 180.0f;
	if (g_draw.m_showUI)
	{
		ImGui::SetNextWindowPos({g_camera.m_width - menuWidth - 10.0f, 10.0f});
		ImGui::SetNextWindowSize({menuWidth, g_camera.m_height - 20.0f});

		ImGui::Begin("Tools", &g_draw.m_showUI, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);

		if (ImGui::BeginTabBar("ControlTabs", ImGuiTabBarFlags_None))
		{
			if (ImGui::BeginTabItem("Controls"))
			{
				ImGui::SliderInt("Vel Iters", &s_settings.m_velocityIterations, 0, 50);
				ImGui::SliderInt("Pos Iters", &s_settings.m_positionIterations, 0, 50);
				ImGui::SliderFloat("Hertz", &s_settings.m_hertz, 5.0f, 120.0f, "%.0f hz");

				ImGui::Separator();

				ImGui::Checkbox("Sleep", &s_settings.m_enableSleep);
				ImGui::Checkbox("Warm Starting", &s_settings.m_enableWarmStarting);

				ImGui::Separator();

				ImGui::Checkbox("Shapes", &s_settings.m_drawShapes);
				ImGui::Checkbox("Joints", &s_settings.m_drawJoints);
				ImGui::Checkbox("AABBs", &s_settings.m_drawAABBs);
				ImGui::Checkbox("Contact Points", &s_settings.m_drawContactPoints);
				ImGui::Checkbox("Contact Normals", &s_settings.m_drawContactNormals);
				ImGui::Checkbox("Contact Impulses", &s_settings.m_drawContactImpulse);
				ImGui::Checkbox("Friction Impulses", &s_settings.m_drawFrictionImpulse);
				ImGui::Checkbox("Center of Masses", &s_settings.m_drawCOMs);
				ImGui::Checkbox("Statistics", &s_settings.m_drawStats);

				ImVec2 button_sz = ImVec2(-1, 0);
				if (ImGui::Button("Pause (P)", button_sz))
				{
					s_settings.m_pause = !s_settings.m_pause;
				}

				if (ImGui::Button("Single Step (O)", button_sz))
				{
					s_settings.m_singleStep = !s_settings.m_singleStep;
				}

				if (ImGui::Button("Restart (R)", button_sz))
				{
					RestartTest();
				}

				if (ImGui::Button("Quit", button_sz))
				{
					glfwSetWindowShouldClose(g_mainWindow, GL_TRUE);
				}

				ImGui::EndTabItem();
			}

			ImGuiTreeNodeFlags leafNodeFlags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
			leafNodeFlags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;

			ImGuiTreeNodeFlags nodeFlags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;

			if (ImGui::BeginTabItem("Tests"))
			{
				int categoryIndex = 0;
				const char* category = g_sampleEntries[categoryIndex].category;
				int i = 0;
				while (i < g_sampleCount)
				{
					bool categorySelected = strcmp(category, g_sampleEntries[s_settings.m_sampleIndex].category) == 0;
					ImGuiTreeNodeFlags nodeSelectionFlags = categorySelected ? ImGuiTreeNodeFlags_Selected : 0;
					bool nodeOpen = ImGui::TreeNodeEx(category, nodeFlags | nodeSelectionFlags);

					if (nodeOpen)
					{
						while (i < g_sampleCount && strcmp(category, g_sampleEntries[i].category) == 0)
						{
							ImGuiTreeNodeFlags selectionFlags = 0;
							if (s_settings.m_sampleIndex == i)
							{
								selectionFlags = ImGuiTreeNodeFlags_Selected;
							}
							ImGui::TreeNodeEx((void*)(intptr_t)i, leafNodeFlags | selectionFlags, "%s", g_sampleEntries[i].name);
							if (ImGui::IsItemClicked())
							{
								s_selection = i;
							}
							++i;
						}
						ImGui::TreePop();
					}
					else
					{
						while (i < g_sampleCount && strcmp(category, g_sampleEntries[i].category) == 0)
						{
							++i;
						}
					}

					if (i < g_sampleCount)
					{
						category = g_sampleEntries[i].category;
						categoryIndex = i;
					}
				}
				ImGui::EndTabItem();
			}
			ImGui::EndTabBar();
		}

		ImGui::End();

		s_sample->UpdateUI();
	}
}

//
int main(int, char**)
{
#if defined(_WIN32)
	// Enable memory-leak reports
	_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG | _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
	//_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
	//_CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
	//_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
	//_CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
	{
		// Get the current bits
		// int tmp = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);

		// Clear the upper 16 bits and OR in the desired frequency
		// tmp = (tmp & 0x0000FFFF) | _CRTDBG_CHECK_EVERY_16_DF;

		// Set the new bits
		//_CrtSetDbgFlag(tmp);

		//_CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_DELAY_FREE_MEM_DF | _CRTDBG_CHECK_CRT_DF | _CRTDBG_LEAK_CHECK_DF);
		//_CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_DELAY_FREE_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG));
		//_CrtSetDbgFlag(_CRTDBG_DELAY_FREE_MEM_DF | _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG));
	}
	_CrtSetAllocHook(MyAllocHook);
	//_CrtSetBreakAlloc(196);
#endif

	char buffer[128];

	s_settings.Load();
	SortTests();

	glfwSetErrorCallback(glfwErrorCallback);

	g_camera.m_width = s_settings.m_windowWidth;
	g_camera.m_height = s_settings.m_windowHeight;

	if (glfwInit() == 0)
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		return -1;
	}

#if __APPLE__
	const char* glslVersion = "#version 150";
#else
	const char* glslVersion = nullptr;
#endif

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// MSAA
	glfwWindowHint(GLFW_SAMPLES, 4);

	sprintf(buffer, "Solver2D");

	if (GLFWmonitor* primaryMonitor = glfwGetPrimaryMonitor())
	{
#ifdef __APPLE__
		glfwGetMonitorContentScale(primaryMonitor, &s_framebufferScale, &s_framebufferScale);
#else
		glfwGetMonitorContentScale(primaryMonitor, &s_windowScale, &s_windowScale);
#endif
	}

	bool fullscreen = false;
	if (fullscreen)
	{
		g_mainWindow = glfwCreateWindow(int(1920 * s_windowScale), int(1080 * s_windowScale), buffer, glfwGetPrimaryMonitor(), nullptr);
	}
	else
	{
		g_mainWindow =
			glfwCreateWindow(int(g_camera.m_width * s_windowScale), int(g_camera.m_height * s_windowScale), buffer, nullptr, nullptr);
	}

	if (g_mainWindow == nullptr)
	{
		fprintf(stderr, "Failed to open GLFW g_mainWindow.\n");
		glfwTerminate();
		return -1;
	}

#ifdef __APPLE__
	glfwGetWindowContentScale(g_mainWindow, &s_framebufferScale, &s_framebufferScale);
#else
	glfwGetWindowContentScale(g_mainWindow, &s_windowScale, &s_windowScale);
#endif

	glfwMakeContextCurrent(g_mainWindow);

	// Load OpenGL functions using glad
	if (!gladLoadGL())
	{
		fprintf(stderr, "Failed to initialize glad\n");
		glfwTerminate();
		return -1;
	}

	printf("GL %d.%d\n", GLVersion.major, GLVersion.minor);
	printf("OpenGL %s, GLSL %s\n", glGetString(GL_VERSION), glGetString(GL_SHADING_LANGUAGE_VERSION));

	glfwSetWindowSizeCallback(g_mainWindow, ResizeWindowCallback);
	glfwSetKeyCallback(g_mainWindow, KeyCallback);
	glfwSetCharCallback(g_mainWindow, CharCallback);
	glfwSetMouseButtonCallback(g_mainWindow, MouseButtonCallback);
	glfwSetCursorPosCallback(g_mainWindow, MouseMotionCallback);
	glfwSetScrollCallback(g_mainWindow, ScrollCallback);

	g_draw.Create();
	CreateUI(g_mainWindow, glslVersion);

	s_settings.m_sampleIndex = S2_CLAMP(s_settings.m_sampleIndex, 0, g_sampleCount - 1);
	s_selection = s_settings.m_sampleIndex;
	s_sample = g_sampleEntries[s_settings.m_sampleIndex].createFcn(s_settings);

	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);

	float frameTime = 0.0;

	int32_t frame = 0;

	while (!glfwWindowShouldClose(g_mainWindow))
	{
		double time1 = glfwGetTime();

		if (glfwGetKey(g_mainWindow, GLFW_KEY_Z) == GLFW_PRESS)
		{
			// Zoom out
			g_camera.m_zoom = S2_MIN(1.005f * g_camera.m_zoom, 20.0f);
		}
		else if (glfwGetKey(g_mainWindow, GLFW_KEY_X) == GLFW_PRESS)
		{
			// Zoom in
			g_camera.m_zoom = S2_MAX(0.995f * g_camera.m_zoom, 0.02f);
		}

		glfwGetWindowSize(g_mainWindow, &g_camera.m_width, &g_camera.m_height);
		g_camera.m_width = int(g_camera.m_width / s_windowScale);
		g_camera.m_height = int(g_camera.m_height / s_windowScale);

		int bufferWidth, bufferHeight;
		glfwGetFramebufferSize(g_mainWindow, &bufferWidth, &bufferHeight);
		glViewport(0, 0, bufferWidth, bufferHeight);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		double cursorPosX = 0, cursorPosY = 0;
		glfwGetCursorPos(g_mainWindow, &cursorPosX, &cursorPosY);
		ImGui_ImplGlfw_CursorPosCallback(g_mainWindow, cursorPosX / s_windowScale, cursorPosY / s_windowScale);
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui_ImplGlfw_CursorPosCallback(g_mainWindow, cursorPosX / s_windowScale, cursorPosY / s_windowScale);

		ImGuiIO& io = ImGui::GetIO();
		io.DisplaySize.x = float(g_camera.m_width);
		io.DisplaySize.y = float(g_camera.m_height);
		io.DisplayFramebufferScale.x = bufferWidth / float(g_camera.m_width);
		io.DisplayFramebufferScale.y = bufferHeight / float(g_camera.m_height);

		ImGui::NewFrame();

		ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
		ImGui::SetNextWindowSize(ImVec2(float(g_camera.m_width), float(g_camera.m_height)));
		ImGui::SetNextWindowBgAlpha(0.0f);
		ImGui::Begin("Overlay", nullptr,
					 ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_AlwaysAutoResize |
						 ImGuiWindowFlags_NoScrollbar);
		ImGui::End();

		if (g_draw.m_showUI)
		{
			const SampleEntry& entry = g_sampleEntries[s_settings.m_sampleIndex];
			sprintf(buffer, "%s : %s", entry.category, entry.name);
			s_sample->DrawTitle(buffer);
		}

		s_sample->Step(s_settings);

		g_draw.Flush();

		UpdateUI();

		// ImGui::ShowDemoWindow();

		// if (g_draw.m_showUI)
		{
			sprintf(buffer, "%.1f ms", 1000.0f * frameTime);

			ImGui::Begin("Overlay", nullptr,
						ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_AlwaysAutoResize |
							ImGuiWindowFlags_NoScrollbar);
			ImGui::SetCursorPos(ImVec2(5.0f, g_camera.m_height - 20.0f));
			ImGui::TextColored(ImColor(153, 230, 153, 255), buffer);
			ImGui::End();
		}

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(g_mainWindow);

		if (s_selection != s_settings.m_sampleIndex)
		{
			s_settings.m_sampleIndex = s_selection;
			delete s_sample;
			s_sample = g_sampleEntries[s_settings.m_sampleIndex].createFcn(s_settings);
			g_camera.ResetView();
		}

		glfwPollEvents();

		// Limit frame rate to 60Hz
		double time2 = glfwGetTime();
		double targetTime = time1 + 1.0f / 60.0f;
		int loopCount = 0;
		while (time2 < targetTime)
		{
#if defined(_WIN32)
			Sleep(0.0f);
#endif
			time2 = glfwGetTime();
			++loopCount;
		}

		frameTime = (float)(time2 - time1);
		// if (frame % 17 == 0)
		//{
		//	printf("loop count = %d, frame time = %.1f\n", loopCount, 1000.0f * frameTime);
		// }
		++frame;
	}

	delete s_sample;
	s_sample = nullptr;

	g_draw.Destroy();

	DestroyUI();
	glfwTerminate();

	s_settings.Save();

#if defined(_WIN32)
	_CrtDumpMemoryLeaks();
#endif

	return 0;
}