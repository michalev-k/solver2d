// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/types.h"

struct Settings
{
	void Save();
	void Load();

	int sampleIndex = 0;
	int windowWidth = 1920;
	int windowHeight = 1080;
	float hertz = 60.0f;
	float timeStep = 1.0f / 60.0f;
	int primaryIterations = 4;
	int secondaryIterations = 2;
	int relaxationIterations = 2;
	int groupCount = 4;
	int multiSteps = 1;
	int textLine = 0;
	int textIncrement = 18;
	bool enabledSolvers[s2_solverTypeCount] = {};
	bool enableWarmStarting = true;
	bool solveJointsWithContacts = false;
	bool drawShapes = true;
	bool drawJoints = true;
	bool drawAABBs = false;
	bool drawContactPoints = false;
	bool drawContactNormals = false;
	bool drawContactImpulse = false;
	bool drawFrictionImpulse = false;
	bool drawMass = false;
	bool drawStats = false;
	bool pause = false;
	bool singleStep = false;
	bool restart = false;
};
