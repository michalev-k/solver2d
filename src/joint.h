// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "solver2d/id.h"
#include "solver2d/types.h"

#include "pool.h"

#include <stdint.h>

typedef struct s2DebugDraw s2DebugDraw;
typedef struct s2StepContext s2StepContext;
typedef struct s2World s2World;

typedef enum s2JointType
{
	s2_revoluteJoint,
	s2_mouseJoint,
} s2JointType;

/// A joint edge is used to connect bodies and joints together
/// in a joint graph where each body is a node and each joint
/// is an edge. A joint edge belongs to a doubly linked list
/// maintained in each attached body. Each joint has two joint
/// nodes, one for each attached body.
typedef struct s2JointEdge
{
	int32_t bodyIndex;
	int32_t prevKey;
	int32_t nextKey;
} s2JointEdge;

typedef struct s2MouseJoint
{
	s2Vec2 targetA;
	float stiffness;
	float damping;
	float beta;

	// Solver shared
	s2Vec2 impulse;
	float maxForce;
	float gamma;

	// Solver temp
	s2Vec2 rB;
	s2Vec2 localCenterB;
	float invMassB;
	float invIB;
	s2Mat22 mass;
	s2Vec2 C;
} s2MouseJoint;

typedef struct s2RevoluteJoint
{
	// Solver shared
	s2Vec2 impulse;
	float motorImpulse;
	float lowerImpulse;
	float upperImpulse;
	bool enableMotor;
	float maxMotorTorque;
	float motorSpeed;
	bool enableLimit;
	float referenceAngle;
	float lowerAngle;
	float upperAngle;

	// Solver temp
	s2Vec2 rA;
	s2Vec2 rB;
	s2Vec2 localCenterA;
	s2Vec2 localCenterB;
	float invMassA;
	float invMassB;
	float invIA;
	float invIB;
	s2Mat22 K;
	float angle;
	float axialMass;
} s2RevoluteJoint;

/// The base joint class. Joints are used to constraint two bodies together in
/// various fashions. Some joints also feature limits and motors.
typedef struct s2Joint
{
	s2Object object;

	s2JointType type;

	s2JointEdge edges[2];

	s2Vec2 localAnchorA;
	s2Vec2 localAnchorB;

	union
	{
		s2MouseJoint mouseJoint;
		s2RevoluteJoint revoluteJoint;
	};

	bool collideConnected;
} s2Joint;

void s2InitVelocityConstraints(s2Joint* joint, s2StepContext* data);
void s2SolveVelocityConstraints(s2Joint* joint, s2StepContext* data);

// This returns true if the position errors are within tolerance.
bool s2SolvePositionConstraints(s2Joint* joint, s2StepContext* data);

void s2DrawJoint(s2DebugDraw* draw, s2World* world, s2Joint* joint);