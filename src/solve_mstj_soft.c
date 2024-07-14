#include "allocate.h"
#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include <stdbool.h>

typedef struct SplitBody {
	s2Vec2 dv;
	float dw;
	int32_t prev;
} SplitBody;

typedef struct SplitConstraint {
	SplitBody el[2];
} SplitConstraint;

typedef struct SolveTestContext {
	SplitConstraint* splitConstraints;

	s2ContactConstraint* constraints;
	int constraintCount;

	s2Joint* joints;
	int jointCount;

	s2Body* bodies;
	int bodyCount;

	s2StepContext *context;

} SolveMstjSoftContext;

void SplitBodies(s2World* world, SolveMstjSoftContext *sctx) {
	s2ContactConstraint* constraints = sctx->constraints;
	int constraintCount = sctx->constraintCount;

	s2Joint* joints = sctx->joints;
	int jointCount = sctx->jointCount;

	s2Body* bodies = sctx->bodies;	
	int bodyCapacity = sctx->bodyCount;
	
	// reset splitLists for bodies
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectValid(&body->object) == false)
		{
			continue;
		}
		body->splitList = S2_NULL_INDEX;
		body->splitCount = 0;

		body->splitListJoints = S2_NULL_INDEX;
		body->splitJointsCount = 0;
	}

	SplitConstraint* solveDatas =
		s2AllocateStackItem(world->stackAllocator, (constraintCount + jointCount) * sizeof(SplitConstraint), "solve data");


	for (int i = 0; i < constraintCount; ++i) {
		s2ContactConstraint* constraint = constraints + i;	
		constraint->indexA = constraint->contact->edges[0].bodyIndex;
		constraint->indexB = constraint->contact->edges[1].bodyIndex;

		s2Body* bodyA = bodies + constraint->indexA;
		s2Body* bodyB = bodies + constraint->indexB;

		SplitConstraint *solveData = solveDatas + i;

		int keyA = i << 1;
		int keyB = keyA + 1;

		solveData->el[0].prev = bodyA->splitList;
		bodyA->splitList = keyA;

		solveData->el[1].prev = bodyB->splitList;
		bodyB->splitList = keyB;
		
		bodyA->splitCount++;
		bodyB->splitCount++;
	}

	for (int i = 0; i < jointCount; i++) {
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object)) { continue; }

		s2Body* bodyA = bodies + joint->edges[0].bodyIndex;
		s2Body* bodyB = bodies + joint->edges[1].bodyIndex;

		int32_t solveBodyIndex = constraintCount + i;
		joint->solveBodyIndex = solveBodyIndex;

		SplitConstraint *solveData = solveDatas + solveBodyIndex;

		int keyA = (i + constraintCount) << 1;
		int keyB = keyA + 1;

		if (sctx->context->solveJointsWithContacts)
		{
			solveData->el[0].prev = bodyA->splitList;
			bodyA->splitList = keyA;

			solveData->el[1].prev = bodyB->splitList;
			bodyB->splitList = keyB;
		}
		else
		{
			solveData->el[0].prev = bodyA->splitListJoints;
			bodyA->splitListJoints = keyA;

			solveData->el[1].prev = bodyB->splitListJoints;
			bodyB->splitListJoints = keyB;
		}

		bodyA->splitJointsCount++;
		bodyB->splitJointsCount++;
	}

	sctx->splitConstraints = solveDatas;
	return;
}

void s2MergeSplitContactVelocities(s2World* world, SolveMstjSoftContext *sctx)
{
	int constraintCount = sctx->constraintCount;

	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectValid(&body->object) == false)
		{
			continue;
		}

		if (body->type == s2_staticBody)
		{
			continue;
		}

		if (body->splitList == S2_NULL_INDEX) {
			continue;
		}

		SplitConstraint* solveDatas = sctx->splitConstraints;

		s2Vec2 dv = S2_ZERO_INIT;
		float dw = 0;

		int32_t splitKey = body->splitList;
		while (splitKey != S2_NULL_INDEX) {
			int id = splitKey >> 1;
			int elid = splitKey & 1;

			SplitConstraint *solveData = solveDatas + id;
			SplitBody *sb = &solveData->el[elid];

			dv = s2Add(dv, sb->dv);
			dw += sb->dw;

			sb->dv = s2Vec2_zero;
			sb->dw = 0;

			splitKey = sb->prev;
		}

		body->linearVelocity = s2Add(body->linearVelocity, dv);
		body->angularVelocity += dw;
	}
}

void s2MergeSplitJointVelocities(s2World* world, SolveMstjSoftContext *sctx)
{
	int constraintCount = sctx->constraintCount;

	s2Body* bodies = sctx->bodies;
	int bodyCount = sctx->bodyCount;
	for (int i = 0; i < bodyCount; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectValid(&body->object) == false)
		{
			continue;
		}

		if (body->type == s2_staticBody)
		{
			continue;
		}

		if (body->splitListJoints == S2_NULL_INDEX) {
			continue;
		}

		SplitConstraint* solveDatas = sctx->splitConstraints;

		s2Vec2 dv = S2_ZERO_INIT;
		float dw = 0;

		int32_t splitKey = body->splitListJoints;
		while (splitKey != S2_NULL_INDEX) {
			int id = splitKey >> 1;
			int elid = splitKey & 1;

			SplitConstraint *solveData = solveDatas + id;
			SplitBody *sb = &solveData->el[elid];

			dv = s2Add(dv, sb->dv);
			dw += sb->dw;

			sb->dv = s2Vec2_zero;
			sb->dw = 0;

			splitKey = sb->prev;
		}

		body->linearVelocity = s2Add(body->linearVelocity, dv);
		body->angularVelocity += dw;
	}
}

static void s2PrepareContacts_MSTJ_Soft(s2World* world, SolveMstjSoftContext *sctx, float h, float hertzDynamic, float hertzStatic)
{
	s2ContactConstraint* constraints = sctx->constraints;
	int constraintCount = sctx->constraintCount;
	bool warmStart = sctx->context->warmStart;

	s2Body* bodies = sctx->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;

		s2Contact* contact = constraint->contact;
		const s2Manifold* manifold = &contact->manifold;
		int pointCount = manifold->pointCount;
		S2_ASSERT(0 < pointCount && pointCount <= 2);
		int indexA = contact->edges[0].bodyIndex;
		int indexB = contact->edges[1].bodyIndex;

		constraint->indexA = indexA;
		constraint->indexB = indexB;
		constraint->normal = manifold->normal;
		constraint->friction = contact->friction;
		constraint->pointCount = pointCount;

		s2Body* bodyA = bodies + indexA;
		s2Body* bodyB = bodies + indexB;

		float mA = bodyA->invMass * bodyA->splitCount;
		float iA = bodyA->invI    * bodyA->splitCount;
		float mB = bodyB->invMass * bodyB->splitCount;
		float iB = bodyB->invI    * bodyB->splitCount;

		s2Rot qA = bodyA->rot;
		s2Rot qB = bodyB->rot;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(constraint->normal);

		for (int j = 0; j < pointCount; ++j)
		{
			const s2ManifoldPoint* mp = manifold->points + j;
			s2ContactConstraintPoint* cp = constraint->points + j;

			if (warmStart)
			{
				cp->normalImpulse = mp->normalImpulse;
				cp->tangentImpulse = mp->tangentImpulse;
			}
			else
			{
				cp->normalImpulse = 0.0f;
				cp->tangentImpulse = 0.0f;
			}

			cp->localAnchorA = s2Sub(mp->localAnchorA, bodyA->localCenter);
			cp->localAnchorB = s2Sub(mp->localAnchorB, bodyB->localCenter);
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			cp->separation = mp->separation;
			cp->adjustedSeparation = mp->separation - s2Dot(s2Sub(rB, rA), normal);

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);
			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);
			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
			cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

			float unsplitKNormal = (bodyA->invMass + bodyA->invI * rnA * rnA) + (bodyB->invMass + bodyB->invI * rnB * rnB);
			cp->splitKFactor = unsplitKNormal / kNormal;
		
			// soft contact
			// should use the substep not the full time step			
			// Stiffer for dynamic vs static
			float contactHertz = (mA == 0.0f || mB == 0.0f) ? hertzStatic : hertzDynamic;
			
			const float zeta = 10.0f;
			float omega = 2.0f * s2_pi * contactHertz;
			float c = h * omega * (2.0f * zeta + h * omega);
			cp->biasCoefficient = omega / (2.0f * zeta + h * omega);
			cp->impulseCoefficient = 1.0f / (1.0f + c);
			cp->massCoefficient = c * cp->impulseCoefficient;
		}
	}
}

void s2WarmStartContacts_MSTJ_Soft(s2World* world, SolveMstjSoftContext *sctx)
{
	s2ContactConstraint* constraints = sctx->constraints;
	int constraintCount = sctx->constraintCount;

	s2Body* bodies = sctx->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;
			
		int pointCount = constraint->pointCount;
		S2_ASSERT(0 < pointCount && pointCount <= 2);

		s2Body* bodyA = bodies + constraint->indexA;
		s2Body* bodyB = bodies + constraint->indexB;

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;

		s2Rot qA = bodyA->rot;
		s2Rot qB = bodyB->rot;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);
		
		SplitConstraint *solveData = sctx->splitConstraints + i;
		SplitBody *splitBodyA = &solveData->el[0];
		SplitBody *splitBodyB = &solveData->el[1];

		s2Vec2 dvA = splitBodyA->dv;
		float  dwA = splitBodyA->dw;
		s2Vec2 dvB = splitBodyB->dv;
		float  dwB = splitBodyB->dw;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Current anchors (to support TGS)
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			s2Vec2 P = s2Add(s2MulSV(cp->normalImpulse, normal), s2MulSV(cp->tangentImpulse, tangent));

			dvA = s2MulSub(dvA, mA, P);
			dwA -= iA * s2Cross(rA, P);
			
			dvB = s2MulAdd(dvB, mB, P);
			dwB += iB * s2Cross(rB, P);
		}

		splitBodyA->dv = dvA;
		splitBodyA->dw = dwA;
		splitBodyB->dv = dvB;
		splitBodyB->dw = dwB;
	}
}
static void s2SolveContacts_MSTJ_Soft(s2World* world, SolveMstjSoftContext *sctx, bool useBias, float baum, float inv_h)
{
	s2ContactConstraint* constraints = sctx->constraints;
	int constraintCount = sctx->constraintCount;

	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;
		
		s2Body* bodyA = bodies + constraint->indexA;
		s2Body* bodyB = bodies + constraint->indexB;

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;

		int pointCount = constraint->pointCount;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);
		float friction = constraint->friction;

		SplitConstraint *solveData = sctx->splitConstraints + i;
		SplitBody *splitBodyA = &solveData->el[0];
		SplitBody *splitBodyB = &solveData->el[1];

		s2Vec2 dvA = splitBodyA->dv;
		float  dwA = splitBodyA->dw;
		s2Vec2 dvB = splitBodyB->dv;
		float  dwB = splitBodyB->dw;

		s2Vec2 dcA = bodyA->deltaPosition;
		s2Rot qA = bodyA->rot;
		s2Vec2 dcB = bodyB->deltaPosition;
		s2Rot qB = bodyB->rot;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Vec2 vA = bodyA->linearVelocity;
			float wA = bodyA->angularVelocity;
			s2Vec2 vB = bodyB->linearVelocity;
			float wB = bodyB->angularVelocity;

			vA =  s2Add(vA, dvA);
			wA += dwA;
			vB =  s2Add(vB, dvB);
			wB += dwB;

			// anchor points
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// compute current separation
			s2Vec2 ds = s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA));
			float s = s2Dot(ds, normal) + cp->adjustedSeparation;

			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (s > 0.0f)
			{
				// Speculative
				bias = s * inv_h;
			}
			else if (useBias)
			{
				bias = S2_MAX(cp->biasCoefficient * s, -s2_maxBaumgarteVelocity);
				massScale = cp->massCoefficient;
				impulseScale = cp->impulseCoefficient;
			}


			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			float vn = s2Dot(s2Sub(vrB, vrA), normal);

			// Compute normal impulse
			float impulse = -cp->normalMass * massScale * (vn + bias) - impulseScale * cp->normalImpulse;

			// Clamp the accumulated impulse
			float newImpulse = S2_MAX(cp->normalImpulse + impulse, 0.0f);
			impulse = newImpulse - cp->normalImpulse;
			cp->normalImpulse = newImpulse;

			// // Apply contact impulse
			s2Vec2 P = s2MulSV(impulse, normal);
				
			dvA = s2MulSub(dvA, mA, P);
			dwA -= iA * s2Cross(rA, P);
			
			dvB = s2MulAdd(dvB, mB, P);
			dwB += iB * s2Cross(rB, P);
		}

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;
			
			s2Vec2 vA = bodyA->linearVelocity;
			float wA = bodyA->angularVelocity;
			s2Vec2 vB = bodyB->linearVelocity;
			float wB = bodyB->angularVelocity;
			
			vA =  s2Add(vA, dvA);
			wA += dwA;
			vB =  s2Add(vB, dvB);
			wB += dwB;
			
			// Current anchor points
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);
			
			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			s2Vec2 dv = s2Sub(vrB, vrA);
			
			// Compute tangent force
			float vt = s2Dot(dv, tangent);
			float lambda = cp->tangentMass * (-vt);
			
			// Clamp the accumulated force
			float maxFriction = friction * cp->normalImpulse;
			float newImpulse = S2_CLAMP(cp->tangentImpulse + lambda, -maxFriction, maxFriction);
			lambda = newImpulse - cp->tangentImpulse;
			cp->tangentImpulse = newImpulse;
			
			// Apply contact impulse
			s2Vec2 P = s2MulSV(lambda, tangent);
			
			dvA = s2MulSub(dvA, mA, P);
			dwA -= iA * s2Cross(rA, P);
			
			dvB = s2MulAdd(dvB, mB, P);
			dwB += iB * s2Cross(rB, P);
		}

		splitBodyA->dv = dvA;
		splitBodyA->dw = dwA;
		splitBodyB->dv = dvB;
		splitBodyB->dw = dwB;
	}
	
}


void s2WarmStartRevoluteTEST(s2Joint* base, SolveMstjSoftContext *sctx)
{
	s2StepContext *context = sctx->context;
	S2_ASSERT(base->type == s2_revoluteJoint);

	int32_t indexA = base->edges[0].bodyIndex;
	int32_t indexB = base->edges[1].bodyIndex;
	S2_ASSERT(0 <= indexA && indexA < context->bodyCapacity);
	S2_ASSERT(0 <= indexB && indexB < context->bodyCapacity);

	s2Body* bodyA = context->bodies + indexA;
	s2Body* bodyB = context->bodies + indexB;
	S2_ASSERT(bodyA->object.index == bodyA->object.next);
	S2_ASSERT(bodyB->object.index == bodyB->object.next);

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Rot qA = bodyA->rot;
	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;

	s2Rot qB = bodyB->rot;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
	s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	float axialImpulse = joint->motorImpulse + joint->lowerImpulse - joint->upperImpulse;
	s2Vec2 P = {joint->impulse.x, joint->impulse.y};

	vA = s2MulSub(vA, mA, P);
	wA -= iA * (s2Cross(rA, P) + axialImpulse);

	vB = s2MulAdd(vB, mB, P);
	wB += iB * (s2Cross(rB, P) + axialImpulse);

	SplitConstraint *solveData = sctx->splitConstraints + base->solveBodyIndex;
	SplitBody *splitBodyA = &solveData->el[0];
	SplitBody *splitBodyB = &solveData->el[1];

	splitBodyA->dv = s2Sub(vA, bodyA->linearVelocity);
	splitBodyA->dw = wA - bodyA->angularVelocity;
	splitBodyB->dv = s2Sub(vB, bodyB->linearVelocity);
	splitBodyB->dw = wB - bodyB->angularVelocity;
}


void s2PrepareRevolute_TEST(s2Joint* base, s2StepContext* context, float h, float hertz, bool warmStart)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	int32_t indexA = base->edges[0].bodyIndex;
	int32_t indexB = base->edges[1].bodyIndex;
	S2_ASSERT(0 <= indexA && indexA < context->bodyCapacity);
	S2_ASSERT(0 <= indexB && indexB < context->bodyCapacity);

	s2Body* bodyA = context->bodies + indexA;
	s2Body* bodyB = context->bodies + indexB;
	S2_ASSERT(bodyA->object.index == bodyA->object.next);
	S2_ASSERT(bodyB->object.index == bodyB->object.next);

	s2RevoluteJoint* joint = &base->revoluteJoint;
	joint->localAnchorA = s2Sub(base->localOriginAnchorA, bodyA->localCenter);
	joint->invMassA = bodyA->invMass;
	joint->invIA = bodyA->invI;

	joint->localAnchorB = s2Sub(base->localOriginAnchorB, bodyB->localCenter);
	joint->invMassB = bodyB->invMass;
	joint->invIB = bodyB->invI;

	joint->centerDiff0 = s2Sub(bodyB->position, bodyA->position);

	// J = [-I -r1_skew I r2_skew]
	// r_skew = [-ry; rx]

	// Matlab
	// K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x]
	//     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB]

	float mA = bodyA->invMass * bodyA->splitJointsCount;
	float iA = bodyA->invI    * bodyA->splitJointsCount;
	float mB = bodyB->invMass * bodyB->splitJointsCount;
	float iB = bodyB->invI    * bodyB->splitJointsCount;

	s2Rot qA = bodyA->rot;
	s2Rot qB = bodyB->rot;

	s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
	s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

	s2Mat22 K;
	K.cx.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
	K.cy.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
	K.cx.y = K.cy.x;
	K.cy.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
	
	joint->pivotMass = s2GetInverse22(K);

	{
		const float zeta = 1.0f;
		float omega = 2.0f * s2_pi * hertz;
		joint->biasCoefficient = omega / (2.0f * zeta + h * omega);
		float c = h * omega * (2.0f * zeta + h * omega);
		joint->impulseCoefficient = 1.0f / (1.0f + c);
		joint->massCoefficient = c * joint->impulseCoefficient;
	}

	joint->axialMass = iA + iB;
	bool fixedRotation;
	if (joint->axialMass > 0.0f)
	{
		joint->axialMass = 1.0f / joint->axialMass;
		fixedRotation = false;
	}
	else
	{
		fixedRotation = true;
	}

	if (joint->enableLimit == false || fixedRotation || warmStart == false)
	{
		joint->lowerImpulse = 0.0f;
		joint->upperImpulse = 0.0f;
	}

	if (joint->enableMotor == false || fixedRotation || warmStart == false)
	{
		joint->motorImpulse = 0.0f;
	}

	if (warmStart == false)
	{
		joint->impulse = s2Vec2_zero;
	}
}

void s2SolveRevolute_TEST(s2Joint* base, SolveMstjSoftContext *sctx, float h, float inv_h, bool useBias)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Body* bodyA = sctx->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = sctx->bodies + base->edges[1].bodyIndex;

	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	bool fixedRotation = (iA + iB == 0.0f);

	// Solve motor constraint.
	if (joint->enableMotor && fixedRotation == false)
	{
		float Cdot = wB - wA - joint->motorSpeed;
		float impulse = -joint->axialMass * Cdot;
		float oldImpulse = joint->motorImpulse;
		float maxImpulse = h * joint->maxMotorTorque;
		joint->motorImpulse = S2_CLAMP(joint->motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = joint->motorImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	if (joint->enableLimit && fixedRotation == false)
	{
		float jointAngle = s2RelativeAngle(bodyB->rot, bodyA->rot) - joint->referenceAngle;

		// Lower limit
		{
			float C = jointAngle - joint->lowerAngle;
			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (C > 0.0f)
			{
				// speculation
				bias = C * inv_h;
			}
			else if (useBias)
			{
				bias = joint->biasCoefficient * C;
				massScale = joint->massCoefficient;
				impulseScale = joint->impulseCoefficient;
			}

			float Cdot = wB - wA;
			float impulse = -joint->axialMass * massScale * (Cdot + bias) - impulseScale * joint->lowerImpulse;
			float oldImpulse = joint->lowerImpulse;
			joint->lowerImpulse = S2_MAX(joint->lowerImpulse + impulse, 0.0f);
			impulse = joint->lowerImpulse - oldImpulse;

			wA -= iA * impulse;
			wB += iB * impulse;
		}

		// Upper limit
		// Note: signs are flipped to keep C positive when the constraint is satisfied.
		// This also keeps the impulse positive when the limit is active.
		{
			float C = joint->upperAngle - jointAngle;

			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (C > 0.0f)
			{
				// speculation
				bias = C * inv_h;
			}
			else if (useBias)
			{
				bias = joint->biasCoefficient * C;
				massScale = joint->massCoefficient;
				impulseScale = joint->impulseCoefficient;
			}

			float Cdot = wA - wB;
			float impulse = -joint->axialMass * massScale * (Cdot + bias) - impulseScale * joint->lowerImpulse;
			float oldImpulse = joint->upperImpulse;
			joint->upperImpulse = S2_MAX(joint->upperImpulse + impulse, 0.0f);
			impulse = joint->upperImpulse - oldImpulse;

			wA += iA * impulse;
			wB -= iB * impulse;
		}
	}

	// Solve point-to-point constraint
	{
		// Update anchors for TGS solvers.
		// Anchors are wastfully recomputed for PGS solvers or relax stages.
		s2Rot qA = bodyA->rot;
		s2Rot qB = bodyB->rot;
		s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

		s2Vec2 Cdot = s2Sub(s2Add(vB, s2CrossSV(wB, rB)), s2Add(vA, s2CrossSV(wA, rA)));

		s2Vec2 bias = s2Vec2_zero;
		float massScale = 1.0f;
		float impulseScale = 0.0f;
		if (useBias)
		{
			s2Vec2 dcA = bodyA->deltaPosition;
			s2Vec2 dcB = bodyB->deltaPosition;

			s2Vec2 separation = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);
			bias = s2MulSV(joint->biasCoefficient, separation);
			massScale = joint->massCoefficient;
			impulseScale = joint->impulseCoefficient;
		}
		
		#if S2_FRESH_PIVOT_MASS == 1
		s2Mat22 K;
		K.cx.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
		K.cy.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
		K.cx.y = K.cy.x;
		K.cy.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
		s2Vec2 b = s2Solve22(K, s2Add(Cdot, bias));
		#else
		s2Vec2 b = s2MulMV(joint->pivotMass, s2Add(Cdot, bias));
		#endif

		s2Vec2 impulse;
		impulse.x = -massScale * b.x - impulseScale * joint->impulse.x;
		impulse.y = -massScale * b.y - impulseScale * joint->impulse.y;
		joint->impulse.x += impulse.x;
		joint->impulse.y += impulse.y;

		vA = s2MulSub(vA, mA, impulse);
		wA -= iA * s2Cross(rA, impulse);
		vB = s2MulAdd(vB, mB, impulse);
		wB += iB * s2Cross(rB, impulse);
	}

	SplitConstraint *solveData = sctx->splitConstraints + base->solveBodyIndex;
	SplitBody *splitBodyA = &solveData->el[0];
	SplitBody *splitBodyB = &solveData->el[1];

	splitBodyA->dv = s2Sub(vA, bodyA->linearVelocity);
	splitBodyA->dw = wA - bodyA->angularVelocity;
	splitBodyB->dv = s2Sub(vB, bodyB->linearVelocity);
	splitBodyB->dw = wB - bodyB->angularVelocity;
}


void s2PrepareJoints_MSTJ_Soft(SolveMstjSoftContext *sctx, float h, float hertz, bool warmStart)
{
	s2StepContext* context = sctx->context;

	s2Joint *joints = sctx->joints;
	int jointCapacity = sctx->jointCount;

	for (int i = 0; i < jointCapacity; ++i) {
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object)) { continue; };
	
		switch (joint->type)
		{
			case s2_mouseJoint:
				s2PrepareMouse(joint, context);
				break;

			case s2_revoluteJoint:
				s2PrepareRevolute_TEST(joint, context, h, hertz, warmStart);
				break;

			default:
				S2_ASSERT(false);
		}	
	}
}

void s2SolveJoints_MSTJ_Soft(SolveMstjSoftContext *sctx, float h, float inv_h, bool useBias)
{
	s2StepContext* context = sctx->context;

	s2Joint *joints = sctx->joints;
	int jointCapacity = sctx->jointCount;

	for (int i = 0; i < jointCapacity; ++i) {
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object)) { continue; };

		switch (joint->type)
		{
			case s2_mouseJoint:
				if (useBias)
				{
					s2SolveMouse(joint, context);
				}
				break;

			case s2_revoluteJoint:
				s2SolveRevolute_TEST(joint, sctx, h, inv_h, useBias);
				break;

			default:
				S2_ASSERT(false);
		}
	}
}

void s2WarmStartJoints_MSTJ_Soft(SolveMstjSoftContext *sctx)
{
	s2StepContext* context = sctx->context;

	s2Joint *joints = sctx->joints;
	int jointCapacity = sctx->jointCount;

	for (int i = 0; i < jointCapacity; ++i) {
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object)) { continue; };

		switch (joint->type)
		{
			case s2_mouseJoint:
				s2WarmStartMouse(joint, context);
				break;

			case s2_revoluteJoint:
				s2WarmStartRevoluteTEST(joint, sctx);
				break;

			default:
				S2_ASSERT(false);
		}
	}
}

// mass-split temporal jacobi with soft-contacts
// in short -> mstj_soft 

void s2Solve_MSTJ_Soft(s2World* world, s2StepContext* context)
{
	s2Contact* contacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;

	s2ContactConstraint* constraints =
		s2AllocateStackItem(world->stackAllocator, contactCapacity * sizeof(s2ContactConstraint), "constraint");

	int constraintCount = 0;
	for (int i = 0; i < contactCapacity; ++i)
	{
		s2Contact* contact = contacts + i;
		if (s2IsFree(&contact->object)) { continue; }

		if (contact->manifold.pointCount == 0) { continue; }

		constraints[constraintCount].contact = contact;
		constraints[constraintCount].contact->manifold.constraintIndex = constraintCount;
		constraintCount += 1;
	}

	float h = context->h;
	float inv_h = context->inv_h;
	
	int substepCount = context->iterations;
	int innerSolveCount = context->extraIterations;
	int relaxCount = context->relaxationIterations;
	
	bool warmStart = context->warmStart;
	bool solveJointsWithContacts = context->solveJointsWithContacts;

	SolveMstjSoftContext sctx;
	sctx.constraints = constraints;
	sctx.constraintCount = constraintCount;
	sctx.joints = world->joints;
	sctx.jointCount = world->jointPool.capacity;
	sctx.bodies = world->bodies;
	sctx.bodyCount = world->bodyPool.capacity;
	sctx.context = context;

	// Note: s2_contactHertz by itself is not enough
	// objects will sink just a bit too much, and this will lead e.g the card-house to collapse
	// with a slightly bigger value it remains stable
	float contactHertzDynamic = S2_MIN(s2_contactHertz * 1.2f, 0.25f * inv_h);
	float contactHertzStatic = S2_MIN(s2_contactHertz * 2.0f, 0.5f * inv_h);

	// if you increase the jointHertz the ragdoll-parts will seperate less
	// the joint-grid-sample, if you pull it hard enough will circle forever though
	// so increasing should be done with caution
	float jointHertz = S2_MIN(s2_jointHertz * 1.0f, 0.125f * inv_h);

	// Note: I have made the body-splitting here seperate for demonstration purposes 
	// Every body already knows about its constraints. In an actual implementation you would design around this.
	SplitBodies(world, &sctx);

	s2PrepareContacts_MSTJ_Soft(world, &sctx, h, contactHertzDynamic, contactHertzStatic);
	s2PrepareJoints_MSTJ_Soft(&sctx, h, jointHertz, context->warmStart);

	for (int substep = 0; substep < substepCount; ++substep)
	{
		s2IntegrateVelocities(world, h);		
		
		if (warmStart) {
			s2WarmStartJoints_MSTJ_Soft(&sctx);

			if (!solveJointsWithContacts) {
				s2MergeSplitJointVelocities(world, &sctx);
			}

			s2WarmStartContacts_MSTJ_Soft(world, &sctx);
			s2MergeSplitContactVelocities(world, &sctx);
		}

		for (int extra = 0; extra < innerSolveCount; ++extra)
		{		
			s2SolveJoints_MSTJ_Soft(&sctx, h, inv_h, true);

			if (!solveJointsWithContacts) {
				s2MergeSplitJointVelocities(world, &sctx);
			}
			
			s2SolveContacts_MSTJ_Soft(world, &sctx, true, s2_baumgarte, inv_h);
			s2MergeSplitContactVelocities(world, &sctx);
		}

		s2IntegratePositions(world, h);	

		for (int rel = 0; rel < relaxCount; rel++)
		{
			s2SolveJoints_MSTJ_Soft(&sctx, h, inv_h, false);

			if (!solveJointsWithContacts) {
				s2MergeSplitJointVelocities(world, &sctx);
			}

			s2SolveContacts_MSTJ_Soft(world, &sctx, false, s2_baumgarte, inv_h);
			s2MergeSplitContactVelocities(world, &sctx);		
		}
	}
	
	s2FinalizePositions(world);
	s2StoreContactImpulses(constraints, constraintCount);
	
	s2FreeStackItem(world->stackAllocator, sctx.splitConstraints);
	s2FreeStackItem(world->stackAllocator, constraints);
}