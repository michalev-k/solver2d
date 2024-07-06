#include "allocate.h"
#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include <stdbool.h>

typedef struct s2SplitBody {
	s2Vec2 deltaV;
	float deltaW;
	int32_t next;
} s2SplitBody;

typedef struct s2ContactGroup {
	s2ContactConstraint* constraints;
	int constraintCount;
} s2ContactGroup;

typedef struct s2SplitContext {
	s2ContactGroup *groups;
	int groupsCount;
	s2SplitBody *splitBodies;
	int splitBodiesCount;
} s2SplitContext;

void s2MergeSplitVelocities(s2World* world, s2SplitContext *sctx)
{
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

		if (body->splitCount == 0) {
			continue;
		}

		s2Vec2 v = body->linearVelocity;
		float w = body->angularVelocity;

		int32_t splitKey = body->splitList;
		while (splitKey != S2_NULL_INDEX)
		{
			s2SplitBody *splitBody = sctx->splitBodies + splitKey;

			v = s2Add(v, splitBody->deltaV);
			w += splitBody->deltaW;

			splitBody->deltaV = s2Vec2_zero;
			splitBody->deltaW = 0;

			splitKey = splitBody->next;
		}

		body->linearVelocity = v;
		body->angularVelocity = w;
	}
}

int s2FindSplitBody(s2Body *body, int groupSplitBodyOffset, s2SplitBody *splitBodies, int *splitBodyCount)
{
	int32_t splitIndex = body->splitList;

	if (splitIndex == S2_NULL_INDEX)
	{
		S2_ASSERT(body->splitCount == 0);

		splitIndex = (*splitBodyCount)++;
		s2SplitBody *splitBody = splitBodies + splitIndex;
		splitBody->next = S2_NULL_INDEX;

		body->splitList = splitIndex;
		body->splitCount = 1;
		return splitIndex;
	}

	while (true)
	{
		if (splitIndex >= groupSplitBodyOffset)
		{
			// this splitIndex must have been created for this group 
			return splitIndex;
		}
	
		s2SplitBody *splitBody = splitBodies + splitIndex;
		if (splitBody->next == S2_NULL_INDEX)
		{
			break;
		}

		splitIndex = splitBody->next;
	}

	s2SplitBody *prev = splitBodies + splitIndex;
	
	splitIndex = (*splitBodyCount)++;
	s2SplitBody *splitBody = splitBodies + splitIndex;
	splitBody->next = S2_NULL_INDEX;

	prev->next = splitIndex;
	body->splitCount++;
	return splitIndex;
}

s2SplitContext s2SplitIntoGroups(s2World* world, s2ContactConstraint* constraints, int constraintCount, int groupsCount)
{
	s2Body* bodies = world->bodies;	
	int bodyCapacity = world->bodyPool.capacity;
	
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
	}

	int groupSize = (constraintCount / groupsCount) + 1;

	s2ContactGroup* groups =
		s2AllocateStackItem(world->stackAllocator, groupsCount * sizeof(s2ContactGroup), "contact groups");

	s2SplitBody* splitBodies =
	s2AllocateStackItem(world->stackAllocator, 2 * constraintCount * sizeof(s2SplitBody), "split bodies");

	int splitBodiesCount = 0;

	for (int g = 0; g < groupsCount; g++) {
		s2ContactGroup *group = &groups[g];
		*group = (s2ContactGroup){0};

		int constraintOffset = g * groupSize;
		int remaining = constraintCount - constraintOffset;
		
		group->constraints = constraints + constraintOffset;
		group->constraintCount = (remaining < groupSize) ? remaining : groupSize;

		int groupSplitBodyOffset = splitBodiesCount;

		for (int i = 0; i < group->constraintCount; ++i) {
			int constraintIndex = constraintOffset + i;
			S2_ASSERT(constraintIndex < constraintCount);

			s2ContactConstraint* constraint = constraints + constraintIndex;
			constraint->indexA = constraint->contact->edges[0].bodyIndex;
			constraint->indexB = constraint->contact->edges[1].bodyIndex;

			s2Body* bodyA = bodies + constraint->indexA;
			s2Body* bodyB = bodies + constraint->indexB;

			int splitIndexA = s2FindSplitBody(bodyA, groupSplitBodyOffset, splitBodies, &splitBodiesCount);
			int splitIndexB = s2FindSplitBody(bodyB, groupSplitBodyOffset, splitBodies, &splitBodiesCount);

			constraint->splitIndexA = splitIndexA;
			constraint->splitIndexB = splitIndexB;
		}
	}

	s2SplitContext sctx;
	sctx.groupsCount = groupsCount;
	sctx.groups = groups;
	sctx.splitBodies = splitBodies;
	sctx.splitBodiesCount = splitBodiesCount;

	return sctx;
}

void s2PrepareContacts_SPLIT_MASSES(s2World* world, s2SplitContext *sctx, bool warmStart)
{
	s2Body* bodies = world->bodies;

	for (int g = 0; g < sctx->groupsCount; g++)
	{
		s2ContactGroup *group = sctx->groups + g;

		for (int i = 0; i < group->constraintCount; ++i)
		{
			s2ContactConstraint* constraint = group->constraints + i;

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

			// this is where the mass splitting happens
			float mA = bodyA->invMass * bodyA->splitCount;
			float iA = bodyA->invI    * bodyA->splitCount;
			float mB = bodyB->invMass * bodyB->splitCount;
			float iB = bodyB->invI    * bodyB->splitCount;

			s2Rot qA = bodyA->rot;
			s2Rot qB = bodyB->rot;

			s2Vec2 normal = constraint->normal;
			s2Vec2 tangent = s2RightPerp(normal);

			for (int j = 0; j < pointCount; ++j)
			{
				const s2ManifoldPoint* mp = manifold->points + j;
				s2ContactConstraintPoint* cp = constraint->points + j;

				cp->localAnchorA = s2Sub(mp->localAnchorA, bodyA->localCenter);
				cp->localAnchorB = s2Sub(mp->localAnchorB, bodyB->localCenter);
			
				s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
				s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);
				cp->rA0 = rA;
				cp->rB0 = rB;

				cp->separation = mp->separation;
				cp->adjustedSeparation = mp->separation - s2Dot(s2Sub(rB, rA), normal);

				cp->biasCoefficient = mp->separation > 0.0f ? 1.0f : 0.0f;

				float rtA = s2Cross(rA, tangent);
				float rtB = s2Cross(rB, tangent);
				float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
				cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

				float rnA = s2Cross(rA, normal);
				float rnB = s2Cross(rB, normal);
				float kNormal = (mA + iA * rnA * rnA) + (mB + iB * rnB * rnB);
				cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;

				float unsplitKNormal = (bodyA->invMass + bodyA->invI * rnA * rnA) + (bodyB->invMass + bodyB->invI * rnB * rnB);
				cp->splitKFactor = unsplitKNormal / kNormal;

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

			}
		}
	}
}

void s2WarmStartContacts_SPLIT_MASSES(s2World* world, s2SplitContext *sctx)
{
	s2Body* bodies = world->bodies;

	for (int g = 0; g < sctx->groupsCount; ++g)
	{
		s2ContactGroup *group = sctx->groups + g;

		for (int i = 0; i < group->constraintCount; ++i)
		{
			s2ContactConstraint* constraint = group->constraints + i;
			
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
			
			s2SplitBody *splitBodyA = sctx->splitBodies + constraint->splitIndexA;
			s2SplitBody *splitBodyB = sctx->splitBodies + constraint->splitIndexB;

			s2Vec2 dvA = splitBodyA->deltaV;
			float  dwA = splitBodyA->deltaW;
			s2Vec2 dvB = splitBodyB->deltaV;
			float  dwB = splitBodyB->deltaW;

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

			splitBodyA->deltaV = dvA;
			splitBodyA->deltaW = dwA;
			splitBodyB->deltaV = dvB;
			splitBodyB->deltaW = dwB;
		}
	}

	s2MergeSplitVelocities(world, sctx);
}



static void s2SolveContacts_SPLIT_MASSES(s2World* world, s2SplitContext *sctx, float inv_h)
{
	s2Body* bodies = world->bodies;

	for (int g = 0; g < sctx->groupsCount; ++g)
	{
		s2ContactGroup *group = sctx->groups + g;

		for (int i = 0; i < group->constraintCount; ++i)
		{
			s2ContactConstraint* constraint = group->constraints + i;
		
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

			s2SplitBody *splitBodyA = sctx->splitBodies + constraint->splitIndexA;
			s2SplitBody *splitBodyB = sctx->splitBodies + constraint->splitIndexB;

			s2Vec2 dvA = splitBodyA->deltaV;
			float  dwA = splitBodyA->deltaW;
			s2Vec2 dvB = splitBodyB->deltaV;
			float  dwB = splitBodyB->deltaW;

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

				float bias = 0.0f;
				if (cp->separation > 0.0f)
				{
					// Speculative
					bias = cp->separation * inv_h;
				}
				else
				{
					// bias must be adjusted by how much the normal-mass changed due to the split
					bias = S2_MAX(cp->splitKFactor * s2_baumgarte * inv_h * S2_MIN(0.0f, cp->separation + s2_linearSlop), -s2_maxBaumgarteVelocity);
				}

				// static anchors
				s2Vec2 rA = cp->rA0;
				s2Vec2 rB = cp->rB0;

				// Relative velocity at contact
				s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
				s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
				float vn = s2Dot(s2Sub(vrB, vrA), normal);

				// Compute normal impulse
				float impulse = -cp->normalMass * (vn + bias);

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
			
				// static anchors
				s2Vec2 rA = cp->rA0;
				s2Vec2 rB = cp->rB0;
			
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

			splitBodyA->deltaV = dvA;
			splitBodyA->deltaW = dwA;
			splitBodyB->deltaV = dvB;
			splitBodyB->deltaW = dwB;
		}
	}

	s2MergeSplitVelocities(world, sctx);
}

void s2Solve_SPLIT_MASSES(s2World* world, s2StepContext* context)
{
	s2Contact* contacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;

	s2Joint* joints = world->joints;
	int jointCapacity = world->jointPool.capacity;

	s2ContactConstraint* constraints =
		s2AllocateStackItem(world->stackAllocator, contactCapacity * sizeof(s2ContactConstraint), "constraint");

	int constraintCount = 0;
	for (int i = 0; i < contactCapacity; ++i)
	{
		s2Contact* contact = contacts + i;
		if (s2IsFree(&contact->object))
		{
			continue;
		}

		if (contact->manifold.pointCount == 0)
		{
			continue;
		}

		constraints[constraintCount].contact = contact;
		constraints[constraintCount].contact->manifold.constraintIndex = constraintCount;
		constraintCount += 1;
	}

	int iterations = context->iterations;
	float h = context->dt;
	float inv_h = context->inv_dt;

	s2IntegrateVelocities(world, h);

	s2SplitContext sctx = s2SplitIntoGroups(world, constraints, constraintCount, context->groupCount);

	s2PrepareContacts_SPLIT_MASSES(world, &sctx, context->warmStart);

	if (context->warmStart)
	{
		s2WarmStartContacts_SPLIT_MASSES(world, &sctx);
	}

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}

		s2PrepareJoint(joint, context, context->warmStart);

		if (context->warmStart)
		{
			s2WarmStartJoint(joint, context);
		}
	}

	for (int iter = 0; iter < iterations; ++iter)
	{
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJoint_Baumgarte(joint, context, h, inv_h, true);
		}

		s2SolveContacts_SPLIT_MASSES(world, &sctx, inv_h);
	}

	// Update positions from velocity
	s2IntegratePositions(world, h);
	s2FinalizePositions(world);

	s2StoreContactImpulses(constraints, constraintCount);

	s2FreeStackItem(world->stackAllocator, sctx.splitBodies);
	s2FreeStackItem(world->stackAllocator, sctx.groups);
	s2FreeStackItem(world->stackAllocator, constraints);
}