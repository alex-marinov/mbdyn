/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cfloat>
#include <limits>

#include "modelns.h"
#include "strnode.h"
#include "matvecexp.h"
#include "Rot.hh"

enum IDX_t {
	IDX1 = 1,
	IDX2 = 2,
	IDX3 = 3,
	NORM,
	SQUARE
};

static const char *
IDX2str(IDX_t IDX)
{
	switch (IDX) {
	case NORM:
		return "";

	case SQUARE:
		return "2";

	case IDX1:
		return "<1>";

	case IDX2:
		return "<2>";

	case IDX3:
		return "<3>";
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

enum when_t {
	CURR = 0,
	PREV = 1
};

static const char *
when2str(when_t when)
{
	switch (when) {
	case CURR:
		return "";

	case PREV:
		return "_prev";

	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/*
 * Computes the position of a structural node
 */
template <IDX_t IDX, when_t when>
static int
position(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[2]);
	ASSERT(dm != 0);

	unsigned uLabel = unsigned((*arg1)());

	const StructNode *pNode = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel);
	if (pNode == 0) {
		silent_cerr("model::position" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel << "): "
				"unable to find StructNode(" << uLabel << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 x;
	switch (when) {
	case CURR:
		x = pNode->GetXCurr();
		break;

	case PREV:
		x = pNode->GetXPrev();
		break;
	}

	switch (IDX) {
	case NORM:
		*out = x.Norm();
		break;

	case SQUARE:
		*out = x.Dot();
		break;

	default:
		*out = x(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the distance between two structural nodes
 */
template <IDX_t IDX, when_t when>
static int
distance(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_INT);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);
	ASSERT((*arg2)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel1 = unsigned((*arg1)());
	unsigned uLabel2 = unsigned((*arg2)());

	const StructNode *pNode1 =  (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("model::distance" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const StructNode *pNode2 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("model::distance" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 d;
	switch (when) {
	case CURR:
		d = pNode2->GetXCurr() - pNode1->GetXCurr();
		break;

	case PREV:
		d = pNode2->GetXPrev() - pNode1->GetXPrev();
		break;
	}

	switch (IDX) {
	case NORM:
		*out = d.Norm();
		break;

	case SQUARE:
		*out = d.Dot();
		break;

	default:
		*out = d(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the components of the direction between two structural nodes
 */
template <IDX_t IDX, when_t when>
static int
unitvec(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_INT);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);
	ASSERT((*arg2)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel1 = unsigned((*arg1)());
	unsigned uLabel2 = unsigned((*arg2)());

	const StructNode *pNode1 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("model::unitvec" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const StructNode *pNode2 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("model::unitvec" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 d;
	switch (when) {
	case CURR:
		d = pNode2->GetXCurr() - pNode1->GetXCurr();
		break;

	case PREV:
		d = pNode2->GetXPrev() - pNode1->GetXPrev();
		break;
	}

	doublereal dd = d.Norm();
	if (dd <= std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("model::unitvec" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"null distance"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	switch (IDX) {
	default:
		*out = d(IDX)/dd;
		break;
	}

	return 0;
}

/*
 * Computes the orientation of a structural node
 */
template <IDX_t IDX, when_t when>
static int
angle(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[2]);
	ASSERT(dm != 0);

	unsigned uLabel = unsigned((*arg1)());

	const StructNode *pNode = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel);
	if (pNode == 0) {
		silent_cerr("model::angle" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel << "): "
				"unable to find StructNode(" << uLabel << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 phi;
	switch (when) {
	case CURR:
		phi = RotManip::VecRot(pNode->GetRCurr());
		break;

	case PREV:
		phi = RotManip::VecRot(pNode->GetRPrev());
		break;
	}

	switch (IDX) {
	case NORM:
		*out = phi.Norm();
		break;

	default:
		*out = phi(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the orientation between two structural nodes
 */
template <IDX_t IDX, when_t when>
static int
anglerel(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_INT);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);
	ASSERT((*arg2)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel1 = unsigned((*arg1)());
	unsigned uLabel2 = unsigned((*arg2)());

	const StructNode *pNode1 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("model::anglerel" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const StructNode *pNode2 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("model::anglerel" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 phi;
	switch (when) {
	case CURR:
		phi = RotManip::VecRot(pNode1->GetRCurr().MulTM(pNode2->GetRCurr()));
		break;

	case PREV:
		phi = RotManip::VecRot(pNode1->GetRPrev().MulTM(pNode2->GetRPrev()));
		break;
	}

	switch (IDX) {
	case NORM:
		*out = phi.Norm();
		break;

	default:
		*out = phi(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the velocity of a structural node
 */
template <IDX_t IDX, when_t when>
static int
velocity(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[2]);
	ASSERT(dm != 0);

	unsigned uLabel = unsigned((*arg1)());

	const StructNode *pNode = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel);
	if (pNode == 0) {
		silent_cerr("model::velocity" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel << "): "
				"unable to find StructNode(" << uLabel << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 v;
	switch (when) {
	case CURR:
		v = pNode->GetVCurr();
		break;
	
	case PREV:
		v = pNode->GetVPrev();
		break;
	}
	
	switch (IDX) {
	case NORM:
		*out = v.Norm();
		break;

	case SQUARE:
		*out = v.Dot();
		break;

	default:
		*out = v(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the relative velocity between two structural nodes
 */
template <IDX_t IDX, when_t when>
static int
vrel(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_INT);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);
	ASSERT((*arg2)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel1 = unsigned((*arg1)());
	unsigned uLabel2 = unsigned((*arg2)());

	const StructNode *pNode1 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("model::vrel" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const StructNode *pNode2 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("model::vrel" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 v;
	switch (when) {
	case CURR:
		v = pNode2->GetVCurr() - pNode1->GetVCurr();
		break;

	case PREV:
		v = pNode2->GetVCurr() - pNode1->GetVPrev();
		break;
	}

	switch (IDX) {
	case NORM:
		*out = v.Norm();
		break;

	case SQUARE:
		*out = v.Dot();
		break;

	default:
		*out = v(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the angular velocity of a structural node
 */
template <IDX_t IDX, when_t when>
static int
angvel(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[2]);
	ASSERT(dm != 0);

	unsigned uLabel = unsigned((*arg1)());

	const StructNode *pNode = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel);
	if (pNode == 0) {
		silent_cerr("model::angvel" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel << "): "
				"unable to find StructNode(" << uLabel << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 w;
	switch (when) {
	case CURR:
		w = pNode->GetWCurr();
		break;

	case PREV:
		w = pNode->GetWPrev();
		break;
	}

	switch (IDX) {
	case NORM:
		*out = w.Norm();
		break;

	case SQUARE:
		*out = w.Dot();
		break;

	default:
		*out = w(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the relative angular velocity between two structural nodes
 */
template <IDX_t IDX, when_t when>
static int
angvrel(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_INT);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);
	ASSERT((*arg2)() >= 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel1 = unsigned((*arg1)());
	unsigned uLabel2 = unsigned((*arg2)());

	const StructNode *pNode1 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("model::angvrel" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const StructNode *pNode2 = (*dm)()->pFindNode<const StructNode, Node::STRUCTURAL>(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("model::angvrel" << when2str(when) << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 w;
	switch (when) {
	case CURR:
		w = pNode2->GetWCurr() - pNode1->GetWCurr();
		break;

	case PREV:
		w = pNode2->GetWPrev() - pNode1->GetWPrev();
		break;
	}

	switch (IDX) {
	case NORM:
		*out = w.Norm();
		break;

	case SQUARE:
		*out = w.Dot();
		break;

	default:
		*out = w(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the value of a drive
 */
static int
drive(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_REAL);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);
	ASSERT((*arg1)() >= 0);

	MathParser::MathArgReal_t *arg2 = dynamic_cast<MathParser::MathArgReal_t *>(args[2]);
	ASSERT(arg2 != 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel = unsigned((*arg1)());

	const DriveCaller *pDC = (*dm)()->GetMBDynParser().GetDrive(uLabel);
	if (pDC == 0) {
		silent_cerr("model::drive(" << uLabel << ") not available"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal val = (*arg2)();
	*out = pDC->dGet(val);

	return 0;
}

/*
 * Computes the value of a scalar function
 */
static int
model_sf(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);
	ASSERT(args[2]->Type() == MathParser::AT_INT);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);
	doublereal v = (*arg1)();

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);
	int order = (*arg2)();
	ASSERT(order >= 0);

	ModelNameSpace::MathArgSF *sf = dynamic_cast<ModelNameSpace::MathArgSF *>(args[3]);
	ASSERT(sf != 0);

	if (order < 0 || order > 1) {
		silent_cerr("model::sf invalid derivative order " << order << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (order == 0) {
		*out = (*(*sf)())(v);

	} else {
		const DifferentiableScalarFunction *dsf = dynamic_cast<const DifferentiableScalarFunction *>((*sf)());
		if (dsf == 0) {
			silent_cerr("model::sf "
				"order=" << order << " only allowed "
				"with differentiable scalar functions"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		*out = dsf->ComputeDiff(v, order);
	}

	return 0;
}

/*
 * Computes the value of a node's private data
 */
static int
model_node(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 2);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_STRING);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);
	ASSERT(args[4]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgInt_t *arg_label = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg_label != 0);
	int iLabel = (*arg_label)();

	MathParser::MathArgString_t *arg_val = dynamic_cast<MathParser::MathArgString_t *>(args[2]);
	ASSERT(arg_val != 0);
	std::string v = (*arg_val)();

	ModelNameSpace::MathArgNode *arg_type = dynamic_cast<ModelNameSpace::MathArgNode *>(args[3]);
	ASSERT(arg_type != 0);
	Node::Type type = (*arg_type)();

	ModelNameSpace::MathArgDM *arg_dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[4]);
	ASSERT(arg_dm != 0);
	const DataManager *pDM = (*arg_dm)();

	if (iLabel < 0) {
		silent_cerr("model::node: invalid node label " << iLabel << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	unsigned uLabel = unsigned(iLabel);

	const Node *pNode = pDM->pFindNode(type, uLabel);
	if (pNode == 0) {
		silent_cerr("model::node: unable to find " << psNodeNames[type] << "(" << uLabel << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	int idx = pNode->iGetPrivDataIdx(v.c_str());
	if (idx == 0) {
		silent_cerr("model::node: " << psNodeNames[type] << "(" << pNode->GetLabel() << "): invalid private data \"" << v << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	*out = pNode->dGetPrivData(idx);

	return 0;
}

/*
 * Computes the value of an element's private data
 */
template <bool unique>
static int
model_elem(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 - unique + 2);

	int idx = 0;
	ASSERT(args[idx]->Type() == MathParser::AT_REAL);
	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[idx]);
	ASSERT(out != 0);

	unsigned uLabel(-1);
	if (!unique) {
		++idx;
		ASSERT(args[idx]->Type() == MathParser::AT_INT);
		MathParser::MathArgInt_t *arg_label = dynamic_cast<MathParser::MathArgInt_t *>(args[idx]);
		ASSERT(arg_label != 0);
		int iLabel = (*arg_label)();
		if (iLabel < 0) {
			silent_cerr("model::element: invalid element label " << iLabel << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		uLabel = unsigned(iLabel);
	}

	++idx;
	ASSERT(args[idx]->Type() == MathParser::AT_STRING);
	MathParser::MathArgString_t *arg_val = dynamic_cast<MathParser::MathArgString_t *>(args[idx]);
	ASSERT(arg_val != 0);
	std::string v = (*arg_val)();

	++idx;
	ASSERT(args[idx]->Type() == MathParser::AT_PRIVATE);
	ModelNameSpace::MathArgElem *arg_elem = dynamic_cast<ModelNameSpace::MathArgElem *>(args[idx]);
	ASSERT(arg_elem != 0);
	Elem::Type type = (*arg_elem)();

	++idx;
	ASSERT(args[idx]->Type() == MathParser::AT_PRIVATE);
	ModelNameSpace::MathArgDM *arg_dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[idx]);
	ASSERT(arg_dm != 0);
	const DataManager *pDM = (*arg_dm)();

	const Elem *pElem = pDM->pFindElem(type, uLabel);
	if (pElem == 0) {
		silent_cerr("model::element: unable to find " << psElemNames[type] << "(" << uLabel << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	idx = pElem->iGetPrivDataIdx(v.c_str());
	if (idx == 0) {
		silent_cerr("model::element: " << psElemNames[type] << "(" << pElem->GetLabel() << "): invalid private data \"" << v << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	*out = pElem->dGetPrivData(idx);

	return 0;
}

/*
 * Computes the value of the current simulation entity's private data
 */
static int
model_curr(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1 + 2);

	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	ASSERT(args[1]->Type() == MathParser::AT_STRING);
	MathParser::MathArgString_t *arg_name = dynamic_cast<MathParser::MathArgString_t *>(args[1]);
	ASSERT(arg_name != 0);
	const std::string& name = (*arg_name)();

	ASSERT(args[2]->Type() == MathParser::AT_PRIVATE);
	ModelNameSpace::MathArgMNS *arg_mns = dynamic_cast<ModelNameSpace::MathArgMNS *>(args[2]);
	ASSERT(arg_mns != 0);
	const ModelNameSpace *pMNS = (*arg_mns)();

	TypedValue value;
	if (!pMNS->GetCurrData(name, value)) {
		silent_cerr("model::current(" << name << ") not defined" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	*out = value.GetReal();

	return 0;
}

ModelNameSpace::ModelNameSpace(const DataManager *pDM)
: MathParser::NameSpace("model"), pDM(pDM)
{
	MathParser::MathFunc_t *f;

	// position
	f = new MathParser::MathFunc_t;
	f->fname = "position";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<NORM, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// position2
	f = new MathParser::MathFunc_t;
	f->fname = "position2";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<SQUARE, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xposition
	f = new MathParser::MathFunc_t;
	f->fname = "xposition";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yposition
	f = new MathParser::MathFunc_t;
	f->fname = "yposition";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zposition
	f = new MathParser::MathFunc_t;
	f->fname = "zposition";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// distance
	f = new MathParser::MathFunc_t;
	f->fname = "distance";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<NORM, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// distance2
	f = new MathParser::MathFunc_t;
	f->fname = "distance2";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<SQUARE, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xdistance
	f = new MathParser::MathFunc_t;
	f->fname = "xdistance";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// ydistance
	f = new MathParser::MathFunc_t;
	f->fname = "ydistance";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zdistance
	f = new MathParser::MathFunc_t;
	f->fname = "zdistance";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xunitvec
	f = new MathParser::MathFunc_t;
	f->fname = "xunitvec";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = unitvec<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yunitvec
	f = new MathParser::MathFunc_t;
	f->fname = "yunitvec";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = unitvec<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zunitvec
	f = new MathParser::MathFunc_t;
	f->fname = "zunitvec";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = unitvec<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angle
	f = new MathParser::MathFunc_t;
	f->fname = "angle";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angle<NORM, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xangle
	f = new MathParser::MathFunc_t;
	f->fname = "xangle";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angle<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yangle
	f = new MathParser::MathFunc_t;
	f->fname = "yangle";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angle<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zangle
	f = new MathParser::MathFunc_t;
	f->fname = "zangle";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angle<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// anglerel
	f = new MathParser::MathFunc_t;
	f->fname = "anglerel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<NORM, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xanglerel
	f = new MathParser::MathFunc_t;
	f->fname = "xanglerel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yanglerel
	f = new MathParser::MathFunc_t;
	f->fname = "yanglerel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zanglerel
	f = new MathParser::MathFunc_t;
	f->fname = "zanglerel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// velocity
	f = new MathParser::MathFunc_t;
	f->fname = "velocity";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<NORM, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// velocity2
	f = new MathParser::MathFunc_t;
	f->fname = "velocity2";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<SQUARE, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xvelocity
	f = new MathParser::MathFunc_t;
	f->fname = "xvelocity";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yvelocity
	f = new MathParser::MathFunc_t;
	f->fname = "yvelocity";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zvelocity
	f = new MathParser::MathFunc_t;
	f->fname = "zvelocity";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// vrel
	f = new MathParser::MathFunc_t;
	f->fname = "vrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<NORM, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// vrel2
	f = new MathParser::MathFunc_t;
	f->fname = "vrel2";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<SQUARE, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xvrel
	f = new MathParser::MathFunc_t;
	f->fname = "xvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yvrel
	f = new MathParser::MathFunc_t;
	f->fname = "yvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zvrel
	f = new MathParser::MathFunc_t;
	f->fname = "zvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angvel
	f = new MathParser::MathFunc_t;
	f->fname = "angvel";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<NORM, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angvel
	f = new MathParser::MathFunc_t;
	f->fname = "angvel2";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<SQUARE, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xangvel
	f = new MathParser::MathFunc_t;
	f->fname = "xangvel";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yangvel
	f = new MathParser::MathFunc_t;
	f->fname = "yangvel";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zangvel
	f = new MathParser::MathFunc_t;
	f->fname = "zangvel";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angvrel
	f = new MathParser::MathFunc_t;
	f->fname = "angvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<NORM, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angvrel2
	f = new MathParser::MathFunc_t;
	f->fname = "angvrel2";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<SQUARE, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xangvrel
	f = new MathParser::MathFunc_t;
	f->fname = "xangvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<IDX1, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yangvrel
	f = new MathParser::MathFunc_t;
	f->fname = "yangvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<IDX2, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zangvrel
	f = new MathParser::MathFunc_t;
	f->fname = "zangvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<IDX3, CURR>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// position_prev
	f = new MathParser::MathFunc_t;
	f->fname = "position_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<NORM, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// position2_prev
	f = new MathParser::MathFunc_t;
	f->fname = "position2_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<SQUARE, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xposition_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xposition_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yposition_prev
	f = new MathParser::MathFunc_t;
	f->fname = "yposition_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zposition_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zposition_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// distance_prev
	f = new MathParser::MathFunc_t;
	f->fname = "distance_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<NORM, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// distance2_prev
	f = new MathParser::MathFunc_t;
	f->fname = "distance2_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<SQUARE, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xdistance_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xdistance_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// ydistance_prev
	f = new MathParser::MathFunc_t;
	f->fname = "ydistance_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zdistance_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zdistance_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xunitvec_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xunitvec_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = unitvec<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yunitvec_prev
	f = new MathParser::MathFunc_t;
	f->fname = "yunitvec_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = unitvec<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zunitvec_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zunitvec_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = unitvec<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angle_prev
	f = new MathParser::MathFunc_t;
	f->fname = "angle_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angle<NORM, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xangle_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xangle_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angle<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yangle_prev
	f = new MathParser::MathFunc_t;
	f->fname = "yangle_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angle<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zangle_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zangle_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angle<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// anglerel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "anglerel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<NORM, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xanglerel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xanglerel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yanglerel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "yanglerel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zanglerel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zanglerel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// velocity_prev
	f = new MathParser::MathFunc_t;
	f->fname = "velocity_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<NORM, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// velocity2_prev
	f = new MathParser::MathFunc_t;
	f->fname = "velocity2_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<SQUARE, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xvelocity_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xvelocity_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yvelocity_prev
	f = new MathParser::MathFunc_t;
	f->fname = "yvelocity_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zvelocity_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zvelocity_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// vrel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "vrel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<NORM, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// vrel2_prev
	f = new MathParser::MathFunc_t;
	f->fname = "vrel2_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<SQUARE, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xvrel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xvrel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yvrel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "yvrel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zvrel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zvrel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angvel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "angvel_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<NORM, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angvel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "angvel2_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<SQUARE, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xangvel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xangvel_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yangvel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "yangvel_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zangvel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zangvel_prev";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = angvel<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angvrel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "angvrel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<NORM, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// angvrel2_prev
	f = new MathParser::MathFunc_t;
	f->fname = "angvrel2_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<SQUARE, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// xangvrel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "xangvrel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<IDX1, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// yangvrel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "yangvrel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<IDX2, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// zangvrel_prev
	f = new MathParser::MathFunc_t;
	f->fname = "zangvrel_prev";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = angvrel<IDX3, PREV>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// drive
	f = new MathParser::MathFunc_t;
	f->fname = "drive";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgReal_t;
	f->args[3] = new MathArgDM;
	f->f = drive;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// current simulation entity function
	f = new MathParser::MathFunc_t;
	f->fname = "current";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgString_t;
	f->args[2] = new MathArgMNS(this);
	f->f = model_curr;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// scalar functions
	sf_func.fname = "sf";
	sf_func.args.resize(1 + 2 + 1);
	sf_func.args[0] = new MathParser::MathArgReal_t;
	sf_func.args[1] = new MathParser::MathArgReal_t;
	sf_func.args[2] = new MathParser::MathArgInt_t(0, MathParser::AF_OPTIONAL);
	sf_func.args[3] = new MathArgSF;
	sf_func.f = model_sf;
	sf_func.t = 0;

	// node functions
	node_func.fname = "node";
	node_func.args.resize(1 + 2 + 2);
	node_func.args[0] = new MathParser::MathArgReal_t;
	node_func.args[1] = new MathParser::MathArgInt_t;
	node_func.args[2] = new MathParser::MathArgString_t;
	node_func.args[3] = new MathArgNode;
	node_func.args[4] = new MathArgDM;
	node_func.f = model_node;
	node_func.t = 0;

	// element functions
	elem_func.fname = "element";
	elem_func.args.resize(1 + 2 + 2);
	elem_func.args[0] = new MathParser::MathArgReal_t;
	elem_func.args[1] = new MathParser::MathArgInt_t;
	elem_func.args[2] = new MathParser::MathArgString_t;
	elem_func.args[3] = new MathArgElem;
	elem_func.args[4] = new MathArgDM;
	elem_func.f = model_elem<false>;
	elem_func.t = 0;

	// unique element functions
	unique_elem_func.fname = "uniqueElement";
	unique_elem_func.args.resize(1 + 1 + 2);
	unique_elem_func.args[0] = new MathParser::MathArgReal_t;
	unique_elem_func.args[1] = new MathParser::MathArgString_t;
	unique_elem_func.args[2] = new MathArgElem;
	unique_elem_func.args[3] = new MathArgDM;
	unique_elem_func.f = model_elem<true>;
	unique_elem_func.t = 0;
}

ModelNameSpace::~ModelNameSpace(void)
{
	for (funcType::iterator f = func.begin(); f != func.end(); ++f) {
		for (MathParser::MathArgs::iterator i = f->second->args.begin();
			i != f->second->args.end(); ++i)
		{
			delete *i;
		}

		delete f->second;
	}

	for (MathParser::MathArgs::iterator i = sf_func.args.begin();
		i != sf_func.args.end(); ++i)
	{
		delete *i;
	}

	for (MathParser::MathArgs::iterator i = node_func.args.begin();
		i != node_func.args.end(); ++i)
	{
		delete *i;
	}

	for (MathParser::MathArgs::iterator i = elem_func.args.begin();
		i != elem_func.args.end(); ++i)
	{
		delete *i;
	}

	for (MathParser::MathArgs::iterator i = unique_elem_func.args.begin();
		i != unique_elem_func.args.end(); ++i)
	{
		delete *i;
	}
}

bool
ModelNameSpace::IsFunc(const std::string& fname) const
{
	return GetFunc(fname) != 0;
}

MathParser::MathFunc_t*
ModelNameSpace::GetFunc(const std::string& fname) const
{
	static const std::string sf_prefix("sf::");
	if (fname.compare(0, sf_prefix.size(), sf_prefix) == 0) {
		std::string sfname(fname, sf_prefix.size());
		const BasicScalarFunction *sf = pDM->GetMBDynParser().GetScalarFunction(sfname);

		if (sf == 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		(*dynamic_cast<MathArgSF *>(sf_func.args[3]))() = sf;

		return const_cast<MathParser::MathFunc_t *>(&sf_func);
	}

	static const std::string node_prefix = "node::";
	if (fname.compare(0, node_prefix.size(), node_prefix) == 0) {
		std::string type(fname, node_prefix.size());
		const Node::Type t = str2nodetype(type.c_str());
		if (t == Node::UNKNOWN) {
			silent_cerr("ModelNameSpace::GetFunc(" << fname << "): unable to find node type \"" << type << "\"" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		(*dynamic_cast<MathParser::MathArgInt_t *>(node_func.args[1]))() = -1;
		(*dynamic_cast<MathArgNode *>(node_func.args[3]))() = t;

		return const_cast<MathParser::MathFunc_t *>(&node_func);
	}

	static const std::string elem_prefix = "element::";
	if (fname.compare(0, elem_prefix.size(), elem_prefix) == 0) {
		std::string type(fname, elem_prefix.size());
		Elem::Type t = str2elemtype(type.c_str());
		if (t == Elem::UNKNOWN) {
			silent_cerr("ModelNameSpace::GetFunc(" << fname << "): unable to find element type \"" << type << "\"" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (pDM->GetElemDataStructure(t).bIsUnique()) {
			(*dynamic_cast<MathArgElem *>(unique_elem_func.args[2]))() = t;

			return const_cast<MathParser::MathFunc_t *>(&unique_elem_func);

		} else {
			(*dynamic_cast<MathParser::MathArgInt_t *>(elem_func.args[1]))() = -1;
			(*dynamic_cast<MathArgElem *>(elem_func.args[3]))() = t;

			return const_cast<MathParser::MathFunc_t *>(&elem_func);
		}
	}

	funcType::const_iterator i = func.find(fname);

	if (i != func.end()) {
		return i->second;
	}

	return 0;
}

TypedValue 
ModelNameSpace::EvalFunc(MathParser::MathFunc_t *f, const MathParser::MathArgs& args) const
{
	for (unsigned i = 0; i != args.size(); i++) {
		if (args[i]->Type() == MathParser::AT_PRIVATE) {
			MathArgDM *dm = dynamic_cast<MathArgDM *>(args[i]);
			if (dm) {
				(*dm)() = pDM;
			}
		}
	}

	f->f(args);

	switch (args[0]->Type()) {
	case MathParser::AT_VOID:
		return TypedValue(0);

	case MathParser::AT_INT:
		return TypedValue((*dynamic_cast<MathParser::MathArgInt_t *>(args[0]))());

	case MathParser::AT_REAL:
		return TypedValue((*dynamic_cast<MathParser::MathArgReal_t *>(args[0]))());

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

bool
ModelNameSpace::PushCurrData(const std::string& name, const TypedValue& value)
{
	return currData.insert(currDataType::value_type(name, value)).second;
}

bool
ModelNameSpace::PopCurrData(const std::string& name)
{
	currDataType::iterator i = currData.find(name);
	if (i == currData.end()) {
		return false;
	}

	currData.erase(i);

	return true;
}

bool
ModelNameSpace::GetCurrData(const std::string& name, TypedValue& value) const
{
	currDataType::const_iterator i = currData.find(name);
	if (i == currData.end()) {
		return false;
	}

	value = i->second;

	return true;
}

