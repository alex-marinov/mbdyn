/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "modelns.h"
#include "strnode.h"
#include <matvecexp.h>
#include <Rot.hh>

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

	throw ErrGeneric();
}

/*
 * Computes the position of a structural node
 */
template <IDX_t IDX>
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

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[2]);
	ASSERT(dm != 0);

	unsigned uLabel = unsigned((*arg1)());

	StructNode *pNode = (*dm)()->pFindStructNode(uLabel);
	if (pNode == 0) {
		silent_cerr("position" << IDX2str(IDX)
				<< "(" << uLabel << "): "
				"unable to find StructNode(" << uLabel << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	switch (IDX) {
	case NORM:
		*out = pNode->GetXCurr().Norm();
		break;

	case SQUARE:
		*out = pNode->GetXCurr().Dot();
		break;

	default:
		*out = pNode->GetXCurr()(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the distance between two structural nodes
 */
template <IDX_t IDX>
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

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel1 = unsigned((*arg1)());
	unsigned uLabel2 = unsigned((*arg2)());

	StructNode *pNode1 = (*dm)()->pFindStructNode(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("distance" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	StructNode *pNode2 = (*dm)()->pFindStructNode(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("distance" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	Vec3 d = pNode2->GetXCurr() - pNode1->GetXCurr();

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
 * Computes the orientation between two structural nodes
 */
template <IDX_t IDX>
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

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel1 = unsigned((*arg1)());
	unsigned uLabel2 = unsigned((*arg2)());

	StructNode *pNode1 = (*dm)()->pFindStructNode(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("angle" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	StructNode *pNode2 = (*dm)()->pFindStructNode(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("angle" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	Vec3 phi(RotManip::VecRot(pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()));

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
template <IDX_t IDX>
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

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[2]);
	ASSERT(dm != 0);

	unsigned uLabel = unsigned((*arg1)());

	StructNode *pNode = (*dm)()->pFindStructNode(uLabel);
	if (pNode == 0) {
		silent_cerr("velocity" << IDX2str(IDX)
				<< "(" << uLabel << "): "
				"unable to find StructNode(" << uLabel << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	switch (IDX) {
	case NORM:
		*out = pNode->GetVCurr().Norm();
		break;

	case SQUARE:
		*out = pNode->GetVCurr().Dot();
		break;

	default:
		*out = pNode->GetVCurr()(IDX);
		break;
	}

	return 0;
}

/*
 * Computes the relative velocity between two structural nodes
 */
template <IDX_t IDX>
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

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(arg2 != 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel1 = unsigned((*arg1)());
	unsigned uLabel2 = unsigned((*arg2)());

	StructNode *pNode1 = (*dm)()->pFindStructNode(uLabel1);
	if (pNode1 == 0) {
		silent_cerr("vrel" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel1 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	StructNode *pNode2 = (*dm)()->pFindStructNode(uLabel2);
	if (pNode2 == 0) {
		silent_cerr("vrel" << IDX2str(IDX)
				<< "(" << uLabel1 << "," << uLabel2 << "): "
				"unable to find StructNode(" << uLabel2 << ")"
				<< std::endl);
		throw ErrGeneric();
	}

	Vec3 d = pNode2->GetVCurr() - pNode1->GetVCurr();

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

	MathParser::MathArgReal_t *arg2 = dynamic_cast<MathParser::MathArgReal_t *>(args[2]);
	ASSERT(arg2 != 0);

	ModelNameSpace::MathArgDM *dm = dynamic_cast<ModelNameSpace::MathArgDM *>(args[3]);
	ASSERT(dm != 0);

	unsigned uLabel = unsigned((*arg1)());

	const DriveCaller *pDC = (*dm)()->GetMBDynParser().GetDrive(uLabel);
	if (pDC == 0) {
		silent_cerr("model namespace: drive " << uLabel << " not available" << std::endl);
		throw ErrGeneric();
	}

	*out = pDC->dGet((*arg2)());

	return 0;
}

/*
 * Computes the value of a scalar function
 */
static int
model_sf(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);
	ASSERT(args[2]->Type() == MathParser::AT_INT);
	ASSERT(args[3]->Type() == MathParser::AT_PRIVATE);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);
	doublereal v = (*arg1)();

	MathParser::MathArgInt_t *arg2 = dynamic_cast<MathParser::MathArgInt_t*>(args[2]);
	ASSERT(arg2 != 0);
	int order = (*arg2)();
	ASSERT(order >= 0);

	ModelNameSpace::MathArgSF *sf = dynamic_cast<ModelNameSpace::MathArgSF*>(args[3]);
	ASSERT(sf != 0);

	if (order == 0) {
		*out = (*(*sf)())(v);

	} else {
		const DifferentiableScalarFunction *dsf = dynamic_cast<const DifferentiableScalarFunction *>((*sf)());
		if (dsf == 0) {
			silent_cerr("model namespace scalar function: "
				"order=" << order << " only allowed "
				"with differentiable scalar functions"
				<< std::endl);
			throw ErrGeneric();
		}

		*out = dsf->ComputeDiff(v, order);
	}

	return 0;
}

ModelNameSpace::ModelNameSpace(DataManager *pdm)
: MathParser::NameSpace("model"), pDM(pdm)
{
	MathParser::MathFunc_t *f;

	// position
	f = new MathParser::MathFunc_t;
	f->fname = "position";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<NORM>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// position2
	f = new MathParser::MathFunc_t;
	f->fname = "position2";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<SQUARE>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// xposition
	f = new MathParser::MathFunc_t;
	f->fname = "xposition";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX1>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// yposition
	f = new MathParser::MathFunc_t;
	f->fname = "yposition";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX2>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// zposition
	f = new MathParser::MathFunc_t;
	f->fname = "zposition";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = position<IDX3>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// distance
	f = new MathParser::MathFunc_t;
	f->fname = "distance";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<NORM>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// distance2
	f = new MathParser::MathFunc_t;
	f->fname = "distance2";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<SQUARE>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// xdistance
	f = new MathParser::MathFunc_t;
	f->fname = "xdistance";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX1>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// ydistance
	f = new MathParser::MathFunc_t;
	f->fname = "ydistance";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX2>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// zdistance
	f = new MathParser::MathFunc_t;
	f->fname = "zdistance";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = distance<IDX3>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// anglerel
	f = new MathParser::MathFunc_t;
	f->fname = "anglerel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<NORM>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// xanglerel
	f = new MathParser::MathFunc_t;
	f->fname = "xanglerel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX1>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// yanglerel
	f = new MathParser::MathFunc_t;
	f->fname = "yanglerel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX2>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// zanglerel
	f = new MathParser::MathFunc_t;
	f->fname = "zanglerel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = anglerel<IDX3>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// velocity
	f = new MathParser::MathFunc_t;
	f->fname = "velocity";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<NORM>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// velocity2
	f = new MathParser::MathFunc_t;
	f->fname = "velocity2";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<SQUARE>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// xvelocity
	f = new MathParser::MathFunc_t;
	f->fname = "xvelocity";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX1>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// yvelocity
	f = new MathParser::MathFunc_t;
	f->fname = "yvelocity";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX2>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// zvelocity
	f = new MathParser::MathFunc_t;
	f->fname = "zvelocity";
	f->args.resize(1 + 1 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathArgDM;
	f->f = velocity<IDX3>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// vrel
	f = new MathParser::MathFunc_t;
	f->fname = "vrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<NORM>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// vrel2
	f = new MathParser::MathFunc_t;
	f->fname = "vrel2";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<SQUARE>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// xvrel
	f = new MathParser::MathFunc_t;
	f->fname = "xvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX1>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// yvrel
	f = new MathParser::MathFunc_t;
	f->fname = "yvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX2>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// zvrel
	f = new MathParser::MathFunc_t;
	f->fname = "zvrel";
	f->args.resize(1 + 2 + 1);
	f->args[0] = new MathParser::MathArgReal_t;
	f->args[1] = new MathParser::MathArgInt_t;
	f->args[2] = new MathParser::MathArgInt_t;
	f->args[3] = new MathArgDM;
	f->f = vrel<IDX3>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("model namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
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
		throw ErrGeneric();
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
}

ModelNameSpace::~ModelNameSpace(void)
{
	for (funcType::iterator f = func.begin(); f != func.end(); f++) {
		for (MathParser::MathArgs::iterator i = f->second->args.begin();
			i != f->second->args.end();
			i++)
		{
			delete *i;
		}

		delete f->second;
	}

	for (MathParser::MathArgs::iterator i = sf_func.args.begin();
		i != sf_func.args.end();
		i++)
	{
		delete *i;
	}
}

bool
ModelNameSpace::IsFunc(const char* const s) const
{
	return GetFunc(s) != 0;
}

MathParser::MathFunc_t*
ModelNameSpace::GetFunc(const char* const s) const
{
	if (strncmp(s, "sf::", STRLENOF("sf::")) == 0) {
		const BasicScalarFunction *sf = pDM->GetMBDynParser().GetScalarFunction(&s[STRLENOF("sf::")]);

		if (sf == 0) {
			throw ErrGeneric();
		}

		(*dynamic_cast<MathArgSF*>(sf_func.args[3]))() = sf;

		return const_cast<MathParser::MathFunc_t*>(&sf_func);
	}

	funcType::const_iterator i = func.find(std::string(s));

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
		return TypedValue((*dynamic_cast<MathParser::MathArgInt_t*>(args[0]))());

	case MathParser::AT_REAL:
		return TypedValue((*dynamic_cast<MathParser::MathArgReal_t*>(args[0]))());

	default:
		throw ErrGeneric();
	}
}

