/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 18:19:39 +0530 (Thu, 16 May 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of OrthoCompactionModule which models an
 *  interface for orthogonal compaction approaches.
 *
 * \author Karsten Klein
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#ifdef _MSC_VER
#pragma once
#endif


#ifndef OGDF_ORTHO_COMPACTION_MODULE_H
#define OGDF_ORTHO_COMPACTION_MODULE_H


#include <ogdf/basic/Module.h>
#include <ogdf/orthogonal/OrthoRep.h>
#include <ogdf/planarity/PlanRep.h>

namespace ogdf {

class OGDF_EXPORT OrthoCompactionModule : public Module{

public:
	//! Constructs the compaction module
	OrthoCompactionModule() {}
	//! Destroys the compaction module
	virtual ~OrthoCompactionModule() {}

	//! call to construct a drawing for an orthogonal representation \a OR of PlanRep \a PG.
	//! Has to be implemented by derived classes.
	virtual void callConstructive(
		PlanRep &PG,
		OrthoRep &OR,
		GridLayout &drawing) = 0;


	//! call to improve a given orthogonal drawing \a drawing. Has
	//! to be implemented by derived classes.
	virtual void callImprovement(
		PlanRep &PG,
		OrthoRep &OR,
		GridLayout &drawing) = 0;
};

} // end namespace ogdf


#endif
