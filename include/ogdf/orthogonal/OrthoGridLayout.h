/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 18:19:39 +0530 (Thu, 16 May 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class OrthoGridLayout which represents an
 *        orthogonal planar grid drawing algorithm.
 *
 * \author Carsten Gutwenger, Sebastian Leipert, Karsten Klein
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


#ifndef OGDF_ORTHO_GRID_LAYOUT_H
#define OGDF_ORTHO_GRID_LAYOUT_H


#include <ogdf/module/GridLayoutModule.h>
#include <ogdf/orthogonal/OrthoRep.h>
#include <ogdf/basic/ModuleOption.h>
#include <ogdf/module/EmbedderModule.h>
#include <ogdf/module/OrthoCompactionModule.h>


namespace ogdf {


	enum OrthoDir;


//---------------------------------------------------------
// OrthoGridLayout
// represents planar orthogonal grid drawing algorithm
//---------------------------------------------------------
class OGDF_EXPORT OrthoGridLayout : public GridLayoutPlanRepModule
{
public:
	// constructor
	OrthoGridLayout();


	// calls planar grid layout algorithm. Input is a planarized representation
	// PG of a connected component of the graph, output is a layout of the
	// (modified) planarized representation in drawing
	//void call(PlanRep &PG, adjEntry adjExternal, Layout &drawing);

	//
	// options

	// the distance from the tight bounding box to the boundary of the drawing
	double margin() const {
		return m_margin;
	}

	void margin(int m) {
		m_margin = m;
	}

	//! Set the option profile, thereby fixing a set of drawing options
	void optionProfile(int i) {m_optionProfile = i;}

	//! Set scaling compaction
	//void scaling(bool b) {m_useScalingCompaction = b;}

	//! Set bound on the number of bends
	void setBendBound(int i) { OGDF_ASSERT(i >= 0); m_bendBound = i; }

	//in planarlayout
	//enum LayoutOptions {optScaling = 2, optProgressive = 4}
	//set generic options by setting field bits,
	//necessary to allow setting over base class pointer
	//bit 0 = not used so far
	//bit 1 = scaling
	//bit 2 = progressive/traditional
	//=> 0 is standard
	virtual void setOptions(int optionField)
	{
		//if (optionField & umlOpScale) m_useScalingCompaction = true;
		//else m_useScalingCompaction = false;
		if (optionField & umlOpProg) m_orthoStyle = 1;
		else m_orthoStyle = 0; //traditional
	}//setOptions

	virtual int getOptions()
	{
		int result = 0;
		//if (m_useScalingCompaction) result += umlOpScale;
		if (m_orthoStyle == 1) result += umlOpProg;

		return result;
	}//getOptions

	//! Sets the module option for the graph embedding algorithm.
	void setEmbedder(EmbedderModule *pEmbedder) {
		m_embedder.set(pEmbedder);
	}

	//currently not operable
	//! Sets fixed embedding option, if set to true no embedder is called
	//  void fixEmbed(bool b) {m_fixEmbed = b;}

	//! Sets fixed embedding option, if set to true no embedder is called
	//  void fixEmbed(bool b) {m_fixEmbed = b;}

	void setImprovementCompactor(OrthoCompactionModule* ocm) {
		m_impCompactor.set(ocm);
	}
	//! Sets fixed embedding option, if set to true no embedder is called
	//  void fixEmbed(bool b) {m_fixEmbed = b;}
	void setconstructiveCompactor(OrthoCompactionModule* ocm) {
		m_conCompactor.set(ocm);
	}

protected:
	void classifyEdges(PlanRepUML &PG, adjEntry &adjExternal);

	void doCall(
		PlanRep &PG,
		adjEntry adjExternal,
		GridLayout &gridLayout,
		IPoint &boundingBox,
		bool fixEmbedding=false);

private:
	// compute bounding box and move final drawing such that it is 0 aligned
	// respecting margins
	void computeBoundingBox(const PlanRep &PG, GridLayout &drawing);


	// options
	//inherits a double separation from GridLayoutModule which is used for mapping into GraphAttributes
	int m_margin;
	int m_optionProfile;
	//settings for scaling compaction
	//bool m_useScalingCompaction;
	//int m_scalingSteps;
	//mainly used for OrthoShaper traditional/progressive
	int m_orthoStyle;
	int m_bendBound; //!< bounds number of bends per edge in ortho shaper
	ModuleOption<EmbedderModule> m_embedder;  //!< The planar embedder module.
	ModuleOption<OrthoCompactionModule> m_conCompactor; // The module used for the compaction steps.
	ModuleOption<OrthoCompactionModule> m_impCompactor; // The module used for the compaction steps.
	//bool m_fixEmbed; //!< Embedder is only called if set to false
};


} // end namespace ogdf


#endif
