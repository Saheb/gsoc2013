/*
 * $Revision: 2523 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-07-03 00:29:27 +0530 (Tue, 03 Jul 2012) $
 ***************************************************************/

/** \file
 * \brief Declaration of class PlanarityGrid which implements an
 *        energy function where the energy of a layout depends
 *        on the number of crossings.
 *
 * Uses the UniformGris Class to compute the number of crossings.
 *
 * \author Rene Weiskircher
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

#ifndef OGDF_PLANARITYGRID_H
#define OGDF_PLANARITYGRID_H


#include <ogdf/internal/energybased/EnergyFunction.h>
#include <ogdf/internal/energybased/UniformGrid.h>


namespace ogdf {


class PlanarityGrid: public EnergyFunction {
public:
	//initializes data structures to speed up later computations
	PlanarityGrid(GraphAttributes &AG);
	~PlanarityGrid();
	// computes energy of initial layout and stores it in m_energy
	void computeEnergy();
private:
	// computes energy of candidate
	void compCandEnergy();
	// changes internal data if candidate is taken
	void internalCandidateTaken();
#ifdef OGDF_DEBUG
		virtual void printInternalData() const;
#endif
	const GraphAttributes &m_layout; //The current layout
	UniformGrid *m_currentGrid; //stores grid for current layout
	UniformGrid *m_candidateGrid; //stores grid for candidate layout
}; // class Planarity


}// namespace ogdf

#endif
