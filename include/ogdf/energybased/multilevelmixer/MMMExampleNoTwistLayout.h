/*
 * $Revision: 3432 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-22 15:50:23 +0530 (Mon, 22 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief useable example of the Modular Multilevel Mixer
 *
 * \author Gereon Bartel
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

#ifndef OGDF_MMM_EXAMPLE_NO_TWIST_LAYOUT_H
#define OGDF_MMM_EXAMPLE_NO_TWIST_LAYOUT_H

#include <ogdf/module/LayoutModule.h>
#include <ogdf/internal/energybased/MultilevelGraph.h>

namespace ogdf {

/** \brief An example Layout using the Modular Mutlievel Mixer
 *
 * This example is tuned to reduce twists in the final drawing. Use this layout
 * or a variation of it if many twists occur. LocalBiconnectedMerger and
 * BarycenterPlacer are used as merging and placement strategies. The
 * FastMultipoleEmbedder is for force calculation.
 *
 * For an easy variation of the Modular Multilevel Mixer copy the code in call.
 */
class OGDF_EXPORT MMMExampleNoTwistLayout : public LayoutModule
{
public:

	//! Constructor
	MMMExampleNoTwistLayout();

	//! calculates a drawing for the Graph GA
	void call(GraphAttributes &GA);

	void call(GraphAttributes &GA, GraphConstraints & GC) { call(GA); }

	//! calculates a drawing for the Graph MLG
	void call(MultilevelGraph &MLG);

private:

};

} // namespace ogdf

#endif

