/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 23:37:31 +0530 (Thu, 04 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class FixedEmbeddingInserterUML.
 *
 * \author Carsten Gutwenger
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

#ifndef OGDF_FIX_EMB_INSERTER_UML_H
#define OGDF_FIX_EMB_INSERTER_UML_H


#include <ogdf/module/UMLEdgeInsertionModule.h>
#include <ogdf/planarity/RemoveReinsertType.h>


namespace ogdf {


	//! Edge insertion module that inserts each edge optimally into a fixed embedding.
	class OGDF_EXPORT FixedEmbeddingInserterUML : public UMLEdgeInsertionModule
	{
	public:
		//! Creates an instance of variable embedding edge inserter with default settings.
		FixedEmbeddingInserterUML();

		//! Creates an instance of fixed embedding edge inserter with the same settings as \a inserter.
		FixedEmbeddingInserterUML(const FixedEmbeddingInserterUML &inserter);

		//! Destructor.
		~FixedEmbeddingInserterUML() { }

		//! Returns a new instance of the fixed embedding inserter with the same option settings.
		virtual UMLEdgeInsertionModule *clone() const;

		//! Assignment operator. Copies option settings only.
		FixedEmbeddingInserterUML &operator=(const FixedEmbeddingInserterUML &inserter);


		/**
		 *  @name Optional parameters
		 *  @{
		 */

		//! Sets the remove-reinsert postprocessing method.
		void removeReinsert(RemoveReinsertType rrOption) {
			m_rrOption = rrOption;
		}

		//! Returns the current setting of the remove-reinsert postprocessing method.
		RemoveReinsertType removeReinsert() const {
			return m_rrOption;
		}


		//! Sets the option <i>percentMostCrossed</i> to \a percent.
		/**
		 * This option determines the portion of most crossed edges used if the remove-reinsert
		 * method is set to #rrMostCrossed. This portion is number of edges * percentMostCrossed() / 100.
		 */
		void percentMostCrossed(double percent) {
			m_percentMostCrossed = percent;
		}

		//! Returns the current setting of option percentMostCrossed.
		double percentMostCrossed() const {
			return m_percentMostCrossed;
		}

		//! Sets the option <i>keepEmbedding</i> to \a keep.
		/**
		 * This option determines if the planar embedding of the planarized representation \a PG passed to the call-method
		 * is preserved, or if always a new embedding is computed. If <i>keepEmbedding</i> is set to true,
		 * \a PG must always be planarly embedded.
		 */
		void keepEmbedding(bool keep) {
			m_keepEmbedding = keep;
		}

		//! Returns the current setting of option <i>keepEmbedding</i>.
		bool keepEmbeding() const {
			return m_keepEmbedding;
		}

		//! @}

	private:
		//! Implements the algorithm call.
		ReturnType doCall(
			PlanRepLight              &pr,
			const Array<edge>         &origEdges,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<__uint32> *pEdgeSubgraph);

		RemoveReinsertType m_rrOption; //!< The remove-reinsert method.
		double m_percentMostCrossed;   //!< The portion of most crossed edges considered.
		bool m_keepEmbedding;
	};

} // end namespace ogdf

#endif
