/*
 * $Revision: 3386 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-10 17:30:02 +0530 (Wed, 10 Apr 2013) $
 ***************************************************************/

/*!\file
 * \author Matthias Elf
 *
 * \par License:
 * This file is part of ABACUS - A Branch And CUt System
 * Copyright (C) 1995 - 2003
 * University of Cologne, Germany
 *
 * \par
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * \par
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * \see http://www.gnu.org/copyleft/gpl.html
 */

#ifndef ABA__POOLSLOTREF_INC
#define ABA__POOLSLOTREF_INC

#include <ogdf/abacus/poolslotref.h>

namespace abacus {


template <class BaseType, class CoType>
ostream &operator<<(ostream &out, const PoolSlotRef<BaseType, CoType> &slot)
{
	if (slot.conVar())
		slot.conVar()->print(out);
	return out;
}


template<class BaseType, class CoType>
void PoolSlotRef<BaseType, CoType>::printDifferentVersionError() const
{
	Logger::ilout(Logger::LL_MINOR) << "PoolSlotRef::con: Version of reference to slot " << version_
	 << " and version of slot " << slot_->version() << " differ." << endl;
}


template<class BaseType, class CoType>
void PoolSlotRef<BaseType, CoType>::slot(PoolSlot<BaseType, CoType> *s)
{
	ConVar *cv;

	if(slot_) {
		cv = slot_->conVar();
		if (cv && version_ == slot_->version())
			cv->removeReference();
	}

	slot_    = s;
	version_ = s->version();
	cv       = slot_->conVar();
	if(cv) cv->addReference();
}

} // namespace abacus

#endif // PoolSlotRef_INC
