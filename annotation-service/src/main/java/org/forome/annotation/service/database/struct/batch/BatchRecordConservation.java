/*
 *  Copyright (c) 2020. Vladimir Ulitin, Partners Healthcare and members of Forome Association
 *
 *  Developed by Vladimir Ulitin and Michael Bouzinier
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 * 	 http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

package org.forome.annotation.service.database.struct.batch;

import org.forome.annotation.struct.Interval;
import org.forome.annotation.struct.Position;
import org.forome.annotation.utils.bits.ShortBits;

public class BatchRecordConservation {

	private final Interval interval;
	private final byte[] bytes;
	private final int offsetBytes;

	protected BatchRecordConservation(Interval interval, byte[] bytes, int offsetBytes) {
		this.interval = interval;
		this.bytes = bytes;
		this.offsetBytes = offsetBytes;
	}

	public float getGerpN(Position position) {
		int ioffset = offsetBytes + (position.value - interval.start) * 2;
		return (float) ShortBits.fromByteArray(bytes, ioffset) / 1000.0f;
	}

	public float getGerpRS(Position position) {
		int ioffset = offsetBytes + (position.value - interval.start) * 2 + 2;
		return (float) ShortBits.fromByteArray(bytes, ioffset) / 1000.0f;
	}

	protected int getLengthBytes() {
		return (2 + 2) * (interval.end - interval.start + 1);
	}
}
