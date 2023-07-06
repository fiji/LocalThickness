/*
 * #%L
 * Fiji distribution of ImageJ for the life sciences.
 * %%
 * Copyright (C) 2006 - 2020 Fiji developers.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package sc.fiji.localThickness;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/* Bob Dougherty August 1, 2007

Input: 3D Local Thickness map (32-bit stack)
Output: Same as input with border voxels corrected for "jaggies." Non-background voxels
adjacent to background voxels are have their local thickness values replaced by the average of
their non-background neighbors that do not border background points.

August 10.  Version 3 This version also multiplies the local thickness by 2 to conform with the
official definition of local thickness.

 License:
	Copyright (c)  2007, OptiNav, Inc.
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:

		Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.
		Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.
		Neither the name of OptiNav, Inc. nor the names of its contributors
	may be used to endorse or promote products derived from this software
	without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
	EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
public class Clean_Up_Local_Thickness implements PlugInFilter {

	private ImagePlus imp;
	private ImagePlus resultImage;

	public float[][] s, sNew;
	public boolean runSilent = false;

	@Override
	public int setup(final String arg, final ImagePlus imp) {
		this.imp = imp;
		return DOES_32;
	}

	@Override
	public void run(final ImageProcessor ip) {
		resultImage = null;

		final ImageStack stack = imp.getStack();
		final int w = stack.getWidth();
		final int h = stack.getHeight();
		final int d = imp.getStackSize();
		
		// Create 32 bit floating point stack for output, sNew.
		final ImageStack newStack = new ImageStack(w, h);
		sNew = new float[d][];
		for (int k = 0; k < d; k++) {
			final ImageProcessor ipk = new FloatProcessor(w, h);
			newStack.addSlice(null, ipk);
			sNew[k] = (float[]) ipk.getPixels();
		}
		// Create reference to input data
		s = new float[d][];
		for (int k = 0; k < d; k++)
			s[k] = (float[]) stack.getPixels(k + 1);
		// First set the output array to flags:
		// 0 for a background point
		// -1 for a non-background point that borders a background point
		// s (input data) for an interior non-background point
		for (int k = 0; k < d; k++) {
			for (int j = 0; j < h; j++) {
				final int wj = w * j;
				for (int i = 0; i < w; i++) {
					sNew[k][i + wj] = setFlag(i, j, k, w, h, d);
				} // i
			} // j
		} // k
			// Process the surface points. Initially set results to negative values
			// to be able to avoid including them in averages of for subsequent
			// points.
			// During the calculation, positve values in sNew are interior
			// non-background
			// local thicknesses. Negative values are surface points. In this case the
			// value might be -1 (not processed yet) or -result, where result is the
			// average of the neighboring interior points. Negative values are
			// excluded from
			// the averaging.
		for (int k = 0; k < d; k++) {
			for (int j = 0; j < h; j++) {
				final int wj = w * j;
				for (int i = 0; i < w; i++) {
					final int ind = i + wj;
					if (sNew[k][ind] == -1) {
						sNew[k][ind] = -averageInteriorNeighbors(i, j, k, w, h, d);
					}
				} // i
			} // j
		} // k
			// Fix the negative values and double the results
		for (int k = 0; k < d; k++) {
			for (int j = 0; j < h; j++) {
				final int wj = w * j;
				for (int i = 0; i < w; i++) {
					final int ind = i + wj;
					sNew[k][ind] = Math.abs(sNew[k][ind]);
				} // i
			} // j
		} // k
		IJ.showStatus("Clean Up Local Thickness complete");

		final String title = stripExtension(imp.getTitle());
		final int slices = imp.getNSlices();
		final int channels = imp.getNChannels();
		final int frames = imp.getNFrames();

		resultImage = IJ.createHyperStack(title + "_CL", w, h, channels, slices, frames, 32);
		resultImage.setStack(newStack, channels, slices, frames);
		resultImage.getProcessor().setMinAndMax(0, 2 * imp.getProcessor().getMax());

		if (!runSilent) {
			resultImage.show();
			IJ.run("Fire");
		}
	}

	float setFlag(final int i, final int j, final int k, final int w, final int h, final int d) {
		if (s[k][i + w * j] == 0) return 0;
		// change 1
		if (look(i, j, k - 1, w, h, d) == 0) return -1;
		if (look(i, j, k + 1, w, h, d) == 0) return -1;
		if (look(i, j - 1, k, w, h, d) == 0) return -1;
		if (look(i, j + 1, k, w, h, d) == 0) return -1;
		if (look(i - 1, j, k, w, h, d) == 0) return -1;
		if (look(i + 1, j, k, w, h, d) == 0) return -1;
		// change 1 before plus
		if (look(i, j + 1, k - 1, w, h, d) == 0) return -1;
		if (look(i, j + 1, k + 1, w, h, d) == 0) return -1;
		if (look(i + 1, j - 1, k, w, h, d) == 0) return -1;
		if (look(i + 1, j + 1, k, w, h, d) == 0) return -1;
		if (look(i - 1, j, k + 1, w, h, d) == 0) return -1;
		if (look(i + 1, j, k + 1, w, h, d) == 0) return -1;
		// change 1 before minus
		if (look(i, j - 1, k - 1, w, h, d) == 0) return -1;
		if (look(i, j - 1, k + 1, w, h, d) == 0) return -1;
		if (look(i - 1, j - 1, k, w, h, d) == 0) return -1;
		if (look(i - 1, j + 1, k, w, h, d) == 0) return -1;
		if (look(i - 1, j, k - 1, w, h, d) == 0) return -1;
		if (look(i + 1, j, k - 1, w, h, d) == 0) return -1;
		// change 3, k+1
		if (look(i + 1, j + 1, k + 1, w, h, d) == 0) return -1;
		if (look(i + 1, j - 1, k + 1, w, h, d) == 0) return -1;
		if (look(i - 1, j + 1, k + 1, w, h, d) == 0) return -1;
		if (look(i - 1, j - 1, k + 1, w, h, d) == 0) return -1;
		// change 3, k-1
		if (look(i + 1, j + 1, k - 1, w, h, d) == 0) return -1;
		if (look(i + 1, j - 1, k - 1, w, h, d) == 0) return -1;
		if (look(i - 1, j + 1, k - 1, w, h, d) == 0) return -1;
		if (look(i - 1, j - 1, k - 1, w, h, d) == 0) return -1;
		return s[k][i + w * j];
	}

	float averageInteriorNeighbors(final int i, final int j, final int k, final int w, final int h, final int d) {
		int n = 0;
		float sum = 0;
		// change 1
		float value = lookNew(i, j, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i, j, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i, j - 1, k, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i, j + 1, k, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j, k, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i + 1, j, k, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		// change 1 before plus
		value = lookNew(i, j + 1, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i, j + 1, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i + 1, j - 1, k, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i + 1, j + 1, k, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i + 1, j, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		// change 1 before minus
		value = lookNew(i, j - 1, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i, j - 1, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j - 1, k, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j + 1, k, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i + 1, j, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		// change 3, k+1
		value = lookNew(i + 1, j + 1, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i + 1, j - 1, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j + 1, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j - 1, k + 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		// change 3, k-1
		value = lookNew(i + 1, j + 1, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i + 1, j - 1, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j + 1, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		value = lookNew(i - 1, j - 1, k - 1, w, h, d);
		if (value > 0) {
			n++;
			sum += value;
		}
		if (n > 0) return sum / n;
		return s[k][i + w * j];
	}

	float look(final int i, final int j, final int k, final int w, final int h, final int d) {
		if ((i < 0) || (i >= w)) return -1;
		if ((j < 0) || (j >= h)) return -1;
		if ((k < 0) || (k >= d)) return -1;
		return s[k][i + w * j];
	}

	// A positive result means this is an interior, non-background, point.
	float lookNew(final int i, final int j, final int k, final int w, final int h, final int d) {
		if ((i < 0) || (i >= w)) return -1;
		if ((j < 0) || (j >= h)) return -1;
		if ((k < 0) || (k >= d)) return -1;
		return sNew[k][i + w * j];
	}

	// Modified from ImageJ code by Wayne Rasband
	String stripExtension(String name) {
		if (name != null) {
			final int dotIndex = name.lastIndexOf(".");
			if (dotIndex >= 0) name = name.substring(0, dotIndex);
		}
		return name;
	}

	public ImagePlus getResultImage() {
		return resultImage;
	}

	/**
	 * Remove references to instance variables to allow garbage collection
	 */
	public void purge() {
		s = null;
		sNew = null;
		resultImage = null;
		imp = null;
	}
}
