/*
 * #%L
 * Fiji distribution of ImageJ for the life sciences.
 * %%
 * Copyright (C) 2007 - 2015 Fiji
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

/* Bob Dougherty August 10, 2006

Input: 3D Distance map (32-bit stack)
Output: Distance ridge resulting from a local scan of the distance map.  Overwrites the input.
Note: Non-background points that are not part of the distance ridge are assiged a VERY_SMALL_VALUE.
This is used for subsequent processing by other plugins to find the local thickness.
Reference: T. Holdegrand and P. Ruegsegger, "A new method for the model-independent assessment of
thickness in three-dimensional images," Journal of Microscopy, Vol. 185 Pt. 1, January 1997 pp 67-75.

Version 1: August 10-11, 2006.  Subtracts 0.5 from the distances.
Version 1.01: September 6, 2006.  Corrected some typos in the comments.
Version 1.01: Sept. 7, 2006.  More tiny edits.
Version 2: Sept. 25, 2006.  Creates a separate image stack for symmetry.
                            Temporary version that is very conservative.
                            Admittedly does not produce much impovement on real images.
Version 3: Sept. 30, 2006.  Ball calculations based on grid points.  Should be much more accurate.
Version 3.1 Oct. 1, 2006.  Faster scanning of search points.



 License:
	Copyright (c)  2006, OptiNav, Inc.
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
public class Distance_Ridge implements PlugInFilter {

	private ImagePlus imp;
	private ImagePlus resultImage;

	public float[][] data;
	public int w, h, d;
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
		w = stack.getWidth();
		h = stack.getHeight();
		d = imp.getStackSize();
		// Create 32 bit floating point stack for output, s. Will also use it for g
		// in Transormation 1.
		final ImageStack newStack = new ImageStack(w, h);
		final float[][] sNew = new float[d][];
		for (int k = 0; k < d; k++) {
			final ImageProcessor ipk = new FloatProcessor(w, h);
			newStack.addSlice(null, ipk);
			sNew[k] = (float[]) ipk.getPixels();
		}
		// Create reference to input data
		final float[][] s = new float[d][];
		for (int k = 0; k < d; k++)
			s[k] = (float[]) stack.getPixels(k + 1);
		// Do it
		int k1, j1, i1, dz, dy, dx;
		boolean notRidgePoint;
		float[] sk1;
		float[] sk, skNew;
		int sk0Sq, sk0SqInd, sk1Sq;
		// Find the largest distance in the data
		IJ.showStatus("Distance Ridge: scanning the data");
		float distMax = 0;
		for (int k = 0; k < d; k++) {
			sk = s[k];
			for (int j = 0; j < h; j++) {
				for (int i = 0; i < w; i++) {
					final int ind = i + w * j;
					if (sk[ind] > distMax) distMax = sk[ind];
				}
			}
		}
		final int rSqMax = (int) (distMax * distMax + 0.5f) + 1;
		final boolean[] occurs = new boolean[rSqMax];
		for (int i = 0; i < rSqMax; i++)
			occurs[i] = false;
		for (int k = 0; k < d; k++) {
			sk = s[k];
			for (int j = 0; j < h; j++) {
				for (int i = 0; i < w; i++) {
					final int ind = i + w * j;
					occurs[(int) (sk[ind] * sk[ind] + 0.5f)] = true;
				}
			}
		}
		int numRadii = 0;
		for (int i = 0; i < rSqMax; i++) {
			if (occurs[i]) numRadii++;
		}
		// Make an index of the distance-squared values
		final int[] distSqIndex = new int[rSqMax];
		final int[] distSqValues = new int[numRadii];
		int indDS = 0;
		for (int i = 0; i < rSqMax; i++) {
			if (occurs[i]) {
				distSqIndex[i] = indDS;
				distSqValues[indDS++] = i;
			}
		}
		// Build template
		// The first index of the template is the number of nonzero components
		// in the offest from the test point to the remote point. The second
		// index is the radii index (of the test point). The value of the template
		// is the minimum square radius of the remote point required to cover the
		// ball of the test point.
		IJ.showStatus("Distance Ridge: creating search templates");
		final int[][] rSqTemplate = createTemplate(distSqValues);
		int numCompZ, numCompY, numCompX, numComp;
		for (int k = 0; k < d; k++) {
			IJ.showStatus("Distance Ridge: processing slice " + (k + 1) + "/" + (d +
				1));
			// IJ.showProgress(k/(1.*d));
			sk = s[k];
			skNew = sNew[k];
			for (int j = 0; j < h; j++) {
				for (int i = 0; i < w; i++) {
					final int ind = i + w * j;
					if (sk[ind] > 0) {
						notRidgePoint = false;
						sk0Sq = (int) (sk[ind] * sk[ind] + 0.5f);
						sk0SqInd = distSqIndex[sk0Sq];
						for (dz = -1; dz <= 1; dz++) {
							k1 = k + dz;
							if ((k1 >= 0) && (k1 < d)) {
								sk1 = s[k1];
								if (dz == 0) {
									numCompZ = 0;
								}
								else {
									numCompZ = 1;
								}
								for (dy = -1; dy <= 1; dy++) {
									j1 = j + dy;
									if ((j1 >= 0) && (j1 < h)) {
										if (dy == 0) {
											numCompY = 0;
										}
										else {
											numCompY = 1;
										}
										for (dx = -1; dx <= 1; dx++) {
											i1 = i + dx;
											if ((i1 >= 0) && (i1 < w)) {
												if (dx == 0) {
													numCompX = 0;
												}
												else {
													numCompX = 1;
												}
												numComp = numCompX + numCompY + numCompZ;
												if (numComp > 0) {
													sk1Sq = (int) (sk1[i1 + w * j1] * sk1[i1 + w * j1] +
														0.5f);
													if (sk1Sq >= rSqTemplate[numComp - 1][sk0SqInd])
														notRidgePoint = true;
												}
											} // if in grid for i1
											if (notRidgePoint) break;
										} // dx
									} // if in grid for j1
									if (notRidgePoint) break;
								} // dy
							} // if in grid for k1
							if (notRidgePoint) break;
						} // dz
						if (!notRidgePoint) skNew[ind] = sk[ind];
					} // if not in background
				} // i
			} // j
		} // k
		IJ.showStatus("Distance Ridge complete");

		final String title = stripExtension(imp.getTitle());
		final int slices = imp.getNSlices();
		final int channels = imp.getNChannels();
		final int frames = imp.getNFrames();

		resultImage = IJ.createHyperStack(title  + "_DR", w, h, channels, slices, frames, 32);
		resultImage.setStack(newStack, channels, slices, frames);
		resultImage.getProcessor().setMinAndMax(0, distMax);

		if (!runSilent) {
			resultImage.show();
			IJ.run("Fire");
		}
	}

	// For each offset from the origin, (dx,dy,dz), and each radius-squared,
	// rSq, find the smallest radius-squared, r1Squared, such that a ball
	// of radius r1 centered at (dx,dy,dz) includes a ball of radius
	// rSq centered at the origin. These balls refer to a 3D integer grid.
	// The set of (dx,dy,dz) points considered is a cube center at the origin.
	// The size of the comptued array could be considerably reduced by symmetry,
	// but then the time for the calculation using this array would increase
	// (and more code would be needed).
	int[][] createTemplate(final int[] distSqValues) {
		final int[][] t = new int[3][];
		t[0] = scanCube(1, 0, 0, distSqValues);
		t[1] = scanCube(1, 1, 0, distSqValues);
		t[2] = scanCube(1, 1, 1, distSqValues);
		return t;
	}

	// For a list of r^2 values, find the smallest r1^2 values such
	// that a "ball" of radius r1 centered at (dx,dy,dz) includes a "ball"
	// of radius r centered at the origin. "Ball" refers to a 3D integer grid.
	int[] scanCube(final int dx, final int dy, final int dz,
		final int[] distSqValues)
	{
		final int numRadii = distSqValues.length;
		final int[] r1Sq = new int[numRadii];
		if ((dx == 0) && (dy == 0) && (dz == 0)) {
			for (int rSq = 0; rSq < numRadii; rSq++) {
				r1Sq[rSq] = Integer.MAX_VALUE;
			}
		}
		else {
			final int dxAbs = -Math.abs(dx);
			final int dyAbs = -Math.abs(dy);
			final int dzAbs = -Math.abs(dz);
			for (int rSqInd = 0; rSqInd < numRadii; rSqInd++) {
				final int rSq = distSqValues[rSqInd];
				int max = 0;
				final int r = 1 + (int) Math.sqrt(rSq);
				int scank, scankj;
				int dk, dkji;
				final int iBall;
				int iPlus;
				for (int k = 0; k <= r; k++) {
					scank = k * k;
					dk = (k - dzAbs) * (k - dzAbs);
					for (int j = 0; j <= r; j++) {
						scankj = scank + j * j;
						if (scankj <= rSq) {
							iPlus = ((int) Math.sqrt(rSq - scankj)) - dxAbs;
							dkji = dk + (j - dyAbs) * (j - dyAbs) + iPlus * iPlus;
							if (dkji > max) max = dkji;
						}
					}
				}
				r1Sq[rSqInd] = max;
			}
		}
		return r1Sq;
	} // Modified from ImageJ code by Wayne Rasband

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
}
