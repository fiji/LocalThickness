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
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/* Bob Dougherty 8/8/2006
Saito-Toriwaki algorithm for Euclidian Distance Transformation.
Direct application of Algorithm 1.
Version S1A: lower memory usage.
Version S1A.1 A fixed indexing bug for 666-bin data set
Version S1A.2 Aug. 9, 2006.  Changed noResult value.
Version S1B Aug. 9, 2006.  Faster.
Version S1B.1 Sept. 6, 2006.  Changed comments.
Version S1C Oct. 1, 2006.  Option for inverse case.
                           Fixed inverse behavior in y and z directions.
Version D July 30, 2007.  Multithread processing for step 2.

This version assumes the input stack is already in memory, 8-bit, and
outputs to a new 32-bit stack.  Versions that are more stingy with memory
may be forthcoming.

 License:
	Copyright (c) 2006, OptiNav, Inc.
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
public class EDT_S1D implements PlugInFilter {

	public static final int DEFAULT_THRESHOLD = 128;
	public static final boolean DEFAULT_INVERSE = false;

	private ImagePlus imp;
	private boolean cancelled;

	public byte[][] data;
	public int thresh = DEFAULT_THRESHOLD;
	public boolean inverse = DEFAULT_INVERSE;
	public boolean showOptions = true;
	public boolean runSilent = false;
	private ImagePlus resultImage;

	@Override
	public int setup(final String arg, final ImagePlus imp) {
		this.imp = imp;
		return DOES_8G;
	}

	@Override
	public void run(final ImageProcessor ip) {
		resultImage = null;

		final ImageStack stack = imp.getStack();
		final int w = stack.getWidth();
		final int h = stack.getHeight();
		final int d = imp.getStackSize();
		final int nThreads = Runtime.getRuntime().availableProcessors();

		cancelled = false;
		if (showOptions) {
			if (!getScale()) {
				cancelled = true;
				return;
			}
		}

		// Create references to input data
		data = new byte[d][];
		for (int k = 0; k < d; k++)
			data[k] = (byte[]) stack.getPixels(k + 1);
		// Create 32 bit floating point stack for output, s. Will also use it for g
		// in Transormation 1.
		final ImageStack sStack = new ImageStack(w, h);
		final float[][] s = new float[d][];
		for (int k = 0; k < d; k++) {
			final ImageProcessor ipk = new FloatProcessor(w, h);
			sStack.addSlice(null, ipk);
			s[k] = (float[]) ipk.getPixels();
		}
		float[] sk;
		// Transformation 1. Use s to store g.
		IJ.showStatus("EDT transformation 1/3");
		final Step1Thread[] s1t = new Step1Thread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			s1t[thread] = new Step1Thread(thread, nThreads, w, h, d, thresh, s, data);
			s1t[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				s1t[thread].join();
			}
		}
		catch (final InterruptedException ie) {
			IJ.error("A thread was interrupted in step 1 .");
		}
		// Transformation 2. g (in s) -> h (in s)
		IJ.showStatus("EDT transformation 2/3");
		final Step2Thread[] s2t = new Step2Thread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			s2t[thread] = new Step2Thread(thread, nThreads, w, h, d, s);
			s2t[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				s2t[thread].join();
			}
		}
		catch (final InterruptedException ie) {
			IJ.error("A thread was interrupted in step 2 .");
		}
		// Transformation 3. h (in s) -> s
		IJ.showStatus("EDT transformation 3/3");
		final Step3Thread[] s3t = new Step3Thread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			s3t[thread] = new Step3Thread(thread, nThreads, w, h, d, s, data);
			s3t[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				s3t[thread].join();
			}
		}
		catch (final InterruptedException ie) {
			IJ.error("A thread was interrupted in step 3 .");
		}
		// Find the largest distance for scaling
		// Also fill in the background values.
		float distMax = 0;
		final int wh = w * h;
		float dist;
		for (int k = 0; k < d; k++) {
			sk = s[k];
			for (int ind = 0; ind < wh; ind++) {
				if (((data[k][ind] & 255) < thresh) ^ inverse) {
					sk[ind] = 0;
				}
				else {
					dist = (float) Math.sqrt(sk[ind]);
					sk[ind] = dist;
					distMax = (dist > distMax) ? dist : distMax;
				}
			}
		}

		IJ.showProgress(1.0);
		IJ.showStatus("Done");

		final String title = stripExtension(imp.getTitle());
		final int slices = imp.getNSlices();
		final int channels = imp.getNChannels();
		final int frames = imp.getNFrames();

		resultImage = IJ.createHyperStack(title + "_EDT", w, h, channels, slices, frames, 32);
		resultImage.setStack(sStack, channels, slices, frames);
		resultImage.getProcessor().setMinAndMax(0, distMax);

		if (!runSilent) {
			resultImage.show();
			IJ.run("Fire");
		}
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

	boolean getScale() {
		thresh = (int) Prefs.get("edtS1.thresh", 128);
		inverse = Prefs.get("edtS1.inverse", false);
		final GenericDialog gd = new GenericDialog("EDT...", IJ.getInstance());
		gd.addNumericField("Threshold (1 to 255; value < thresh is background)",
			thresh, 0);
		gd.addCheckbox("Inverse case (background when value >= thresh)", inverse);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		thresh = (int) gd.getNextNumber();
		inverse = gd.getNextBoolean();
		Prefs.set("edtS1.thresh", thresh);
		Prefs.set("edtS1.inverse", inverse);
		return true;
	}

	public boolean gotCancelled() {
		return cancelled;
	}

	class Step1Thread extends Thread {

		int thread, nThreads, w, h, d, thresh;
		float[][] s;
		byte[][] data;

		public Step1Thread(final int thread, final int nThreads, final int w,
			final int h, final int d, final int thresh, final float[][] s,
			final byte[][] data)
		{
			this.thread = thread;
			this.nThreads = nThreads;
			this.w = w;
			this.h = h;
			this.d = d;
			this.thresh = thresh;
			this.data = data;
			this.s = s;
		}

		@Override
		public void run() {
			float[] sk;
			byte[] dk;
			int n = w;
			if (h > n) n = h;
			if (d > n) n = d;
			final int noResult = 3 * (n + 1) * (n + 1);
			final boolean[] background = new boolean[n];
			final boolean nonempty;
			int test, min;
			for (int k = thread; k < d; k += nThreads) {
				IJ.showProgress(k / (1. * d));
				sk = s[k];
				dk = data[k];
				for (int j = 0; j < h; j++) {
					final int wj = w * j;
					for (int i = 0; i < w; i++) {
						background[i] = ((dk[i + wj] & 255) < thresh) ^ inverse;
					}
					for (int i = 0; i < w; i++) {
						min = noResult;
						for (int x = i; x < w; x++) {
							if (background[x]) {
								test = i - x;
								test *= test;
								min = test;
								break;
							}
						}
						for (int x = i - 1; x >= 0; x--) {
							if (background[x]) {
								test = i - x;
								test *= test;
								if (test < min) min = test;
								break;
							}
						}
						sk[i + wj] = min;
					}
				}
			}
		}// run
	}// Step1Thread

	class Step2Thread extends Thread {

		int thread, nThreads, w, h, d;
		float[][] s;

		public Step2Thread(final int thread, final int nThreads, final int w,
			final int h, final int d, final float[][] s)
		{
			this.thread = thread;
			this.nThreads = nThreads;
			this.w = w;
			this.h = h;
			this.d = d;
			this.s = s;
		}

		@Override
		public void run() {
			float[] sk;
			int n = w;
			if (h > n) n = h;
			if (d > n) n = d;
			final int noResult = 3 * (n + 1) * (n + 1);
			final int[] tempInt = new int[n];
			final int[] tempS = new int[n];
			boolean nonempty;
			int test, min, delta;
			for (int k = thread; k < d; k += nThreads) {
				IJ.showProgress(k / (1. * d));
				sk = s[k];
				for (int i = 0; i < w; i++) {
					nonempty = false;
					for (int j = 0; j < h; j++) {
						tempS[j] = (int) sk[i + w * j];
						if (tempS[j] > 0) nonempty = true;
					}
					if (nonempty) {
						for (int j = 0; j < h; j++) {
							min = noResult;
							delta = j;
							for (int y = 0; y < h; y++) {
								test = tempS[y] + delta * delta--;
								if (test < min) min = test;
							}
							tempInt[j] = min;
						}
						for (int j = 0; j < h; j++) {
							sk[i + w * j] = tempInt[j];
						}
					}
				}
			}
		}// run
	}// Step2Thread

	class Step3Thread extends Thread {

		int thread, nThreads, w, h, d;
		float[][] s;
		byte[][] data;

		public Step3Thread(final int thread, final int nThreads, final int w,
			final int h, final int d, final float[][] s, final byte[][] data)
		{
			this.thread = thread;
			this.nThreads = nThreads;
			this.w = w;
			this.h = h;
			this.d = d;
			this.s = s;
			this.data = data;
		}

		@Override
		public void run() {
			int zStart, zStop, zBegin, zEnd;
			final float[] sk;
			int n = w;
			if (h > n) n = h;
			if (d > n) n = d;
			final int noResult = 3 * (n + 1) * (n + 1);
			final int[] tempInt = new int[n];
			final int[] tempS = new int[n];
			boolean nonempty;
			int test, min, delta;
			for (int j = thread; j < h; j += nThreads) {
				IJ.showProgress(j / (1. * h));
				final int wj = w * j;
				for (int i = 0; i < w; i++) {
					nonempty = false;
					for (int k = 0; k < d; k++) {
						tempS[k] = (int) s[k][i + wj];
						if (tempS[k] > 0) nonempty = true;
					}
					if (nonempty) {
						zStart = 0;
						while ((zStart < (d - 1)) && (tempS[zStart] == 0))
							zStart++;
						if (zStart > 0) zStart--;
						zStop = d - 1;
						while ((zStop > 0) && (tempS[zStop] == 0))
							zStop--;
						if (zStop < (d - 1)) zStop++;

						for (int k = 0; k < d; k++) {
							// Limit to the non-background to save time,
							if (((data[k][i + wj] & 255) >= thresh) ^ inverse) {
								min = noResult;
								zBegin = zStart;
								zEnd = zStop;
								if (zBegin > k) zBegin = k;
								if (zEnd < k) zEnd = k;
								delta = k - zBegin;
								for (int z = zBegin; z <= zEnd; z++) {
									test = tempS[z] + delta * delta--;
									if (test < min) min = test;
									// min = (test < min) ? test : min;
								}
								tempInt[k] = min;
							}
						}
						for (int k = 0; k < d; k++) {
							s[k][i + wj] = tempInt[k];
						}
					}
				}
			}
		}// run
	}// Step2Thread
}
