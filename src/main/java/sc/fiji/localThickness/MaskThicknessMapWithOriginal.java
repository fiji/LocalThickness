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
import ij.process.ImageProcessor;

/**
 * An additional image processing step to LocalThickness, which removes pixels
 * from the thickness map, which are background in the original image. This is
 * to avoid volume dilation in the thickness map, i.e. artifacts in the map may
 * affect statistical measures of the sample.
 *
 * @author Richard Domander
 * @author Michael Doube
 */
public class MaskThicknessMapWithOriginal {

	/**
	 * Pixels with values &lt; threshold are considered background
	 */
	public int threshold = EDT_S1D.DEFAULT_THRESHOLD;

	/**
	 * If true, inverts the threshold condition, i.e. value &gt;= threshold is
	 * background
	 */
	public boolean inverse = EDT_S1D.DEFAULT_INVERSE;

	private ImagePlus resultImage = null;

	/**
	 * Creates a copy of the thicknessMap, where "overhanging" pixels have been
	 * removed. Pixel p{x,y} in the map is considered overhanging if the pixel
	 * q{x=p.x,y=p.y} in the original image is background.
	 *
	 * @param original The original, unprocessed 8-bit binary image
	 * @param thicknessMap The 32-bit thickness map image produced from the
	 *          original
	 * @throws NullPointerException if original == null or thicknessMap == null
	 * @throws IllegalArgumentException if original.getBitDepth() != 8 or
	 *           thicknessMap.getBitDepth() != 32
	 * @throws IllegalArgumentException if the dimensions of the images do not
	 *           match
	 * @return A 32-bit thickness map image where the overhanging pixels have been
	 *         set to 0
	 */
	public ImagePlus trimOverhang(final ImagePlus original,
		final ImagePlus thicknessMap)
	{
		if (original == null || thicknessMap == null) {
			throw new NullPointerException("Images must not be null");
		}

		if (original.getBitDepth() != 8 || thicknessMap.getBitDepth() != 32) {
			throw new IllegalArgumentException(
				"One or both of the images have the wrong bit depth");
		}

		final int w = original.getWidth();
		final int h = original.getHeight();
		final int d = original.getImageStackSize();

		if (w != thicknessMap.getWidth() || h != thicknessMap.getHeight() ||
			d != thicknessMap.getImageStackSize())
		{
			throw new IllegalArgumentException(
				"The dimensions of the images do not match");
		}

		resultImage = thicknessMap.duplicate();
		resultImage.setTitle(thicknessMap.getTitle() + "_MASK");

		final ImageStack originalStack = original.getImageStack();
		final ImageStack resultStack = resultImage.getImageStack();

		ImageProcessor originalProcessor;
		ImageProcessor resultProcessor;
		for (int z = 1; z <= d; z++) {
			IJ.showStatus("Masking thickness map...");
			IJ.showProgress(z, d);
			originalProcessor = originalStack.getProcessor(z);
			resultProcessor = resultStack.getProcessor(z);
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					final int value = originalProcessor.get(x, y);
					if ((value < threshold && !inverse) || (value >= threshold &&
						inverse))
					{
						resultProcessor.set(x, y, 0);
					}
				}
			}
		}

		return getResultImage();
	}

	/**
	 * @return The result of the last call to trimOverhang. Returns null if the
	 *         method hasn't been successfully called.
	 */
	public ImagePlus getResultImage() {
		return resultImage;
	}

	/**
	 * Remove references to instance variables to allow garbage collection
	 */
	public void purge() {
		resultImage = null;
	}
}
