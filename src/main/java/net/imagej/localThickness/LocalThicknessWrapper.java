package net.imagej.localThickness;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

/**
 * A class which wraps the different stages of LocalThickness processing.
 *
 * By default LocalThicknessWrapper runs the stages "silently" so that the user doesn't see the intermediate images
 * between the them. Even if the intermediate images are shown, the user won't be able to interrupt the process like in
 * Local_Thickness_Driver (think what happens if an image is closed before WindowManager.getCurrentImage()
 * gets called...).
 *
 * By default the class won't show the options (threshold, inverse) dialog for LocalThickness.
 * Instead it's run with default values defined in EDT_S1D.DEFAULT_INVERSE and EDT_S1D.DEFAULT_THRESHOLD respectively.
 *
 * Mask
 *
 * Calibration
 *
 * @author <a href="mailto:rdomander@rvc.ac.uk">Richard Domander</a>
 */
public class LocalThicknessWrapper implements PlugIn
{
    private static final String DEFAULT_TITLE_SUFFIX = "_LocThk";
    private static final String DEFAULT_TITLE = "ThicknessMap";

    private static final EDT_S1D geometryToDistancePlugin = new EDT_S1D();
    private static final Distance_Ridge distanceRidgePlugin = new Distance_Ridge();
    private static final Local_Thickness_Parallel localThicknessPlugin = new Local_Thickness_Parallel();
    private static final Clean_Up_Local_Thickness thicknessCleaningPlugin = new Clean_Up_Local_Thickness();

    private ImagePlus resultImage = null;
    private boolean showOptions = false;
    private String titleSuffix = DEFAULT_TITLE_SUFFIX;

    // Fields used to set LocalThickness options programmatically
    /**
     * A pixel is considered to be a part of the background if its color < threshold
     */
    public int threshold = EDT_S1D.DEFAULT_THRESHOLD;

    /**
     * Inverts thresholding so that pixels with values >= threshold are considered background.
     */
    public boolean inverse = EDT_S1D.DEFAULT_INVERSE;

    public LocalThicknessWrapper() {
        setSilence(true);
        geometryToDistancePlugin.showOptions = showOptions;
    }

    /**
     * Controls whether the options dialog in EDT_S1D is shown to the user before processing starts.
     *
     * @param show  If true, then the dialog is shown.
     */
    public void setShowOptions(boolean show) {
        showOptions = show;
        geometryToDistancePlugin.showOptions = show;
    }

    /**
     * Sets whether intermediate images are shown to the user, or only the final result
     *
     * @param silent If true, then no intermediate images are shown in the GUI.
     */
    public void setSilence(boolean silent) {
        distanceRidgePlugin.runSilent = silent;
        localThicknessPlugin.runSilent = silent;
        thicknessCleaningPlugin.runSilent = silent;
        geometryToDistancePlugin.runSilent = silent;
    }

    public ImagePlus processImage(ImagePlus image) {
        resultImage = null;
        String originalTitle = Local_Thickness_Driver.stripExtension(image.getTitle());
        if (originalTitle == null || originalTitle.isEmpty()) {
            originalTitle = DEFAULT_TITLE;
        }

        geometryToDistancePlugin.setup("", image);
        if (!showOptions) {
            // set options programmatically
            geometryToDistancePlugin.inverse = inverse;
            geometryToDistancePlugin.thresh = threshold;

            geometryToDistancePlugin.run(null);
        } else {
            // get options from the dialog in EDT_S1D
            geometryToDistancePlugin.run(null);

            if (geometryToDistancePlugin.gotCancelled()) {
                resultImage = null;
                return getResultImage();
            }

            inverse = geometryToDistancePlugin.inverse;
            threshold = geometryToDistancePlugin.thresh;
        }
        resultImage = geometryToDistancePlugin.getResultImage();

        distanceRidgePlugin.setup("", resultImage);
        distanceRidgePlugin.run(null);
        resultImage = distanceRidgePlugin.getResultImage();

        localThicknessPlugin.setup("", resultImage);
        localThicknessPlugin.run(null);
        resultImage = localThicknessPlugin.getResultImage();

        thicknessCleaningPlugin.setup("", resultImage);
        thicknessCleaningPlugin.run(null);
        resultImage = thicknessCleaningPlugin.getResultImage();

        resultImage.setTitle(originalTitle + titleSuffix);
        resultImage.copyScale(image);

        return getResultImage();
    }

    public ImagePlus getResultImage() {
        return resultImage;
    }

    @Override
    public void run(String s) {
        ImagePlus image;

        try {
            image = IJ.getImage();
        } catch(RuntimeException rte) {
            return; // no image currently open
        }

        processImage(image);

        if (resultImage == null) {
            return;
        }

        resultImage.show();
        IJ.run("Fire"); // changes the color palette of the output image
    }

    /**
     * @param imageTitleSuffix  The suffix that's added to the end of the title of the resulting thickness map image
     */
    public void setTitleSuffix(String imageTitleSuffix) {
        titleSuffix = imageTitleSuffix;

        if (titleSuffix == null || titleSuffix.isEmpty()) {
            titleSuffix = DEFAULT_TITLE_SUFFIX;
        }

        assert titleSuffix != null && !this.titleSuffix.isEmpty();
    }
}