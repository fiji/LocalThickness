package sc.fiji.localThickness;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.junit.Before;
import org.junit.Test;
import sc.fiji.localThickness.MaskThicknessMapWithOriginal;

import java.awt.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

/**
 * @author Richard Domander
 */
public class MaskThicknessMapWithOriginalTest
{
    private static final int WHITE = 0xFF;
    private static final int BLACK = 0x00;
    private static final int WIDTH = 10;
    private static final int HALF_WIDTH = WIDTH / 2;
    private static final int HEIGHT = 10;
    private static final int DEPTH = 2;

    private MaskThicknessMapWithOriginal thicknessMasker;

    @Before
    public void setUp() {
        thicknessMasker = new MaskThicknessMapWithOriginal();
        thicknessMasker.threshold = 128;
        thicknessMasker.inverse = false;
    }

    @Test (expected = NullPointerException.class)
    public void testTrimOverhangThrowsNullPointerExceptionWhenOriginalImageIsNull() throws Exception {
        ImagePlus thickness = createTestImage("thickness", WIDTH, HEIGHT, 1, 32, BLACK);
        thicknessMasker.trimOverhang(null, thickness);
    }

    @Test (expected = NullPointerException.class)
    public void testTrimOverhangThrowsNullPointerExceptionWhenThicknessImageIsNull() throws Exception {
        ImagePlus original = createTestImage("thickness", WIDTH, HEIGHT, 1, 8, BLACK);
        thicknessMasker.trimOverhang(original, null);
    }

    @Test (expected = IllegalArgumentException.class)
    public void testTrimOverhangThrowsIllegalArgumentExceptionWhenOriginalImageHasWrongBitDepth() throws Exception {
        ImagePlus floatImage = createTestImage("Float image", WIDTH, HEIGHT, 1, 32, BLACK);
        ImagePlus floatImage2 = createTestImage("Another float image", WIDTH, HEIGHT, 1, 32, BLACK);
        thicknessMasker.trimOverhang(floatImage, floatImage2);
    }

    @Test (expected = IllegalArgumentException.class)
    public void testTrimOverhangThrowsIllegalArgumentExceptionWhenThicknessImageHasWrongBitDepth() throws Exception {
        ImagePlus eightBitImage = createTestImage("8-bit image", WIDTH, HEIGHT, 1, 8, BLACK);
        ImagePlus eightBitImage2 = createTestImage("Another 8-bit image", WIDTH, HEIGHT, 1, 8, BLACK);
        thicknessMasker.trimOverhang(eightBitImage, eightBitImage2);
    }

    @Test (expected = IllegalArgumentException.class)
    public void testTrimOverhangThrowsIllegalArgumentExceptionWhenImageWidthDiffers() throws Exception {
        ImagePlus image = createTestImage("Original", WIDTH, HEIGHT, DEPTH, 8, BLACK);
        ImagePlus tooWide = createTestImage("Too wide", WIDTH + 1, HEIGHT, DEPTH, 32, BLACK);
        thicknessMasker.trimOverhang(image, tooWide);
    }

    @Test (expected = IllegalArgumentException.class)
    public void testTrimOverhangThrowsIllegalArgumentExceptionWhenImageHeightDiffers() throws Exception {
        ImagePlus image = createTestImage("Original", WIDTH, HEIGHT, DEPTH, 8, BLACK);
        ImagePlus tooTall = createTestImage("Too wide", WIDTH, HEIGHT + 1, DEPTH, 32, BLACK);
        thicknessMasker.trimOverhang(image, tooTall);
    }

    @Test (expected = IllegalArgumentException.class)
    public void testTrimOverhangThrowsIllegalArgumentExceptionWhenImageDepthDiffers() throws Exception {
        ImagePlus image = createTestImage("Original", WIDTH, HEIGHT, DEPTH, 8, BLACK);
        ImagePlus tooTall = createTestImage("Too deep", WIDTH, HEIGHT, DEPTH + 1, 32, BLACK);
        thicknessMasker.trimOverhang(image, tooTall);
    }

    @Test
    public void testTrimOverhangMasksPixelsCorrectly() throws Exception {
        // Create an 8-bit test mask that's half black and half white (horizontally),
        // and an all white test thickness image
        ImageStack stack = new ImageStack(WIDTH, HEIGHT);
        Roi left = new Roi(0, 0, HALF_WIDTH, HEIGHT);
        Roi right = new Roi(HALF_WIDTH, 0, HALF_WIDTH, HEIGHT);
        ImageProcessor processor = new ByteProcessor(WIDTH, HEIGHT);
        processor.setColor(WHITE);
        processor.fill(left);
        processor.setColor(BLACK);
        processor.fill(right);
        for (int i = 0; i < DEPTH; i++) {
            stack.addSlice(processor);
        }
        ImagePlus maskImage = new ImagePlus("testMask", stack);

        ImagePlus thicknessImage = createTestImage("testThickness", WIDTH, HEIGHT, DEPTH, 32, WHITE);

        // Test masking
        ImagePlus resultImage = thicknessMasker.trimOverhang(maskImage, thicknessImage);
        assertNotEquals(null, resultImage);
        assertNotEquals(null, thicknessMasker.getResultImage());
        assertEquals("Input mask image must not change", true, areaMatchesColor(maskImage, left.getBounds(), WHITE));
        assertEquals("Input mask image must not change", true, areaMatchesColor(maskImage, right.getBounds(), BLACK));
        assertEquals("Input thickness image must not change", true, areaMatchesColor(thicknessImage,
                new Rectangle(0, 0, WIDTH, HEIGHT), WHITE));
        assertEquals(true, areaMatchesColor(resultImage, left.getBounds(), WHITE));
        assertEquals(true, areaMatchesColor(resultImage, right.getBounds(), BLACK));

        // Test inverted masking
        thicknessMasker.inverse = true;
        resultImage = thicknessMasker.trimOverhang(maskImage, thicknessImage);
        assertEquals(true, areaMatchesColor(resultImage, left.getBounds(), BLACK));
        assertEquals(true, areaMatchesColor(resultImage, right.getBounds(), WHITE));
    }

    /**
     * Checks that the given are in the image matches the color
     *
     * @param image     A 32-bit float image
     * @param bounds    The limits of the area in the image
     * @param color     The expected color
     * @return          True if all the pixels in the area match the color
     */
    private static boolean areaMatchesColor(ImagePlus image, Rectangle bounds, int color) {
        if (image == null || image.getNSlices() < 1) {
            return false;
        }

        final int minX = bounds.x;
        final int maxX = bounds.x + bounds.width;
        final int minY = bounds.y;
        final int maxY = bounds.y + bounds.height;

        if (minX < 0 || minY < 0 || maxX > image.getWidth() || maxY > image.getHeight()) {
            return false;
        }

        ImageStack stack = image.getStack();

        ImageProcessor processor;
        for (int z = 1; z <= image.getNSlices(); z++) {
            processor = stack.getProcessor(z);
            for (int y = minY; y < maxY; y++) {
                for (int x = minX; x < maxX; x++) {
                    int pixelColor = (int)processor.getf(x, y);
                    if (pixelColor != color) {
                        System.out.println("(" + x + "," + y + ") " + pixelColor + " != " + color);
                        return false;
                    }
                }
            }
        }

        return true;
    }

    private static ImagePlus createTestImage(String title, int width, int height, int slices, int bitDepth, int color) {
        ImageStack stack = new ImageStack(width, height);
        ImageProcessor processor;

        switch (bitDepth) {
            case 8:
                processor = new ByteProcessor(width, height);
                break;
            case 32:
                processor = new FloatProcessor(width, height);
                break;
            default:
                return null;
        }

        processor.setColor(color);
        processor.fill();

        for (int i = 0; i < slices; i++) {
            stack.addSlice(processor);
        }

        ImagePlus testImage = new ImagePlus(title, stack);
        return  testImage;
    }
}
