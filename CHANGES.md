## V3.24.06.2
   - CTF addeed to projection protocols
## V3.24.06
   - New protocols:
      - apply_segmentation: Applies a segmentation to a set of tomograms.
      - segment_morphology: Segments a tomogram by means of thresholding, different kind of filters and morphological operations.
   - Protocols updated
      - deep_misalignment: updated model, minor improvements
      - extract_subtomos: improve extract subtomos readability
   - Protocols fixed
      - dose_filter: Dose filter, removing extract stacks
   - More xmippTomo
      - New xmippBase repository to manage shared protocols with Xmipp
      - protocols deprecated: score_transform, denoising_confidence
      - Removed tilt particles object (not used)

## V3.23.11
   - Protocols updated
      - applyAlignmentTS
      - extract_particlesstacks: Added angles from Tilt Series
      - extract_subtomos: Keeping directions in extract subtomos
      - project_subtomograms: Added angles from Tilt Series
      - deep_misalignment_detection: Update deep misalignment detection with new version of subtomo extraction and fiducial size options
      - subtraction_subtomo: Improved implementation with MPI and md convert
      - cltomo:
      - resolution_local_monotomo: Allowing odd-even associated to the full tomogram
   - Protocls fixed
      - extract_particlestacks: Dose fixes
      - deep_misalignment_detection:  Recover deleted methods
      - dose_filter: Dose validation fixed
      - score_coordinates: 
   - More xmippTomo
      - Fixed calculateRotationAngleAndShiftsFromTM
      - Test refactoring

## V3.23.07
   - New protocols
      - project_subtomograms (for obtaining sobtomogram projections)
      - extract_particlestacks (extract from tilt series and SetOfTiltSeriesParticle)
      - denoising_confidence
      - extract_particlestacks
      - deep_misalignment_detection (misalignment detection in tomographic reconstructions from high-contrast regions)
      - peak_high_contrast (for detecting high contrast regions in tomographic reconstruction)
   - Protocols updated 
      - project_top: Subtomo projector compatible with hdf stacks. Test fixed. Projections keep orientation, Ignore alignments in projection
      - subtomo_map_back, score_coordinates: Mapback scales any reference to match the tomogram sampling
      - subtomo_map_back: Mapback waits now when scheduled and not the reference is needed,  works with 3d classes, works with "other tomograms"
      - project_subtomograms: Parallelized protocol 
- Update requirements.txt
- Fix tests: protocol_extract_subtomos, monotomo, crop, resize, xmipptomo


## V3.23.03.6:
 - Subtomo projections: new advance option to ignore alignments
## V3.23.03.6:
 - Map back fix: works with "other tomograms"

## V3.23.03.5:
 - Map back fix: works with 3d classes

## V3.23.03.4:
 - Mapback waits now when scheduled and not the reference is needed (validation added)

## V3.23.03.3:
 - Fit vesicles removed (moved to tomo)

## V3.23.03.2:
 - Mapback scales any reference to match the tomogram sampling.
## V3.23.03.1:
 - Subtomo projector compatible with hdf stacks. Test fixed. Projections keep orientation.

## V3.23.03.0:
 - New protocol and test extract_subtomos
 - New protocol splitTS
 - New protocol score_transform
 - New versions numering (associated to xmipp repository)
 - Updated flexalign


## V3.1.3:
 - Fix score transformations: avoid scores over 180.
 - Fix score transformations: not doing the correct match with score index.

## V3.1.2:
 - Projection subtomograms uses alignment information
 - Alignment consensus: Plot series is angle now. Output is of same type as input.

## V3.1.1:
 - First released version of the plugin documented properly in this file. The plugin contains the following 23 protocols:

   xmipptomo - align transformations ( XmippProtAlignTransform ):
      Protocol to rotate a series of alignments to a common reference defined by a
      Subtomogram Average

   xmipptomo - apply alignment subtomo ( XmippProtApplyTransformSubtomo ):
      Apply alignment matrix and produce a new setOfSubtomograms, with each subtomogram aligned to its reference.

   xmipptomo - apply alignment tilt-series ( XmippProtApplyTransformationMatrixTS ):
      Compute the interpolated tilt-series from its transform matrix.

   xmipptomo - connected components to ROIs ( XmippProtCCroi ):
      This protocol adjust a SetOfCoordinates (which usually will come from a
      connected componnent) to a ROI (region of interest) previously defined

   xmipptomo - cltomo ( XmippProtCLTomo ):
      Averages a set of subtomograms taking into account the missing edge.

   xmipptomo - connected components ( XmippProtConnectedComponents ):
      This protocol takes a set of coordinates and identifies connected
      components among the picked particles.

   xmipptomo - crop tomograms ( XmippProtCropTomograms ):
      Protocol to crop tomograms using xmipp_transform_window.
      The protocol allows to change the size of a tomogram/s, by removing the
      borders defined by the users

   xmipptomo - Filter coordinates by map ( XmippProtFilterCoordinatesByMap ):
      Filter coordinate by map both given a mask or a resolucion map from a tomogram

   xmipptomo - half maps ( XmippProtHalfMapsSubtomo ):
      Create half maps from a SetOfSubtomograms and its alignment

   xmipptomo - local Resolution MonoTomo ( XmippProtMonoTomo ):
      Given a tomogram the protocol assigns local resolutions to each voxel of the tomogram.
      To do that, thje protocol makes use of two half tomograms, called odd and even.
      These tomograms are reconstructed with the same alignment parameter but using the
      half of the data. For instance, the odd/even-images of the tilt series, or much
      better using the odd/even frames of the movies (recommended). The result is a
      tomogram with the values of local resolution.

   xmipptomo - phantom create subtomo ( XmippProtPhantomSubtomo ):
      Create subtomogram phantoms

   xmipptomo - phantom tomograms ( XmippProtPhantomTomo ):
      Create phantom tomograms with phantom particles and its coordinates with the right Scipion transformation matrix

   xmipptomo - resize tilt-series ( XmippProtResizeTiltSeries ):
      Wrapper protocol to Xmipp image resize applied on tilt-series

   xmipptomo - resize tomograms ( XmippProtResizeTomograms ):
      Protocol to to resize tomograms using xmipp_image_resize.
      The protocol allows to change the size of a tomogram/s by means
      of different methods

   xmipptomo - imagej roi ( XmippProtRoiIJ ):
      Tomogram ROI selection in IJ

   xmipptomo - score/filter coordinates ( XmippProtScoreCoordinates ):
      Scoring and (optional) filtering of coordinates based on different scoring
      functions (carbon distance, neighbour distance)

   xmipptomo - subtomo alignment consensus ( XmippProtScoreTransform ):
      Protocol to score a series of alignments stored in a SetOfSubtomograms by
      quaternion distance analysis.

      xmipp_alignmentDistance ranges from 0º to 180º. Therefore, a 0º distance is the best and means alignment is the same.
      The lower the score the more similar is the alignment.

   xmipptomo - split tilt-series ( XmippProtSplitTiltSeries ):
      Wrapper protocol to Xmipp split Odd Even on tilt-series

   xmipptomo - map back subtomos ( XmippProtSubtomoMapBack ):
      This protocol takes a tomogram, a reference subtomogram and a metadata with geometrical parameters
      (x,y,z) and places the reference subtomogram on the tomogram at the designated locations (map back).
      It has different representation options.

   xmipptomo - subtomo projection ( XmippProtSubtomoProject ):

      Project a set of volumes or subtomograms to obtain their X, Y or Z projection of the desired range of slices.

   xmipptomo - subtomo subtraction ( XmippProtSubtractionSubtomo ):
      This protocol subtracts a subtomogram average to a SetOfSubtomograms, which are internally aligned and
      numerically adjusted in order to obtain reliable results. The adjustment and subtraction is perfomed by
      xmipp_volume_subtraction program. A mask can be provided if the user wants to perform the subtraction in a
      determined region.

   xmipptomo - tiltseries FlexAlign ( XmippProtTsFlexAlign ):

      Simple protocol to average TiltSeries movies as basic
      motion correction. It is used mainly for testing purposes.
