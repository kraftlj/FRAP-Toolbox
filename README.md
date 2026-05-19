FRAP-Toolbox: Software for the analysis of Fluorescence Recovery After Photobleaching
============
Lewis J. Kraft, Jacob Dowler, and Anne K. Kenworthy
Vanderbilt University Medical Center, Nashville TN

# FRAP-Toolbox Website Content

This Markdown file reproduces the informational content from the legacy FRAP-Toolbox website so it can live in a GitHub repository. Use the table of contents below to navigate the sections.

## Table of Contents
- [Home](#home)
- [About](#about)
- [Downloads](#downloads)
- [Getting Started](#getting-started)
- [FRAP Models and Their Applications](#frap-models-and-their-applications)
- [Test Data](#test-data)
- [Troubleshooting](#troubleshooting)
- [Developer Topics](#developer-topics)
- [Recent Applications](#recent-applications)
- [Archive](#archive)
- [FAQ](#faq)
- [Bug Reporting](#bug-reporting)
- [Contact](#contact)

## Home
FRAP-Toolbox is software for the analysis of Fluorescence Recovery After Photobleaching (FRAP). The toolbox provides:

- A downloadable application for Windows and macOS systems.
- Documentation to help new users get started quickly.
- Guidance for developers who want to extend the toolbox.
- Examples of recent scientific applications that rely on the included models.

## About
### History
FRAP-Toolbox was developed to make quantitative FRAP analysis accessible to the broader research community. The project grew out of a collaboration between Anne Kenworthy's research group and biomathematicians led by Emmanuele DiBenedetto, funded by the National Science Foundation.

The team's goal is to help researchers extract accurate diffusion coefficients and binding constants from FRAP data. Many consensus models already existed, but the developers wanted to remove the steep learning curve associated with advanced image analysis and non-linear fitting routines.

### Philosophy
Researchers should be able to apply established FRAP analysis methods to their own data without having to rebuild complicated toolchains from scratch.

### Citing FRAP-Toolbox
If you cite FRAP-Toolbox, please use:

```
Kraft, LJ, Dowler J, Kenworthy AK. (2014). Frap-Toolbox: Software for the analysis of Fluorescence Recovery After Photobleaching. http://www.fraptoolbox.com (accessed Month Day, Year)
```

### License
FRAP-Toolbox is open source and released under the [GNU General Public License](http://www.gnu.org/licenses/gpl.html).

### Questions
- See the [FAQ](#faq) for answers to common questions.
- Visit the [Bug Reporting](#bug-reporting) section to report software issues.
- Email [Anne Kenworthy](mailto:anne.kenworthy@vanderbilt.edu) with additional comments or feedback.

## Downloads
### Standalone packages (no MATLAB Compiler Runtime included)

| Client      | File name                     | Size   |
|-------------|------------------------------:|-------:|
| Windows x32 | `DownloadFiles/FRAP-Toolbox_x32_v1.1.zip` | 7.9 MB |
| Windows x64 | `DownloadFiles/FRAP-Toolbox_x64_v1.1.zip` | 18.8 MB |
| macOS x64   | `DownloadFiles/FRAP-Toolbox_MAC_v1.1.zip` | 609 MB |

### Source files
You can also run FRAP-Toolbox with a full MATLAB installation by cloning the source code from GitHub: <https://github.com/kraftlj/FRAP-Toolbox>.

### Installation and documentation
- Installation instructions and test data live in the [Getting Started](#getting-started) section.
- A downloadable PDF user guide is available at `DownloadFiles/UserGuide.pdf`.

### Legacy versions
Previous releases are archived on the [FRAP-Toolbox archive](#archive) page.

## Getting Started
### Installation
You can run FRAP-Toolbox in two ways:

1. Use the MATLAB source files with a full MATLAB installation.
2. Install the standalone application after installing the free MATLAB Compiler Runtime (MCR).

The downloads page lists both the standalone installers and the source code.

#### System requirements
FRAP-Toolbox has been tested on:

- Windows XP (32-bit) and Windows 7 (64-bit).
- macOS 10.9.

Refer to the [MATLAB 2013 system requirements](http://www.mathworks.com/support/sysreq/current_release/) for MATLAB-specific details.

#### Running the source code
1. Open MATLAB.
2. Navigate to the `FRAP-Toolbox` directory that contains the source files.
3. Open and run `MainGUI.m`.

#### Installing the standalone application
##### Windows instructions
1. Download the installer from the [downloads](#downloads) section.
2. Move the `FRAP-Toolbox` folder to a location such as `C:\FRAP-Toolbox`.
3. Install the MATLAB Compiler Runtime by running `MCR_R2013a_win32_installer.exe` and follow the prompts.
4. Edit `classpath.txt` located at `C:\Program Files\MATLAB\MATLAB Compiler Runtime\v81\toolbox\local\classpath.txt`:
   - Grant yourself permission to edit the file (Properties > Security > Edit).
   - Append a new line with `C:\FRAP-Toolbox\loci_tools.jar`.
5. Launch FRAP-Toolbox by double-clicking `FrapToolbox.exe`.

> Note: You can store the `FRAP-Toolbox` folder anywhere, but `classpath.txt` must reference the correct path to `loci_tools.jar`.

##### macOS instructions
1. Download the installer from the [downloads](#downloads) section. Right-click the `FRAPToolbox_Installer` icon and choose **Open**.
2. The installer creates `/Applications/FRAPToolbox` and installs the MCR into `/Applications/MATLAB/`.
3. Move `loci_tools.jar` into `/Applications/FRAPToolbox/`.
4. Edit `/Applications/MATLAB_R2013b/toolbox/local/classpath.txt` and add `/Applications/FRAPToolbox/loci_tools.jar` on a new line. You may need to right-click `MATLAB_R2013b` and choose **Show Package Contents**. Administrative privileges may be required.
5. Launch the app from `/Applications/FRAPToolbox/application/`.

> Note: As on Windows, make sure the path to `loci_tools.jar` in `classpath.txt` matches the actual file location.

### Overview of the software
![FRAP-Toolbox overview](Images/FRAP%20Toolbox%20Overview.jpg)

The main window collects the essential inputs needed to analyze FRAP datasets. Users choose a directory containing raw data files, select one or more datasets, pick the analysis model, define the bleach geometry, indicate the frame number of the bleaching event, provide a background intensity, decide whether to normalize by the whole cell, and choose how many pre-bleach images to use.

From the main window, users can preview the first dataset before proceeding to data fitting.

### Supported image formats
FRAP-Toolbox integrates with [Bio-Formats](http://www.openmicroscopy.org/site/support/bio-formats4/supported-formats.html), a Java library for life-science image file formats. Verified formats include `.lsm` and `.nd2` files from Zeiss and Nikon microscopes.

### Using FRAP-Toolbox
![Main GUI (Figure 1)](Images/Figure1%20MainGUI_Diffusion.jpg)

- The main window lists available datasets and displays the inputs described above.
- Previewing a dataset (Figure 2) loads the image stack, provides a scrubber to browse time points, and overlays the bleaching ROI.
- Proceeding to data analysis loads all datasets through Bio-Formats and opens the analysis window (Figure 3). Users can set initial guesses, parameter bounds, and exclude points.

![Preview window (Figure 2)](Images/Figure2%20ImagePreview_Diffusion.jpg)
![Analysis window (Figure 3)](Images/Figure3%20DataAnalysis_Diffusion_after.jpg)

#### Data presentation
Running the fit produces diagnostic plots (Figure 4) showing initial conditions, FRAP curves, fits, and residuals. Optimized parameters populate a table and can be saved, along with raw and processed data, as tab-delimited text files.

![Results inspection (Figure 4)](Images/Figure4%20ResultsInspection_Diffusion.jpg)

### Considerations for designing FRAP experiments
#### Diffusion model requirements
1. Use a circular bleach ROI; record its center `(x, y)` and radius in pixels.
2. Record the frame number of the first post-bleach image.
3. Measure the mean background intensity using unlabeled controls (often approximated as zero).
4. Capture the full cell within the frame if you plan to normalize by whole-cell intensity.
5. Acquire at least one pre-bleach image for normalization.

> Note: Normalizing by the whole-cell intensity inherently corrects for imaging decay and fluorescence loss due to the bleach. If you skip normalization, account for photobleaching manually. FRAP-Toolbox can model this decay as a single exponential if you collect data until steady state.

#### Reaction 1 and Reaction 2 model requirements
1. Use either a circular bleach ROI (record center and radius) or a user-defined polygon.
2. Record the frame number of the first post-bleach image.
3. Measure the mean background intensity from unlabeled samples.
4. Capture the full cell within the frame when normalizing by whole-cell intensity.
5. Acquire at least one pre-bleach image for normalization.

## FRAP Models and Their Applications
### Diffusion model
The diffusion model simulates recoveries governed by single-component Brownian motion. It provides an analytical equation for extracting diffusion coefficients, under the assumptions of a homogeneous molecular distribution, two-dimensional diffusion (due to complete bleaching through the sample thickness), infinite boundaries, and a single diffusing component.

Let `I(t)` be the mean fluorescence intensity within the bleach region, normalized to the pre-bleach steady state. The diffusion coefficient `D` and the mobile fraction `Mf` are obtained by fitting:

$$
I(t) = I_0 \left( \sum_{m=0}^{20} \frac{-K^m r_e^2}{m! \left[r_e^2 + m (8 D t + r_n^2)\right]} \right) Mf + (1 - Mf) I(0)
$$

`I_0` equals 1 for a normalized FRAP curve, and `r_n` is the nominal bleaching radius. This equation modifies the Axelrod solution to account for a Gaussian laser profile. Parameters `r_e` and `K` derive from fitting the normalized radial post-bleach profile `I(x; t = 0)`:

$$
I(x; t = 0) = I_0 \exp \left( -K \exp \left[-\frac{2 x^2}{r_e^2} \right] \right)
$$

`D` relates to molecular and medium properties. For spherical molecules undergoing Brownian motion:

$$
D = \frac{k_B T}{6 \pi \eta R}
$$

`k_B` is Boltzmann's constant, `T` is absolute temperature, `\eta` is viscosity, and `R` is the particle radius. The mobile fraction `Mf` quantifies the percentage of molecules free to diffuse.

#### Corrections
To correct for photobleaching during imaging, divide by the integrated whole-cell intensity, or approximate the decay as a single exponential after returning to steady state:

$$
I(t) = e^{-k_{decay} t}
$$

`k_{decay}` is the photobleaching rate constant. To correct mobile fractions when fluorescence is lost from the compartment, measure an adjacent ROI and compute:

$$
Mf_{correct} = 1 - \left(I_{adjacent}(\infty) - I(\infty)\right)
$$

#### Curve fitting parameters
- `D`
- `Mf`

#### References
1. Axelrod D, Koppel DE, Schlessinger J, Elson E, Webb WW (1976) Mobility measurement by analysis of fluorescence photobleaching recovery kinetics. *Biophysical Journal* 16: 1055-1069.
2. Braga J, Desterro JM, Carmo-Fonseca M (2004) Intracellular macromolecular mobility measured by FRAP with confocal laser scanning microscopes. *Molecular Biology of the Cell* 15: 4749-4760.
3. Kang M, Day CA, Drake K, Kenworthy AK, DiBenedetto E (2009) A generalization of theory for two-dimensional FRAP applicable to confocal laser scanning microscopes. *Biophysical Journal* 97: 1501-1511.

### Reaction 1 model
The Reaction 1 model uses a single-component exponential:

$$
I(t) = a - b e^{-c t}
$$

Consider molecules free to diffuse (`f1`) or bound in an immobile complex (`c1`):

```
f1 <-> c1
```

Assuming `f1` equilibrates rapidly within the bleach region (`f1 = F_eq`), the concentration of the complex over time follows:

$$
\frac{d c_1}{d t} = k^*_{on} F_{eq} - k_{off} c_1
$$

The FRAP curve becomes:

$$
I(t) = I(\infty) - \left[I(\infty) - I(0)\right] C_{eq} e^{-k_{off} t}
$$

`C_{eq} = (k^*_{on} / k_{off}) F_{eq}`. Comparing with the generic form above gives `a = I(\infty)`, `b = I(\infty) - I(0)`, and `c = k_{off}`.

#### Curve fitting parameters
- `a`
- `b`
- `c`

#### References
1. Mueller F, Wach P, McNally JG (2008) Evidence for a common mode of transcription factor interaction with chromatin revealed by quantitative FRAP. *Biophysical Journal* 94: 3323-3339.
2. Wei X, Henke VG, Strub C, Brown EB, Clapham DE (2003) Real-time imaging of nuclear permeation by EGFP. *Biophysical Journal* 84: 1317-1327.

### Reaction 2 model
The Reaction 2 model uses a two-component exponential:

$$
I(t) = a - b e^{-c t} - d e^{-f t}
$$

Example system with two immobile complexes (`c1`, `c2`):

```
f1 <-> c1
f1 <-> c2
```

Assuming `f1 = F_eq`, the concentrations follow:

$$
\frac{d c_1}{d t} = k^*_{1on} F_{eq} - k_{1off} c_1
$$
$$
\frac{d c_2}{d t} = k^*_{2on} F_{eq} - k_{2off} c_2
$$

The FRAP expression becomes:

$$
I(t) = I(\infty) - \left[I(\infty) - I(0)\right] C_{1eq} e^{-k_{1off} t} - \left[I(\infty) - I(0)\right] C_{2eq} e^{-k_{2off} t}
$$

where:

$$
C_{1eq}^{-1} = 1 + \frac{k_{1off}}{k^*_{1on}} \left(1 + \frac{k^*_{2on}}{k_{2off}}\right)
$$
$$
C_{2eq}^{-1} = 1 + \frac{k_{2off}}{k^*_{2on}} \left(1 + \frac{k^*_{1on}}{k_{1off}}\right)
$$

Comparing terms yields `a = I(\infty)`, `b = [I(\infty) - I(0)] C_{1eq}`, `c = k_{1off}`, `d = [I(\infty) - I(0)] C_{2eq}`, and `f = k_{2off}`.

#### Curve fitting parameters
- `a`
- `b`
- `c`
- `d`
- `f`

#### References
1. Benesh AE, Nambiar R, McConnell RE, Mao S, Tabb DL, et al. (2010) Differential localization and dynamics of class I myosins in the enterocyte microvillus. *Molecular Biology of the Cell* 21: 970-978.
2. Mueller F, Wach P, McNally JG (2008) Evidence for a common mode of transcription factor interaction with chromatin revealed by quantitative FRAP. *Biophysical Journal* 94: 3323-3339.

## Test Data
Fully documented examples and datasets accompany each model.

### Diffusion model workflow
Using a Zeiss LSM 510 confocal microscope, a 1 um circular region in the cytoplasm of HeLa cells expressing Venus or Venus-ATG5 was photobleached.

1. Download the datasets at `DownloadFiles/Diffusion.zip`.
2. Open FRAP-Toolbox and select the **Diffusion** model.
3. Choose **Circle** for the ROI.
4. Set the post-bleach image to `21`.
5. Set the mean background intensity to `0`.
6. Disable whole-cell normalization.
7. Use `10` pre-bleach images for normalization.
8. Select all datasets named `Venus_Cytoplasm_*.lsm`.
9. Click **Next**.
10. When prompted, provide the bleach ROI center `256 23` and radius `9`.
11. Choose **Yes** to calculate a corrected mobile fraction.
12. After loading, open the analysis window and configure parameters as in `Images/Table1%20Curve%20fitting%20parameters%20for%20the%20Diffusion%20test%20data.jpg`.
13. Select **No** when asked to fit averaged data.
14. Ensure all datasets remain selected and press **Run**.
15. Review the output plots and compare optimized parameters with `Images/Table2%20Optimized%20Curve%20fitting%20parameters%20for%20the%20Diffusion%20test%20data.jpg`.
16. Save parameters and datasets using a descriptive tag such as `Venus_Cytoplasm`. FRAP-Toolbox appends `_Diffusion_Fit_Parameters.txt`, `_Diffusion_FRAP_datasets.txt`, and `_Diffusion_Postbleach_profiles.txt` to the chosen tag.

### Reaction 1 model workflow
HeLa cells expressing Venus were photobleached in the nucleus with a user-defined ROI using a Nikon Eclipse Ti confocal microscope.

1. Download the datasets at `DownloadFiles/Reaction 1.zip`.
2. Open FRAP-Toolbox and select the **Reaction 1** model.
3. Choose **User Defined** ROI.
4. Set the post-bleach image to `6`.
5. Set the mean background intensity to `0`.
6. Enable whole-cell normalization.
7. Use `5` pre-bleach images for normalization.
8. Select datasets named `Venus_*.nd2`.
9. Click **Next** and draw the bleach ROI as shown in `Images/SuppFigure5%20Draw%20ROIs%20ATG5.jpg` (panel A).
10. Draw the whole-cell ROI (same figure, panel A).
11. Configure parameters as in `Images/Table3%20Curve%20fitting%20parameters%20for%20the%20Reaction%201%20test%20data.jpg`.
12. Choose **No** when asked to fit averaged data.
13. Ensure all datasets are selected and press **Run**.
14. Review plots and compare optimized parameters with `Images/Table%204%20Optimized%20Curve%20fitting%20parameters%20for%20the%20Reaction%201%20test%20data..jpg`.
15. Save results with a tag such as `Venus_NCTransport`; files receive `_Reaction_Fit_Parameters.txt` and `_Reaction_FRAP_datasets.txt` suffixes.

### Reaction 2 model workflow
HeLa cells expressing Venus-ATG5 were photobleached in the nucleus with a user-defined ROI.

1. Download the datasets at `DownloadFiles/Reaction 2.zip`.
2. Open FRAP-Toolbox and select the **Reaction 2** model.
3. Choose **User Defined** ROI.
4. Set the post-bleach image to `6`.
5. Set the mean background intensity to `0`.
6. Enable whole-cell normalization.
7. Use `5` pre-bleach images for normalization.
8. Select datasets named `Venus-Atg5_*.nd2`.
9. Click **Next** and draw the bleach ROI as shown in `Images/SuppFigure5%20Draw%20ROIs%20ATG5.jpg` (panel B).
10. Draw the whole-cell ROI (panel B).
11. Configure parameters as in `Images/Table%205%20Curve%20fitting%20parameters%20for%20the%20Reaction%202%20test%20data..jpg`.
12. Choose **Yes** when asked to fit averaged data.
13. Confirm all datasets are selected and press **Run**.
14. Review plots and compare optimized parameters with `Images/Table%206%20Optimized%20Curve%20fitting%20parameters%20for%20the%20Reaction%202%20test%20data..jpg`.
15. Save outputs using a tag such as `Venus-Atg5_NCTransport`; files receive `_Reaction_Fit_Parameters.txt` and `_Reaction_FRAP_datasets.txt` suffixes.

## Troubleshooting
FRAP-Toolbox displays warning dialogs when it detects potential problems, such as mismatched acquisition settings during batch processing. If you encounter reproducible issues, gather the information outlined in the [Bug Reporting](#bug-reporting) section and contact the team.

## Developer Topics
Interested in contributing or porting FRAP-Toolbox to Python? Reach out to the team and explore the source code on GitHub: <https://github.com/kraftlj/FRAP-Toolbox>.

### Python port (experimental)
- The Python implementation lives under `frap_toolbox_py/` and currently supports the diffusion workflow.
- The preferred image I/O layer is BioIO with the `bioio-tifffile` plugin. Optional extras add ND2, Bio-Formats, and legacy AICSImageIO readers.
- Large local microscopy fixtures belong in ignored `test-data/`; exploratory porting scripts belong in ignored `scratch/`.
- Install the modern app stack with Python 3.10+:
   ```bash
   python -m pip install -e ".[app,test]"
   ```
- Launch the local browser app:
   ```bash
   frap-toolbox-app
   ```
- Run the diffusion workflow from the command line:
   ```bash
   frap-toolbox /path/to/dataset1.lsm --roi 256 23 9 --post-bleach-frame 21
   ```
- Use `--fit-mode simplified_kang` for the Kang et al. per-curve half-time
  diffusion estimator, or `--fit-mode simplified_kang_global` for a pooled
  global fit of the simplified profile/recovery equations.
- Reaction 1 and Reaction 2 backend fitters are available through
  `--model reaction1` and `--model reaction2`; raw guide parity still needs
  stored user-defined bleach/cell ROI masks for the ND2 workflows.

## Recent Applications
Below are recent publications that applied FRAP-Toolbox algorithms. Each image is available in the `Images/` directory.

- Kraft LJ, Nguyen TA, Vogel SS, Kenworthy AK (2014) Size, stoichiometry, and organization of soluble LC3-associated complexes. *Autophagy* 10. [Learn more](https://www.landesbioscience.com/journals/autophagy/article/28175/)
- Kraft L, Kenworthy A (2012) Imaging protein complex formation in the autophagy pathway. *Journal of Biomedical Optics* 17: 11008. [Learn more](http://biomedicaloptics.spiedigitallibrary.org/article.aspx?articleid=1182818)
- Day CA*, Kraft LJ*, Kang M, Kenworthy AK (2012) Analysis of protein and lipid dynamics using confocal FRAP. *Current Protocols in Cytometry*. [Learn more](http://onlinelibrary.wiley.com/doi/10.1002/0471142956.cy0219s62/full)
- Kang M, Day CA, Kenworthy AK, DiBenedetto E (2012) Simplified equation to extract diffusion coefficients from confocal FRAP data. *Traffic*. [Learn more](http://onlinelibrary.wiley.com/doi/10.1111/tra.12008/abstract)
- Kang M, Day CA, DiBenedetto E, Kenworthy AK (2010) A quantitative approach to analyze binding diffusion kinetics by confocal FRAP. *Biophysical Journal*. [Learn more](http://www.sciencedirect.com/science/article/pii/S0006349510011148)
- Kang M, Day CA, Drake K, Kenworthy AK, DiBenedetto E (2009) Generalized theory for two-dimensional FRAP. *Biophysical Journal*. [Learn more](http://www.sciencedirect.com/science/article/pii/S0006349509011229)
- Day CA, Kenworthy AK (2012) Mechanisms underlying the confined diffusion of cholera toxin B-subunit. *PLoS ONE*. [Learn more](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034923)

## Archive
Older FRAP-Toolbox versions remain available for download:

- Windows x32: `DownloadFiles/FRAP-Toolbox_x32.zip`
- Windows x64: `DownloadFiles/FRAP-Toolbox_x64.zip`
- macOS x64: `DownloadFiles/FRAP-Toolbox_MAC.app.zip`

## FAQ
1. **Is FRAP-Toolbox free to use?** Yes. It is open source under the GNU GPL. Please cite the software if it supports your work.
2. **How do I cite FRAP-Toolbox?** Use the citation listed in the [Citing FRAP-Toolbox](#citing-frap-toolbox) section.
3. **macOS installation warning:** If you see "FRAPToolboxMac can't be opened because it is from an unidentified developer," right-click the installer, choose **Open**, or hold **Control** while clicking to bypass Gatekeeper restrictions.

## Bug Reporting
Before submitting a bug report, make sure you are running the latest release. For reproducible issues, gather:

- Exact error messages (copy and paste the text).
- Version information for FRAP-Toolbox and your operating system.
- Non-working datasets (send files, or request FTP upload instructions if files are large).
- Additional details and screenshots of the settings that triggered the issue.

Email the compiled information to [Lewis Kraft](mailto:kraftlj@gmail.com).

## Contact
**Address**

Vanderbilt University School of Medicine  
Dept. of Molecular Physiology and Biophysics  
1211 Medical Center Drive  
718 Light Hall 0615  
Nashville, TN 37232  
Phone: 615.322.6615

**Email**

- Support: [kraftlj@gmail.com](mailto:kraftlj@gmail.com)
- Correspondence: [anne.kenworthy@vanderbilt.edu](mailto:anne.kenworthy@vanderbilt.edu)
