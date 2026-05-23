# Legacy MATLAB Archive

This directory preserves the original MATLAB FRAP-Toolbox implementation for
historical reproduction, parity checks, and comparison with archived exports.
The modern Python app and CLI are the primary implementation for new analysis
work.

The archived MATLAB entry point is:

```matlab
Main_GUI
```

Run it from this directory, or add this directory to the MATLAB path from the
repository root:

```matlab
addpath(fullfile(pwd, 'legacy', 'matlab'));
Main_GUI;
```

`loci_tools.jar` is kept here with the MATLAB source because the legacy
Bio-Formats helpers expect it beside the archived MATLAB files.

The old website-style README content is preserved in
`legacy-website-content.md`. It may mention paths and distribution channels from
the MATLAB era; the root README is the current project front door.
