# B4 Virtual Morphogenesis Laboratory

A static, browser-based educational simulation for **Computational Morphology — Module B4: Computational Models of Morphogenesis**.

## Quick start

Open `index.html` in a modern desktop browser, or deploy the whole folder to GitHub Pages. No dependencies, build step, API keys or backend are needed.

## Deploy to GitHub Pages

1. Create a public GitHub repository (for example, `b4-morphogenesis-lab`).
2. Upload **all files inside this folder** (`index.html`, `style.css`, `app.js`, and `assets/`) to the repository's root. Or upload the ZIP contents.
3. In the GitHub repository, go to **Settings → Pages**.
4. Under **Build and deployment**, select **Deploy from a branch**, choose `main` and `/(root)`, then save.
5. The live site will appear at `https://YOUR-USERNAME.github.io/b4-morphogenesis-lab/` after deployment.

The application uses relative paths and works in a project subdirectory.

## Student workflow

1. Choose guided colonization, geometry and fate, or viability-kernel discovery.
2. Explore pore diameter, curvature, stiffness, adhesion, motility, proliferation and scaffold geometry.
3. Run to t=300, save a result, and vary a single parameter.
4. Inspect trajectories and the time series. Compare multiple experiments in morphospace.
5. Run a 20-design sweep, adjust the outcome criteria, then export CSV or a report.

## Model summary

- 2D off-lattice stochastic cell agents moving through a simplified spatial scaffold obstacle mask.
- Neighbour short-range repulsion and adhesion, scaffold-weighted motility, directional guidance, crowding and seed-controlled noise.
- Crowding/geometry/stiffness-dependent probabilistic proliferation and differentiation.
- Synthetic extracellular-matrix deposition event dots.
- Outputs: population, spatial-bin coverage, differentiated fraction, average horizontal infiltration and matrix index.
- Parameter sweeps run the **simulation**, not an analytic score lookup.
- The green morphospace ellipse is **only a visual reference heuristic**. Simulated points are colored by user-defined viability criteria.

## Scientific integrity and limitations

The model is an **original illustrative teaching model**, informed by the concepts in Galli et al., *Morphospace Engineering: Morphological Computation in Scaffold Design*, *BioSystems*. The article's CA and ABM discussions were conceptual, not implementations. This app should **never be described as a reproduction of the article's computational results** or as experimentally validated predictions.

Physical sizes and curvature axes are approximate design labels. The graphics do not maintain metrological pore-size scaling and forces are not fitted to cell biology. The chemical cue is not a reaction-diffusion PDE. The viability classification is a user-defined multi-objective threshold and **not a rigorous viability-theory kernel**.

## Paper figures and reuse

Two cropped figures from the instructor-provided manuscript are supplied in `assets/` for optional educational viewing. If the site is to be publicly published, ensure that the figures' reuse complies with the final publication's license or obtain publisher permission, or remove `assets/paper-figure-6.png` and `assets/paper-figure-7.png` and edit the two links in `app.js`.

## Suggested classroom activity

Hypothesis: increasing pore diameter increases infiltration at the expense of clustering. Run pore diameters of 180, 300 and 470 µm with other settings and the seed held constant; save the three trials; repeat with a different seed; discuss differences, assumptions, and experiments required for biological validation.

## Files

- `index.html` — semantic page layout
- `style.css` — responsive styling
- `app.js` — simulation, visualizations, exports, pedagogical guides
- `assets/` — two instructor-supplied manuscript figure crops
