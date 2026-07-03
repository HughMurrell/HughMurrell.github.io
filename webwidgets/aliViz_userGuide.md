# aliViz User Guide

aliViz is a bioinformatics alignment and phylogeny viewer. It supports loading alignments, inferring or importing trees, grouping sequences, clustering (tree-based and projection-based), epitope definition and logos, and 3D structure viewing with epitope coloring.

---

## 1. Loading and view controls

### Choose file
- **Function:** Load an alignment from a **FASTA** or **FASTQ** file.
- **Usage:** Click **Choose file**, select your alignment. The loaded filename appears in the text box beside the button (replacing the initial “No file chosen” placeholder).
- **Reference:** The **first** sequence is always treated as **Reference (REF)**.
- **Subtype:** If **Has subtype** is checked (see below), the **second** sequence is treated as **SubType**. If unchecked, no sequence is designated as subtype unless you add one later by other means.

### Has subtype
- **Function:** Tell aliViz whether the alignment includes a dedicated subtype sequence in position 2 (after reference).
- **Default:** **Unchecked** (many alignments have no subtype).
- **When checked:** Sequence 2 is the subtype; it receives the **[Subtype]** label and is excluded from tree inference (unless it is also the founder).
- **When unchecked:** No subtype index is assigned; subtype labeling and subtype-specific clustering behaviour are disabled.
- **After load:** You can toggle this checkbox. Changing it **clears grouping, tree inference, and clustering** and resets related state so you can re-run those steps with the new subtype setting.

### Character sanitization (FASTA/FASTQ load)
- **Allowed characters:** Uppercase letters **A–Z**, gap **`-`**, and stop codon **`*`** (common in amino acid alignments).
- **Non-standard symbols:** Any other character in a sequence line is **replaced with `X`** (alignment length is unchanged). You receive an alert listing what was replaced (e.g. `? → X (12)`).
- **Nucleotide alignments:** Standard IUPAC nucleotide letters are within A–Z and are kept as-is.

### View Mode: NT / AA
- **Function:** Switch between **Nucleotide (NT)** and **Amino acid (AA)** view.
- **Algorithm:** In AA mode, nucleotides are translated using the standard genetic code (frame can be set if frame selector is enabled). Gaps and invalid codons produce gap characters or `X`.
- **Modes:** **NT** (default), **AA**.

### AA palette
- **Function:** Choose the amino acid color scheme.
- **Options:** *Alignment (IUPAC)* (default) or *Biochemical* (e.g. by property).

### Highlighter
- **Function:** Highlight differences in the alignment with respect to a chosen reference.
- **Options:** **Off** (default), **Founder**, **SubType**, **Reference**. Mismatches against the selected sequence are highlighted.
- **Founder:** Works for both a **consensus** founder (`consensus_of_…` sequence) and a **medoid** founder (an existing sequence marked as founder; see §2).

---

## 2. Grouping and sequences

### Group
- **Function:** Assign sequences to groups using a delimiter and a field in the sequence name.
- **Algorithm:** You specify a **delimiter** (e.g. `_`) and a **field number** (1-based). Each sequence name is split by the delimiter; the value at that field becomes the group label. REF, SubType (if present), and PDB chain sequences are excluded from grouping. Groups are assigned unique IDs and used for coloring and legend.
- **Result:** `state.sequenceGroups` (name → group ID) and `state.groupNames` (group ID → label) are set; sequence names are colored by group in the name panel.

### Sort
- **Function:** Reorder sequences by current grouping (or by name if no grouping), keeping REF first and SubType second (when a subtype exists).
- **Algorithm:** REF and SubType stay at the top; other sequences are sorted by group ID (from `state.sequenceGroups`), then by name within each group. PDB chains and founder are handled in the sort order.

### Prune
- **Function:** Remove sequences that belong to selected groups.
- **Usage:** Open the Prune overlay, select groups to **remove**, then Apply. Removed sequences are deleted from the alignment. Group IDs are renumbered to 0, 1, 2, … for the remaining groups.

### Add Founder
- **Function:** Designate a **founder** sequence for a chosen group (or from the subtype sequence when offered in the group list).
- **Founder source (default: Medoid):**
  - **Consensus of a group:** Builds a **new** consensus sequence (majority character per column, including gaps) from all sequences in the group (excluding REF and subtype). It is inserted into the alignment with a name such as `consensus_of_<groupLabel>` (or `consensus_of_subtype` when using the subtype option in consensus mode).
  - **Medoid sequence of a group (default):** Selects an **existing** sequence in the group that minimizes total **padded Hamming distance** to all other sequences in that group. The sequence **keeps its original name**; it is marked as founder and shown with the **[Founder]** label in the name column and tree.
- **Switching modes:** Replacing a consensus founder removes or overwrites the synthetic `consensus_of_` row; switching to medoid on a group removes a previous synthetic founder if present and tags the medoid sequence instead.
- **Clustering:** The founder is **not** forced into cluster 0. When **no subtype** is present, cluster IDs start at **1** (cluster 0 is reserved for subtype when it exists). A medoid founder participates in clustering like any other sample sequence; the **[Founder]** label is preserved after clustering (name matching tolerates `_cl-*` suffixes added by clustering).
- **Usage:** Group sequences first (for group-based founders). Open **Add Founder**, choose **Founder source**, select a group (or **[SubType]** when available), then **Add Founder**. Tree inference is cleared; re-infer or load a tree afterward.

---

## 3. Tree inference and manipulation

### Infer
- **Function:** Build or load a phylogeny for the current alignment (excluding REF, subtype when present, and PDB chains; **including** founder sequences, including medoid founders).
- **Methods:**
  - **FastTree** (default in the infer dialog): Uses FastTree (via bioWASM) on the alignment. Model options appear for NT and AA.
  - **Neighbor Joining (NJ):** Uses PhyloTools (pairwise distances → NJ tree).
  - **Load inferred tree (Newick file):** Opens a file picker for an existing **Newick** tree. Accepted extensions: **`.nwk`**, **`.newick`**, **`.tree`**, **`.treefile`**, **`.txt`**.
- **Loading an external tree:** Leaf names in the Newick file must match alignment sequence names after normalization (cluster suffixes `_cl-*` are ignored for comparison). **Reference**, **subtype** (if present), and **PDB chains** are excluded from the required leaf set. A medoid founder keeps its sample name in the alignment, so tree leaves should use that name—not a separate `consensus_of_` name. If names do not match, loading is aborted with a message listing missing or extra leaves.
- **Note:** Before inference (NJ/FastTree), any existing clustering is cleared and cluster suffixes (`_cl-*`) are removed from sequence and tree node names.

### Reroot
- **Function:** Reroot the tree on the **Founder** sequence.
- **Targets:** **Founder** (requires a founder to be defined—consensus or medoid).
- **Algorithm:** The founder leaf is located and the tree is rerooted at the edge leading to it so that the founder is the outgroup. Branch lengths and topology are preserved; only the root position changes.

### Ladderize
- **Function:** Reorder the tree and the sequence list so that the tree is ladderized and sequences match leaf order.
- **Modes:**
  - **By Depth** (default): At each node, compute the maximum root-to-leaf distance in each child’s subtree (using branch lengths). Sort children by that depth (ascending). Shallower subtrees appear first.
  - **By Weight:** At each node, sort children by number of leaves (ascending). Lighter subtrees appear first.
- **Result:** Tree drawing order and `state.viewSequences` are updated to follow the new leaf order. Histogram, Cluster, and tree export (Newick **tree** and **svg**) are enabled after ladderizing.

### Histogram
- **Function:** Show a histogram of **root-to-leaf distances** (sum of branch lengths from root to each leaf).
- **Algorithm:** For each leaf, the path from root to leaf is traversed and branch lengths (`node.len`) are summed. Distances are binned (number of bins between 10 and 30, about √n). Bars are drawn for each bin count.
- **Requirement:** Tree must be loaded and ladderized.

### tree (Export tree)
- **Function:** Export the current tree in **Newick** format.
- **Usage:** Click “tree” to download a Newick file containing the current tree with names and branch lengths (if present).

### svg (Export tree as SVG)
- **Function:** Download a **vector (SVG)** figure of the current tree for publications or slides.
- **Requirement:** A tree must be loaded and **ladderized** (same as the Newick export).
- **Usage:** Click **svg** in the phylogeny toolbar. A dialog offers **Geometry** and **Scale**, then **Export** to save `tree.svg`. **Cancel** closes without downloading.

#### Geometry
- **Linear:** Rectangular cladogram (horizontal branches, vertical backbone at internal nodes, diamond tips, dashed connectors to labels).
- **Circular:** Radial layout with **orthogonal** branches (circumferential arc + radial segment), outward-pointing tip diamonds, and dotted connectors to radial labels.

#### Scale
The plot has a **520 px** branch region (≈ **13.76 cm** at 96 px/inch) available for branch lengths; its half-way mark is at **260 px**.

- **Auto (default):** Scales the tree to fit the plot area (linear branch span ≈ 520 px; circular outer branch radius ≈ 260 px, with branch lengths in circular mode drawn at **half** the linear scale per unit).
- **Fixed:** Uses a chosen **branch length per 1 cm** of plot (dropdown: 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001; default **0.001**) so trees from different datasets are comparable. The tree is **not** shrunk to fit; shallow trees use only part of the space. If the tree would extend beyond the fixed plot width, export **aborts** with a message suggesting a **larger** value in the dropdown (more compact drawing) or **Auto** scale.
- **DSI (days since infection):** Scales so an **expected** branch length lands at the **half-way mark (260 px)** of the branch region. You supply two parameters:
  - **DSI (days):** number of days since infection (default **14**).
  - **Mutations/day (MPD):** expected substitution rate per day (default **0.001**).
  - The expected branch length is **DSI × MPD** (e.g. 14 × 0.001 = 0.014), and the scale is set so this value maps to 260 px (i.e. `px per unit = 260 / (DSI × MPD)`). In circular mode the same value maps to a **130 px** radius (half, matching the Fixed circular convention).
  - **No overflow abort:** unlike Fixed scale, DSI **never** aborts. Branches longer than the plot are **clipped at the right border** (so they cannot run into the sequence-name column); nothing is clipped on the left. In circular mode over-long branches are clamped to the outer radius.
  - **Expected-depth marker:** a **gray vertical dotted line** is drawn across the plot at the expected depth (the end of the scale bar, 260 px). Tips to the **right** are mutating **faster** than expected; tips to the **left** are mutating **slower** than expected.

#### Titles, scale bar, and legend (SVG only)
- **Title (two lines):** Full alignment **filename** (no truncation); second line includes **inference method** (or load method), **max tree depth** (3 significant figures), and **last clustering method** (or “No clustering”). Linear layout: titles are **left-aligned** with the scale bar; circular: titles are **centred**.
- **Phylogenetic scale bar:** Drawn on the **tree** panel (top-left for linear, top-right for circular—not in the floating HTML legend). Shows branch-length units; in **fixed** mode the bar is **1 cm** long with a numeric label matching the selected scale (circular bar length is **half** the linear bar for the same scale value). In **DSI** mode the bar is **260 px** long (linear; 130 px circular) and is annotated with the parameters and expected depth, e.g. **`(14 days @ 0.001 = 0.014)`**, updating to match whatever DSI and MPD you enter.
- **Legend panel:** Header **Legend** (not “Color Legend”). **Groups** and **Clusters** as in the app; cluster rows show the **number only** (e.g. `1`, not `Cluster 1`). Page width split **7/8** tree, **1/8** legend.

#### Layout details (unchanged behaviour, for reference)
- Tip markers use **cluster** colors where clustering is active; label and connector colors use **group** / special-sequence rules (magenta for REF, subtype, founder, PDB as in the app).
- Circular plot height grows as needed for long radial labels.

### Clear Tree
- **Function:** Remove the current tree. Clustering state is cleared. Histogram, Cluster, and related buttons are disabled.

---

## 4. Clustering (Cluster button)

Clicking **Cluster** opens a method selector, then the corresponding clustering interface.
Each method clusters sample sequences; **Reference**, **Subtype** (if present), and **PDB chains** are excluded from clustering assignments.
- **Methods (selector):** **None** (default), **Hierarchical (tree-clade clustering)**, **Hierarchical (tree-cut clustering)**, **UMAP**, **MDS**.

### None
- **Function:** Remove active clustering.
- **Effect:** Clears cluster state, strips `_cl-*` suffixes from sequence and tree node names, and updates the legend. **Group** colors and group membership are preserved.

### Auto button behaviour (optimizes selected index)
- **Radio buttons (Calinski vs Ball-Hall-Adapted):** choose which index the dialog’s **Auto** button maximizes while it searches slider values. This does not change the clustering method itself; it only changes the scoring metric used by Auto (the radio choice is read when you press **Auto**). Default is **Calinski–Harabasz**; switch to **Ball-Hall-Adapted** to optimize that index instead.
- **Tree-clade Auto:** searches the branch-length threshold (step to max) to maximize the **selected index** for the current **min leaves per cluster**, then sets the slider to the best threshold.
- **Tree-cut Auto:** performs a full grid search over depth triples (d1, d2, d3) with d1 ≤ d2 ≤ d3 (sampled from the slider range) and sets the sliders to the triple that maximizes the **selected index**.
- **MDS/UMAP eps Auto:** varies the **radius (eps)** from the slider minimum (step) to max and selects the eps that maximizes the **selected index** for the current **min neighbors**.

### 4.1 Hierarchical (tree-clade clustering)
- **Function:** Define clusters by **clades** on the tree: a new clade (and thus a new cluster) starts when the **incoming branch length** to a node exceeds a threshold.
- **Parameters:**
  - **Branch length threshold:** Minimum branch length that starts a new clade (slider from step to max). Leaves are assigned to the clade that contains them.
  - **Min leaves per cluster:** Clusters with fewer than this many leaves are relabelled as noise (`_cl-na`).
- **Algorithm:** DFS from root. When traversing an edge longer than the threshold, increment cluster ID. Assign each leaf to the current cluster. Subtype (if present) is cluster 0; without subtype, numbering starts at **1**. Then apply min-leaves filter, mark small clusters as noise, and renumber clusters 1, 2, 3, …
- **Auto:** Searches the threshold (from step to max) to **maximize the Calinski–Harabasz index** (see §4.6) for the current min-leaves. Sets the slider to the best value.
- **Accept:** Applies the clustering: adds `_cl-<id>` or `_cl-na` to sequence and tree node names, updates the legend and tree colors.
- **Cancel:** Removes all clustering: clears `state.leafClusters`, strips `_cl-*` from names and tree nodes, updates legend and redraws. **Group colors** on sequence names are preserved (group mapping is rekeyed to the stripped names).

### 4.2 Hierarchical (tree-cut clustering)
- **Function:** Cut the tree with **three depth lines** (root-to-node distance). Each region between/above lines defines clusters.
- **Parameters:**
  - **Level 1, 2, 3 Depth:** Sliders (0 to max root-to-leaf distance). Order is enforced: Level 1 ≤ Level 2 ≤ Level 3.
  - **Min leaves per cluster:** As in tree-clade; clusters below this size become noise.
- **Algorithm:** Leaves with depth &lt; Level 1 = cluster 1. BFS from root: when a node’s depth crosses Level 1 (or 2 or 3), all leaves in its subtree with depth ≥ that level are assigned the next cluster ID. Then min-leaves → noise, renumber clusters.
- **Auto:** **Full grid search** over all triples (d1, d2, d3) with d1 ≤ d2 ≤ d3 (sampled from slider range). For each triple, computes the cluster map and the Calinski–Harabasz index; chooses the triple that maximizes it and sets the three sliders.
- **Accept / Cancel:** Same idea as tree-clade: Accept writes cluster tags and updates legend/tree; Cancel clears clustering and strips `_cl-*` while keeping group colors.

### 4.3 UMAP
- **Function:** Reduce pairwise leaf distances to 2D with **UMAP**, then cluster in 2D with DBSCAN (same as MDS after projection).
- **Algorithm:** UMAP (via `umap-js`) is run on the distance matrix (or a derived affinity matrix). Resulting 2D coordinates are then passed to the same DBSCAN + renumbering + Calinski–Harabasz pipeline as MDS.
- **Parameters:** **nNeighbors**, **Spread**, **Min Distance** (UMAP), plus **Radius (eps)** and **Min Neighbors** (DBSCAN). Eps min = step to avoid CH infinity; **Auto** on eps optimizes CH over the eps range.

### 4.4 MDS (Classical Multidimensional Scaling)
- **Function:** Reduce **pairwise leaf distances** (from the tree) to 2D, then cluster points in 2D with **DBSCAN**.
- **Projection:** Pairwise distances are taken from the tree (path length between leaves). **Classical MDS:** D² is double-centered to form **B** = −0.5 · H · D² · H, where **H** = I − (1/n)·1·1ᵀ (identity minus n⁻¹ times the matrix of ones). Top two eigenvalues and eigenvectors of B are computed; coordinates are the eigenvectors scaled by √λᵢ.
- **Clustering:** DBSCAN on the 2D points (see §4.5). The subtype point is excluded from DBSCAN expansion and effectively not part of the density-based clustering (it is treated as an always-isolated reference). Clusters are renumbered by distance from the subtype reference (and, if a founder exists, may be used for ordering). Calinski–Harabasz is computed from the **tree** using the same cluster assignment (so all methods are comparable).
- **Parameters:** **Radius (eps)** and **Min Neighbors** for DBSCAN. Eps slider range is step to max (no 0) to avoid degenerate CH. **Auto** on eps: varies eps from step to max, maximizes CH, updates slider and plot.
- **Apply:** Applies the clustering to the tree (writes cluster tags to names and tree, updates legend). **Close:** Dismisses the overlay without applying.

### 4.5 DBSCAN (used in MDS/UMAP)
- **Algorithm:** Standard DBSCAN. Points within **eps** (Euclidean in 2D) are neighbors. If a point has ≥ **minPts** neighbors, it and all density-reachable points form a cluster. Otherwise it is noise (-1).

### 4.6. Calinski–Harabasz index (tree-based)

Used to score and optimize clusterings (tree-clade, tree-cut, MDS/UMAP). All use the **same** index on the **tree**.

- **Definitions:**  
  - **k** = number of clusters, **n** = number of leaves in those clusters.  
  - **Between-cluster (B):** For each cluster, take the MRCA of its leaves; *d*<sub>MRCA</sub> = branch distance from MRCA to root. Then B = Σ<sub>c</sub> n<sub>c</sub> · *d*<sub>MRCA</sub>² (for k ≥ 2; for k = 1 a synthetic B is used: (n/2)·(halfMax)² where halfMax = half the longest branch in the tree).  
  - **Within-cluster (W):** For each leaf, *d*<sub>leaf</sub> = distance to root; W = Σ over leaves of (*d*<sub>leaf</sub> − *d*<sub>MRCA</sub>)² for that leaf’s cluster.

- **Formula:**  
  **CH** = [ (B/(k−1)) / (W/(n−k)) ] × 1/(2<sup>k</sup>)  
  (for k = 1 the numerator is B and denominator W/(n−1); then the same 1/2<sup>k</sup> factor). The 2<sup>k</sup> term biases against large k.

- **Special cases:** If W ≤ 0, CH is returned as ∞ (best). Noise (−1) is excluded from clustering; remaining cluster IDs are used in the CH computation.

### 4.7. Ball–Hall-Adapted index (tree-based)

This index is displayed alongside the Calinski–Harabasz value in the tree-based hierarchical clustering dialogs.

- **Ball–Hall dispersion for the current clustering (k clusters):**  
  For each cluster c, let MRCA(c) be the most recent common ancestor of that cluster’s leaves, and let `dist(x, y)` be the branch-length distance between nodes. Each leaf contributes the squared distance from the leaf to its cluster’s MRCA:
  - `dist(leaf, MRCA(c)) = dist(leaf, root) − dist(MRCA(c), root)`
  - **BH(k) = Σ<sub>c</sub> Σ<sub>leaf in c</sub> dist(leaf, MRCA(c))²**

- **Ball–Hall dispersion for k=1 (whole tree):**  
  When there is only one cluster, MRCA(c) is the tree root, so:
  - **BH(1) = Σ<sub>leaf</sub> dist(leaf, root)²**

- **Ball–Hall-Adapted value (requested adaptation):**  
  The adapted index compares dispersion for the full tree to dispersion for the current clustering, then applies the same **2<sup>k</sup>** penalty used for Calinski–Harabasz:
  - **BH<sub>adapted</sub> = [ BH(1) / BH(k) ] / 2<sup>k</sup>**

If BH(k) = 0, the adapted index is returned as ∞.

---

## 5. Epitopes

### Load Epitopes
- **Function:** Load epitope definitions from a **CSV** file.
- **CSV format:** One row per epitope. First column = epitope **name**. Remaining columns = regions: either `start:end` or a single position (same start and end). Example:

  `VRC01,197:198,230,276,278:282,365:371,427:428,430,455:463,465,469,471:474`  
  `CAP256,156:163,166:167,169:170,178:179,181:184`
- **Result:** Epitopes are stored in `state.epitopes`. If the alignment was in NT mode, it is converted to AA (and mode switched to AA). Select Epitope and Export Epitopes are enabled.
- **Merge behavior:** Loading new epitopes **merges** with existing ones: same name overwrites, other names are kept.

### Select Epitope
- **Function:** Choose which epitope is active for display (e.g. alignment columns and logo). Option “None” shows all columns.
- **Usage:** Opens a dropdown of loaded epitopes; selection restricts the visible alignment to that epitope’s regions.

### Show Logo
- **Function:** Generate a **sequence logo** for the currently selected epitope (and selected sequences).
- **Requirement:** An epitope must be selected. Logo shows conservation per position in the epitope regions.

### Export Epitopes
- **Function:** Download the current set of epitopes as **epitopes.csv**.
- **Format:** Same CSV as load (name, then region columns as `start:end` or single positions). Only epitopes with at least one region are exported.
- **State:** Disabled when there are no epitopes.

### New epitope (3D)
- **Function:** Define an epitope from the loaded 3D structure by **distance**: residues on an “intrinsic” chain whose minimum distance to an “extrinsic” chain is below a cutoff.
- **Parameters:** Intrinsic chain, extrinsic chain, **distance cutoff (Å)**. Name is auto-generated as `<pdbCode>_<cutoff>A`.
- **Result:** Residues are mapped to alignment columns (via reference); epitope is added to `state.epitopes` and selected. Select Epitope and Export Epitopes are enabled if they were not already.

---

## 6. PDB structure

### Load file / Fetch
- **Function:** Load a PDB/CIF from file or fetch by **PDB ID** from RCSB. Chains are added as sequences (e.g. `xxxx_Chain_A [PDB_A]`) and aligned to the reference in AA mode.
- **Usage:** “Choose file” or enter a 4-character ID and click Fetch.

### View 3D
- **Function:** Open the 3D viewer (3Dmol.js). Chains can be toggled and styled. Epitope residues (from loaded or 3D-defined epitopes) are colored on intrinsic chains; rest can be gray.
- **Chain categories:** The 3D viewer categorizes chains into three groups:
  - **Intrinsic (≥50%)**: Chains present in the alignment panel that align at **≥50%** coverage/identity to the alignment reference.
  - **Extrinsic (<50%)**: Chains present in the alignment panel that align at **<50%** coverage/identity.
  - **Superficial**: Chains **not** present in the alignment panel.
- **Epitope residue coloring** is applied to intrinsic chains; clicking on an active epitope residue initiates a **partial logo** dialogue; you can show/hide chains and set styles separately per category.
- **New epitope:** Define a distance-based epitope from the 3D structure (see §5).
- **Other buttons:** Full screen, Export PNG, Close.

---

## 7. Other controls

### fasta (Download)
- **Function:** Download the **current view** (NT or AA) as FASTA. Filename includes mode, e.g. `alignment_AA.fasta`.

### Clear Alignment
- **Function:** Remove all sequences and reset state (groups, clusters, tree, founder designation, etc.). Returns to initial empty state.

### Help (?)
- **Function:** Toggle the help overlay with short descriptions of each control. For full documentation, open the published **aliViz user guide** (link in the help overlay).

---

## 8. Color legend

- **Groups:** Lists group labels and colors for sequence **names** (from Group).
- **Clusters:** Lists cluster IDs (numeric labels) and colors for tree tips (from any clustering method). Noise appears when applicable.
- **Note:** The floating legend does **not** include an alignment-length or phylogenetic scale bar; phylogenetic scale bars appear only on **exported SVG** tree figures (see §3).
- Cluster colors are removed when clustering is cleared (e.g. **None** or **Cancel**). Group colors are preserved after Cancel by rekeying group membership to the stripped (no `_cl-*`) names.

---

## 9. Summary of algorithms and formulas

| Feature        | Algorithm / formula |
|----------------|---------------------|
| Load sanitize  | A–Z, `-`, `*` kept; other chars → `X`; alert user. |
| Has subtype    | If off, no subtype index; toggling clears group/tree/cluster state. |
| Group          | Split name by delimiter; group = field value; unique IDs. |
| Sort           | REF, SubType (if any) fixed; others by group ID then name. |
| Prune          | Remove sequences in selected groups; renumber group IDs. |
| Founder consensus | Majority per column over selected group; `consensus_of_*` name. |
| Founder medoid | Minimize sum of padded Hamming distances to other group sequences; keep original name; `[Founder]` label. |
| NJ tree        | PhyloTools from alignment (pairwise distances → NJ). |
| FastTree       | bioWASM FastTree on alignment. |
| Load Newick    | Parse Newick; validate leaf names vs alignment; no auto-reroot on subtype. |
| Reroot         | Find founder leaf; reroot on edge to that leaf. |
| Ladderize      | By weight: sort children by leaf count. By depth: sort by max root-to-leaf depth in subtree. |
| Histogram      | Root-to-leaf distance = sum of branch lengths; bin and plot. |
| SVG scale      | Auto: fit to 520 px (linear) / R=260 (circular, ½ px per unit). Fixed: branch length per cm; overflow check. DSI: expected depth (DSI × MPD) → 260 px (linear) / 130 px (circular); no overflow abort (clip right border); dotted expected-depth line. |
| Tree-clade     | DFS; new cluster when edge length > threshold; min-leaves → noise. |
| Tree-cut       | Three depth cutoffs; BFS assign clusters; min-leaves → noise. |
| Cluster None   | Clear `leafClusters`; strip `_cl-*`; keep groups. |
| MDS            | B = −0.5·H·D²·H; eigendecomposition; coords = eigenvectors × √(eigenvalues). |
| UMAP           | External UMAP on distances → 2D. |
| DBSCAN         | eps-neighborhood; minPts; density-connected components; subtype isolated. |
| Calinski–Harabasz | B/(k−1) and W/(n−k) on tree (MRCA/root distances); divide by 2<sup>k</sup>. |
| Epitope CSV    | Rows: name, region1, region2, … (e.g. `start:end` or single position). |
| Distance epitope | Min distance intrinsic→extrinsic &lt; cutoff (Å); map residues to alignment. |

---

*End of aliViz User Guide.*
