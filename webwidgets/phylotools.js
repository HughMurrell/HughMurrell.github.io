/**
 * Phylogenetic Analysis Tools (JC69 + Neighbor Joining)
 */

const PhyloTools = {

    /**
     * 1) Calculate JC69 distance between two nucleotide sequences.
     * Ignores positions where either sequence has a gap or ambiguity.
     */
    calculateJC69: function(seq1, seq2) {
        if (seq1.length !== seq2.length) {
            throw new Error("Sequences must be of equal length (aligned).");
        }

        const validBases = new Set(['A', 'C', 'G', 'T', 'U', 'a', 'c', 'g', 't', 'u']);
        let differences = 0;
        let validSites = 0;

        for (let i = 0; i < seq1.length; i++) {
            const n1 = seq1[i];
            const n2 = seq2[i];

            // Check if both positions contain valid nucleotides
            if (validBases.has(n1) && validBases.has(n2)) {
                validSites++;
                if (n1.toUpperCase() !== n2.toUpperCase()) {
                    differences++;
                }
            }
        }

        if (validSites === 0) return 0; // Avoid division by zero

        const p = differences / validSites;

        // Jukes-Cantor Correction
        // Formula: d = -3/4 * ln(1 - 4/3 * p)
        // If p >= 0.75, the distance is mathematically undefined (saturation)
        if (p >= 0.75) {
            return Infinity; // Represents infinite evolutionary distance
        }

        return -0.75 * Math.log(1 - (4 / 3) * p);
    },

    /**
     * 2) Takes an array of aligned sequences and returns a 2D distance matrix.
     */
    computeDistanceMatrix: function(sequences) {
        const n = sequences.length;
        const matrix = Array.from({ length: n }, () => new Array(n).fill(0));

        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                const dist = this.calculateJC69(sequences[i], sequences[j]);
                matrix[i][j] = dist;
                matrix[j][i] = dist;
            }
        }
        return matrix;
    },

    /**
     * 3) Neighbor-Joining (NJ) Algorithm.
     * Takes a distance matrix and labels, returns a Newick string.
     */
    buildNeighborJoiningTree: function(distMatrix, labels) {
        // Deep copy matrix to avoid modifying the input
        let D = distMatrix.map(row => [...row]);
        // Create a list of active nodes (starting with leaf labels)
        // Format: { name: "Label", id: original_index }
        let clusters = labels.map((label, i) => ({
            newick: label,
            id: i
        }));

        while (clusters.length > 2) {
            const N = clusters.length;

            // 1. Calculate Net Divergence (r) for each node
            const R = new Array(N).fill(0);
            for (let i = 0; i < N; i++) {
                for (let j = 0; j < N; j++) {
                    if (i !== j) R[i] += D[i][j];
                }
            }

            // 2. Calculate Q-matrix and find pair with minimum Q
            let minQ = Infinity;
            let pair = [-1, -1]; // indices in the current D/clusters arrays

            for (let i = 0; i < N; i++) {
                for (let j = i + 1; j < N; j++) {
                    // Q_ij = (N-2)d_ij - r_i - r_j
                    const qVal = (N - 2) * D[i][j] - R[i] - R[j];
                    if (qVal < minQ) {
                        minQ = qVal;
                        pair = [i, j];
                    }
                }
            }

            const [i, j] = pair;

            // 3. Calculate branch lengths from existing nodes i and j to new node u
            // d_iu = 0.5 * d_ij + (1 / (2(N-2))) * (r_i - r_j)
            const dist_ij = D[i][j];
            const val_i_u = 0.5 * dist_ij + (1 / (2 * (N - 2))) * (R[i] - R[j]);
            const val_j_u = dist_ij - val_i_u;

            // 4. Create new Newick string
            // Note: NJ can sometimes produce negative branch lengths due to stochastic noise;
            // standard implementations usually report them as-is.
            const newNodeName = `(${clusters[i].newick}:${val_i_u.toFixed(5)},${clusters[j].newick}:${val_j_u.toFixed(5)})`;

            // 5. Update Distance Matrix
            // Calculate distances from new node u to all other nodes k
            // d_uk = 0.5 * (d_ik + d_jk - d_ij)
            const newDistRow = [];
            for (let k = 0; k < N; k++) {
                if (k !== i && k !== j) {
                    const d_uk = 0.5 * (D[i][k] + D[j][k] - dist_ij);
                    newDistRow.push(d_uk);
                }
            }

            // 6. Manipulate arrays to remove i and j and add u
            // Remove the higher index first to maintain indices of lower items
            const removeFirst = Math.max(i, j);
            const removeSecond = Math.min(i, j);

            // Update Clusters
            // We keep "others" then add "new"
            const nextClusters = clusters.filter((_, idx) => idx !== i && idx !== j);
            nextClusters.push({ newick: newNodeName });
            clusters = nextClusters;

            // Update Matrix D
            // Remove rows
            let nextD = D.filter((_, idx) => idx !== i && idx !== j);
            // Remove columns from remaining rows
            nextD = nextD.map(row => row.filter((_, idx) => idx !== i && idx !== j));
            
            // Add new column (distances to u) to existing rows
            for (let k = 0; k < nextD.length; k++) {
                nextD[k].push(newDistRow[k]);
            }
            // Add new row (distances to u) + 0 for self-distance
            newDistRow.push(0);
            nextD.push(newDistRow);
            
            D = nextD;
        }

        // Final step: Connect the last two remaining clusters
        // For the final two, the branch length is simply the distance between them.
        // It is customary to root it arbitrarily at the midpoint or just split the distance.
        // Here we format as (A:dist, B:dist) effectively treating the connection as a root.
        const distFinal = D[0][1];
        // Standard unrooted newick often looks like (A:len, B:len, C:len);
        // But strictly bifurcating at the end:
        return `(${clusters[0].newick}:${(distFinal/2).toFixed(5)},${clusters[1].newick}:${(distFinal/2).toFixed(5)});`;
    },

    /**
     * 4) Convenience Wrapper
     */
    generatePhylogeny: function(sequences, labels) {
        if (!sequences || !labels || sequences.length !== labels.length) {
            throw new Error("Mismatch between sequences and labels.");
        }
        
        console.log("Calculating distance matrix...");
        const distMatrix = this.computeDistanceMatrix(sequences);
        
        console.log("Building Neighbor-Joining tree...");
        const newick = this.buildNeighborJoiningTree(distMatrix, labels);
        
        return newick;
    }
};

// --- Example Usage ---
// const labels = ["SeqA", "SeqB", "SeqC", "SeqD"];
// const seqs = [
//   "ATGCATGCAT", // Clean
//   "ATGCATGC-T", // Has gap (ignored for dist)
//   "ATCCATGCAT", // 1 diff from A
//   "AGGCATGCAT"  // 1 diff from A
// ];
// console.log(PhyloTools.generatePhylogeny(seqs, labels));