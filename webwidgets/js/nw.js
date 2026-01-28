/**
 * Needleman-Wunsch Alignment Library
 * Hosted at: https://murrellgroup.github.io/WebWidgets/nw.js
 */

function affineNWAlign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const s1arr = s1.split('');
    const s1len = s1arr.length;
    const s2arr = s2.split('');
    const s2len = s2arr.length;
    const M = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const IX = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const IY = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    for (let i = 0; i <= s1len; i++) {
        IX[i][0] = gapOpen + gapExtend * i;
        IY[i][0] = -Infinity;
        M[i][0] = gapOpen + gapExtend * i;
    }
    for (let j = 0; j <= s2len; j++) {
        IX[0][j] = -Infinity;
        IY[0][j] = gapOpen + gapExtend * j;
        M[0][j] = gapOpen + gapExtend * j;
    }
    const traceM = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const traceIX = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const traceIY = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    traceM[0].fill(3);
    traceIX[0].fill(3);
    traceIY[0].fill(3);
    for (let i = 1; i <= s1len; i++) {
        traceM[i][0] = 2;
        traceIX[i][0] = 2;
        traceIY[i][0] = 2;
    }
    for (let i = 1; i <= s1len; i++) {
        for (let j = 1; j <= s2len; j++) {
            const diagCost = s1arr[i - 1] === s2arr[j - 1] ? matchCost : mismatchCost;
            const diagM = M[i - 1][j - 1] + diagCost;
            const IX2M = IX[i - 1][j - 1] + diagCost;
            const IY2M = IY[i - 1][j - 1] + diagCost;
            [M[i][j], traceM[i][j]] = findMax([diagM, IX2M, IY2M]);

            const boundaryGapExtendX = (i === s1len) ? gapExtend / boundaryGapFactor : gapExtend;
            const boundaryGapExtendY = (j === s2len) ? gapExtend / boundaryGapFactor : gapExtend;
            
            const M2IX = M[i - 1][j] + gapOpen;
            const IXextend = IX[i - 1][j] + boundaryGapExtendY;
            [IX[i][j], traceIX[i][j]] = findMax([M2IX, IXextend]);

            const M2IY = M[i][j - 1] + gapOpen;
            const IYextend = IY[i][j - 1] + boundaryGapExtendX;
            [IY[i][j], traceIY[i][j]] = findMax([M2IY, -Infinity, IYextend]);
        }
    }
    const revArr1 = [];
    const revArr2 = [];
    const mats = [traceM, traceIX, traceIY];
    let xI = s1len;
    let yI = s2len;
    let mI = findMax([M[xI][yI], IX[xI][yI], IY[xI][yI]])[1];
    while (xI > 0 && yI > 0) {
        const nextMI = mats[mI - 1][xI][yI];
        if (mI === 1) {
            revArr1.push(s1arr[xI - 1]);
            revArr2.push(s2arr[yI - 1]);
            xI--;
            yI--;
        } else if (mI === 2) {
            revArr1.push(s1arr[xI - 1]);
            revArr2.push('-');
            xI--;
        } else if (mI === 3) {
            revArr1.push('-');
            revArr2.push(s2arr[yI - 1]);
            yI--;
        }
        mI = nextMI;
    }
    while (xI > 0) {
        revArr1.push(s1arr[xI - 1]);
        revArr2.push('-');
        xI--;
    }
    while (yI > 0) {
        revArr1.push('-');
        revArr2.push(s2arr[yI - 1]);
        yI--;
    }
    return [revArr1.reverse().join(''), revArr2.reverse().join('')];
}

function findMax(arr) {
    let maxVal = arr[0];
    let maxIndex = 1;
    for (let i = 1; i < arr.length; i++) {
        if (arr[i] > maxVal) {
            maxVal = arr[i];
            maxIndex = i + 1;
        }
    }
    return [maxVal, maxIndex];
}









/**
 * Constrained NW Aligner
 * Used for aligning gaps between fixed blocks.
 * 
 * isStart: true if this segment is the very beginning of the sequences (Left Flank).
 * isEnd:   true if this segment is the very end of the sequences (Right Flank).
 * 
 * If isEnd is false, boundaryGapFactor is ignored (treated as 1.0) to prevent free gaps 
 * where we should be connecting to a match.
 * M[0][0] is explicitly 0 to ensure we start from a "Match" state unless isStart implies otherwise.
 */
function constrainedNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor, isStart, isEnd) {
    const s1arr = s1.split('');
    const s1len = s1arr.length;
    const s2arr = s2.split('');
    const s2len = s2arr.length;

    // Use boundary factor ONLY if we are at the true end of the alignment
    const useBoundaryFactor = isEnd ? boundaryGapFactor : 1;

    const M = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const IX = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const IY = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));

    // Initialization
    // M[0][0] = 0 implies we are in a match state (or start of seq)
    M[0][0] = 0; 
    IX[0][0] = -Infinity;
    IY[0][0] = -Infinity;

    for (let i = 1; i <= s1len; i++) {
        IX[i][0] = gapOpen + gapExtend * i;
        IY[i][0] = -Infinity;
        M[i][0] = gapOpen + gapExtend * i;
    }
    for (let j = 1; j <= s2len; j++) {
        IX[0][j] = -Infinity;
        IY[0][j] = gapOpen + gapExtend * j;
        M[0][j] = gapOpen + gapExtend * j;
    }

    const traceM = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const traceIX = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const traceIY = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    
    // Traceback initialization
    traceM[0].fill(3);
    traceIX[0].fill(3);
    traceIY[0].fill(3);
    for (let i = 1; i <= s1len; i++) {
        traceM[i][0] = 2;
        traceIX[i][0] = 2;
        traceIY[i][0] = 2;
    }

    for (let i = 1; i <= s1len; i++) {
        for (let j = 1; j <= s2len; j++) {
            const diagCost = s1arr[i - 1] === s2arr[j - 1] ? matchCost : mismatchCost;
            const diagM = M[i - 1][j - 1] + diagCost;
            const IX2M = IX[i - 1][j - 1] + diagCost;
            const IY2M = IY[i - 1][j - 1] + diagCost;
            [M[i][j], traceM[i][j]] = findMax([diagM, IX2M, IY2M]);

            // Apply factor only if isEnd is true (via useBoundaryFactor)
            const boundaryGapExtendX = (i === s1len) ? gapExtend / useBoundaryFactor : gapExtend;
            const boundaryGapExtendY = (j === s2len) ? gapExtend / useBoundaryFactor : gapExtend;
            
            const M2IX = M[i - 1][j] + gapOpen;
            const IXextend = IX[i - 1][j] + boundaryGapExtendY;
            [IX[i][j], traceIX[i][j]] = findMax([M2IX, IXextend]);

            const M2IY = M[i][j - 1] + gapOpen;
            const IYextend = IY[i][j - 1] + boundaryGapExtendX;
            [IY[i][j], traceIY[i][j]] = findMax([M2IY, -Infinity, IYextend]);
        }
    }

    const revArr1 = [];
    const revArr2 = [];
    const mats = [traceM, traceIX, traceIY];
    let xI = s1len;
    let yI = s2len;
    let mI = findMax([M[xI][yI], IX[xI][yI], IY[xI][yI]])[1];

    while (xI > 0 && yI > 0) {
        const nextMI = mats[mI - 1][xI][yI];
        if (mI === 1) {
            revArr1.push(s1arr[xI - 1]);
            revArr2.push(s2arr[yI - 1]);
            xI--;
            yI--;
        } else if (mI === 2) {
            revArr1.push(s1arr[xI - 1]);
            revArr2.push('-');
            xI--;
        } else if (mI === 3) {
            revArr1.push('-');
            revArr2.push(s2arr[yI - 1]);
            yI--;
        }
        mI = nextMI;
    }
    while (xI > 0) {
        revArr1.push(s1arr[xI - 1]);
        revArr2.push('-');
        xI--;
    }
    while (yI > 0) {
        revArr1.push('-');
        revArr2.push(s2arr[yI - 1]);
        yI--;
    }
    return [revArr1.reverse().join(''), revArr2.reverse().join('')];
}

/**
 * Task 1: K-mer Seeded Alignment
 * Uses LIS to find consistent anchors and Erodes edges to prevent boundary issues.
 */
function kmer_seeded_nwalign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const K = 25; 
    const minLen = Math.min(s1.length, s2.length);
    
    // Fallback if sequences are too short to support seeds + erosion
    if (minLen < 3 * K) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 1. Identify Unique K-mers
    const getUniqueKmers = (str) => {
        const counts = new Map();
        const indices = new Map();
        for (let i = 0; i <= str.length - K; i++) {
            const sub = str.substring(i, i + K);
            counts.set(sub, (counts.get(sub) || 0) + 1);
            indices.set(sub, i);
        }
        return { counts, indices };
    };

    const map1 = getUniqueKmers(s1);
    const map2 = getUniqueKmers(s2);

    let matches = [];
    for (const [kmer, count] of map1.counts) {
        if (count === 1 && map2.counts.get(kmer) === 1) {
            matches.push({
                i: map1.indices.get(kmer),
                j: map2.indices.get(kmer),
                len: K
            });
        }
    }

    if (matches.length === 0) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 2. LIS for consistency (Sort by I, LIS on J)
    matches.sort((a, b) => a.i - b.i);
    const lisMatches = getLIS(matches);
    
    // 3. Merge and Erode to create safe anchors
    const anchors = mergeAndErodeAnchors(lisMatches, K);
    
    // 4. Stitch using constrained NW
    return stitchAlignments(s1, s2, anchors, [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]);
}

/**
 * Task 2: Double DP Alignment
 * Uses Sparse DP to chain small seeds, then erodes edges.
 */
function doubleDP_nwalign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const K = 11;
    const minLen = Math.min(s1.length, s2.length);

    if (minLen < 3 * K) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 1. Index s1
    const s1Map = new Map();
    for (let i = 0; i <= s1.length - K; i++) {
        const sub = s1.substring(i, i + K);
        if (!s1Map.has(sub)) s1Map.set(sub, []);
        s1Map.get(sub).push(i);
    }

    // 2. Find Seeds in s2
    let matches = [];
    const maxHits = 25; 
    
    for (let j = 0; j <= s2.length - K; j++) {
        const sub = s2.substring(j, j + K);
        const hits = s1Map.get(sub);
        if (hits && hits.length <= maxHits) {
            for (let i of hits) {
                matches.push({ i, j, len: K });
            }
        }
    }

    if (matches.length === 0) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    matches.sort((a, b) => (a.i - b.i) || (a.j - b.j));

    // 3. Sparse DP Chain
    const scores = new Float64Array(matches.length);
    const parents = new Int32Array(matches.length).fill(-1);
    const lookback = 80;

    for (let cur = 0; cur < matches.length; cur++) {
        const mCurr = matches[cur];
        let maxScore = mCurr.len; 
        
        const startSearch = Math.max(0, cur - lookback);
        for (let prev = cur - 1; prev >= startSearch; prev--) {
            const mPrev = matches[prev];

            // Strict time-forward check
            if (mPrev.i < mCurr.i && mPrev.j < mCurr.j) {
                
                const diagDiff = Math.abs((mCurr.j - mCurr.i) - (mPrev.j - mPrev.i));
                const gapI = mCurr.i - (mPrev.i + mPrev.len);
                const gapJ = mCurr.j - (mPrev.j + mPrev.len);
                const dist = Math.max(0, gapI) + Math.max(0, gapJ);
                
                // Heuristic scoring
                const penalty = (diagDiff * 3.0) + (dist * 0.1); 
                
                const newScore = scores[prev] + mCurr.len - penalty;
                if (newScore > maxScore) {
                    maxScore = newScore;
                    parents[cur] = prev;
                }
            }
        }
        scores[cur] = maxScore;
    }

    let bestIdx = 0;
    let maxS = scores[0];
    for(let i=1; i<matches.length; i++){
        if(scores[i] > maxS){
            maxS = scores[i];
            bestIdx = i;
        }
    }

    const chain = [];
    let curr = bestIdx;
    while(curr !== -1) {
        chain.push(matches[curr]);
        curr = parents[curr];
    }
    chain.reverse();

    // 4. Merge and Erode
    const anchors = mergeAndErodeAnchors(chain, K);

    return stitchAlignments(s1, s2, anchors, [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]);
}

// --- Helpers ---

function getLIS(matches) {
    if (matches.length === 0) return [];
    const tails = []; 
    const parent = new Int32Array(matches.length).fill(-1);
    for (let i = 0; i < matches.length; i++) {
        const val = matches[i].j;
        let left = 0, right = tails.length;
        while (left < right) {
            const mid = (left + right) >>> 1;
            if (matches[tails[mid]].j < val) left = mid + 1;
            else right = mid;
        }
        if (left < tails.length) {
            tails[left] = i;
            parent[i] = (left > 0) ? tails[left - 1] : -1;
        } else {
            tails.push(i);
            parent[i] = (tails.length > 1) ? tails[tails.length - 2] : -1;
        }
    }
    const result = [];
    let curr = tails[tails.length - 1];
    while (curr !== -1) {
        result.push(matches[curr]);
        curr = parent[curr];
    }
    return result.reverse();
}

/**
 * Merges continuous/overlapping blocks on same diagonal, 
 * then erodes ends by 'margin' to ensure safe stitching.
 */
function mergeAndErodeAnchors(matches, margin) {
    if (matches.length === 0) return [];
    
    // Phase 1: Merge
    const merged = [];
    let curr = { ...matches[0] };

    for (let k = 1; k < matches.length; k++) {
        const next = matches[k];
        
        const diagCurr = curr.j - curr.i;
        const diagNext = next.j - next.i;
        
        // Check connectivity/overlap
        const isConnected = next.i <= (curr.i + curr.len);

        if (diagCurr === diagNext && isConnected) {
            const newEnd = Math.max(curr.i + curr.len, next.i + next.len);
            curr.len = newEnd - curr.i;
        } else {
            merged.push(curr);
            curr = { ...next };
        }
    }
    merged.push(curr);

    // Phase 2: Erode
    const eroded = [];
    for (const m of merged) {
        // "Chew back" margin from both sides
        const newLen = m.len - (2 * margin);
        
        if (newLen > 0) {
            eroded.push({
                i: m.i + margin,
                j: m.j + margin,
                len: newLen
            });
        }
    }
    return eroded;
}

function stitchAlignments(s1, s2, anchors, nwParams) {
    let finalS1 = "";
    let finalS2 = "";
    let idx1 = 0;
    let idx2 = 0;

    // nwParams: [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]

    for (let k = 0; k < anchors.length; k++) {
        const anchor = anchors[k];

        // --- GAP REGION ---
        const sub1 = s1.substring(idx1, anchor.i);
        const sub2 = s2.substring(idx2, anchor.j);
        
        // Check if we are at Start or End of entire sequence
        // Left Flank (Start): k === 0 && idx1 === 0
        const isStart = (k === 0 && idx1 === 0);
        // Internal Gap: isStart = false, isEnd = false (since we are connecting to an anchor)
        const isEnd = false;

        // Perform constrained alignment on the gap
        // Even if strings are empty, run it if we need to ensure state transition logic, 
        // though typically empty aligns to empty.
        if (sub1.length > 0 || sub2.length > 0) {
            const [aln1, aln2] = constrainedNWAlign(sub1, sub2, ...nwParams, isStart, isEnd);
            finalS1 += aln1;
            finalS2 += aln2;
        }

        // --- ANCHOR REGION (Perfect Match) ---
        const anchorSeq = s1.substring(anchor.i, anchor.i + anchor.len);
        finalS1 += anchorSeq;
        finalS2 += anchorSeq;

        idx1 = anchor.i + anchor.len;
        idx2 = anchor.j + anchor.len;
    }

    // --- TAIL REGION (Right Flank) ---
    if (idx1 < s1.length || idx2 < s2.length) {
        const sub1 = s1.substring(idx1);
        const sub2 = s2.substring(idx2);
        
        // This is the Right Flank: isStart=false, isEnd=true
        const [aln1, aln2] = constrainedNWAlign(sub1, sub2, ...nwParams, false, true);
        finalS1 += aln1;
        finalS2 += aln2;
    }

    return [finalS1, finalS2];
}

/**
 * Partial Order Alignment (POA) Library
 * Extends the NW logic to align sequences against a growing graph.
 */

// --- Data Structures ---

class POANode {
    constructor(id, char) {
        this.id = id;
        this.char = char;
        this.next = []; // Array of POANode IDs
        this.prev = []; // Array of POANode IDs
        this.seqs = new Set(); // IDs of sequences visiting this node
        this.alignedTo = null; // Used during alignment (traceback ptr)
    }
}

class POAGraph {
    constructor() {
        this.nodes = []; // Storage: index is ID
        this.sequences = []; // Keep track of original sequence strings
        this.nodeCount = 0;
    }

    createNode(char) {
        const node = new POANode(this.nodeCount++, char);
        this.nodes.push(node);
        return node;
    }

    addEdge(fromId, toId) {
        const from = this.nodes[fromId];
        const to = this.nodes[toId];
        if (!from.next.includes(toId)) from.next.push(toId);
        if (!to.prev.includes(fromId)) to.prev.push(fromId);
    }

    // Initialize graph with the first sequence (linear chain)
    initFirstSequence(seq) {
        this.sequences.push(seq);
        let prevNode = null;
        for (let i = 0; i < seq.length; i++) {
            const node = this.createNode(seq[i]);
            node.seqs.add(0);
            if (prevNode) {
                this.addEdge(prevNode.id, node.id);
            }
            prevNode = node;
        }
    }

    /**
     * Topological Sort (Kahn's Algorithm usually, but DFS post-order is fine for DAG)
     * Returns array of Node IDs.
     */
    getTopologicalSort() {
        const visited = new Array(this.nodeCount).fill(false);
        const stack = [];
        
        const visit = (u) => {
            visited[u] = true;
            for (const vId of this.nodes[u].next) {
                if (!visited[vId]) visit(vId);
            }
            stack.push(u);
        };

        // Handle disconnected components if any (though unlikely in this use case)
        for (let i = 0; i < this.nodeCount; i++) {
            if (!visited[i]) visit(i);
        }
        
        return stack.reverse(); // Reverse post-order is topological order
    }

    /**
     * Extracts the "Heaviest Bundle" path.
     * This serves as a temporary "consensus" sequence for K-mer seeding.
     * Returns: { str: string, mapping: Array<NodeID> }
     */
    getConsensusPath() {
        const sorted = this.getTopologicalSort();
        const scores = new Array(this.nodeCount).fill(0);
        const predecessors = new Array(this.nodeCount).fill(-1);

        // Simple weight: number of sequences passing through node
        for (const u of sorted) {
            const weight = this.nodes[u].seqs.size;
            const currentScore = (scores[u] || 0) + weight;
            scores[u] = currentScore;

            for (const v of this.nodes[u].next) {
                if (scores[u] > scores[v]) {
                    scores[v] = scores[u];
                    predecessors[v] = u;
                }
            }
        }

        // Find max score end node
        let maxScore = -1;
        let curr = -1;
        for (let i = 0; i < this.nodeCount; i++) {
            if (this.nodes[i].next.length === 0 && scores[i] > maxScore) {
                maxScore = scores[i];
                curr = i;
            }
        }

        // Backtrack
        const path = [];
        while (curr !== -1) {
            path.push(curr);
            curr = predecessors[curr]; // This logic is slightly flawed for DAG traversal max weight, 
                                       // but standard topological DP works for finding heaviest path.
        }
        
        // Fix: The loop above finds path to specific end. 
        // Better: Standard DP. Score[v] = weight(v) + max(Score[u] for u in prev).
        // Let's redo quickly for correctness.
        const pathScore = new Array(this.nodeCount).fill(0);
        const bestPrev = new Array(this.nodeCount).fill(-1);
        
        for (const u of sorted) {
            const w = this.nodes[u].seqs.size;
            let maxPrevScore = 0;
            let bp = -1;
            
            for (const p of this.nodes[u].prev) {
                if (pathScore[p] > maxPrevScore) {
                    maxPrevScore = pathScore[p];
                    bp = p;
                }
            }
            pathScore[u] = w + maxPrevScore;
            bestPrev[u] = bp;
        }

        // Find global max
        let bestEnd = 0; 
        for(let i=0; i<this.nodeCount; i++) {
            if(pathScore[i] > pathScore[bestEnd]) bestEnd = i;
        }

        const resIds = [];
        let ptr = bestEnd;
        while(ptr !== -1) {
            resIds.push(ptr);
            ptr = bestPrev[ptr];
        }
        resIds.reverse();

        return {
            str: resIds.map(id => this.nodes[id].char).join(''),
            mapping: resIds
        };
    }
}

// --- Sequence to Graph Alignment ---

/**
 * Aligns a new sequence to the POAGraph.
 * 1. Generates "Consensus" string from Graph.
 * 2. Uses K-mer seeding to align NewSeq to Consensus.
 * 3. Maps seeds back to Graph Node IDs.
 * 4. Runs constrained Sequence-to-Graph DP between seeds.
 * 5. Mutates graph (fusing nodes or adding bubbles).
 */
function addSequenceToGraph(graph, sequence, nwParams) {
    const seqIdx = graph.sequences.length;
    graph.sequences.push(sequence);
    
    // 1. Get Consensus (Graph approximation)
    const consensus = graph.getConsensusPath(); // { str, mapping: [nodeId, ...] }
    
    // 2. Find Seeds (Seq vs Consensus)
    // We reuse the logic from the pairwise K-mer aligner, but we just want the anchors.
    // We assume getAnchors is a helper derived from kmer_seeded_nwalign logic.
    const anchors = findSeqToGraphAnchors(sequence, consensus.str);

    // 3. Convert anchors (seqIdx, consIdx) -> (seqIdx, graphNodeId)
    const constraints = anchors.map(a => ({
        seqStart: a.i,
        seqEnd: a.i + a.len,
        nodeStart: consensus.mapping[a.j],
        nodeEnd: consensus.mapping[a.j + a.len - 1]
    }));

    // 4. Perform Spliced Alignment
    // We iterate through gaps between constraints and the constraints themselves.
    let currSeqPos = 0;
    
    // Helper to track the "previous" node ID we attached to in the graph
    // Start with a virtual "source" if needed, but we handle first align manually.
    let prevNodeId = -1; 

    // Sort graph topologically for DP
    const topoOrder = graph.getTopologicalSort(); 
    const topoIndexMap = new Int32Array(graph.nodeCount).fill(-1);
    for(let i=0; i<topoOrder.length; i++) topoIndexMap[topoOrder[i]] = i;

    // Iterate through constraints
    for (let k = 0; k <= constraints.length; k++) {
        let blockSeqEnd, blockNodeTarget;
        let isAnchor = false;

        if (k < constraints.length) {
            // Gap before anchor
            blockSeqEnd = constraints[k].seqStart;
            blockNodeTarget = constraints[k].nodeStart; 
        } else {
            // Final tail
            blockSeqEnd = sequence.length;
            blockNodeTarget = null; // No specific target, just align to end of graph
        }

        // --- Align Segment (Gap) ---
        // We align sequence[currSeqPos ... blockSeqEnd]
        // to the graph region between prevNodeId and blockNodeTarget.
        // This is the hard part: Sequence-to-Graph DP.
        
        if (currSeqPos < blockSeqEnd || (prevNodeId !== -1 && blockNodeTarget !== null)) {
            const subSeq = sequence.substring(currSeqPos, blockSeqEnd);
            
            // Subgraph selection: Logic to define relevant subgraph is complex.
            // Simplification: We align subSeq to the ENTIRE graph but constrained 
            // by the fact that we must start after prevNodeId and end before blockNodeTarget.
            // Performance trick: Extract subgraph between prevNodeId and blockNodeTarget.
            
            const subgraphNodes = getSubgraph(graph, prevNodeId, blockNodeTarget, topoOrder, topoIndexMap);
            
            if (subSeq.length > 0 || subgraphNodes.length > 0) {
                const alignPath = runSeqToGraphDP(subSeq, subgraphNodes, graph, nwParams);
                
                // Mutate Graph based on alignment
                prevNodeId = applyAlignmentToGraph(graph, alignPath, subSeq, prevNodeId, seqIdx);
            }
        }

        // --- Align Anchor (Match) ---
        if (k < constraints.length) {
            const constr = constraints[k];
            const len = constr.seqEnd - constr.seqStart;
            
            // The anchor implies a direct walk along the consensus path
            // We just fuse the sequence to these nodes.
            // Note: The consensus path nodes are guaranteed to be connected linearly? 
            // Usually yes, but let's be safe.
            
            // Find path of nodes in consensus mapping
            const startIdxInCons = consensus.mapping.indexOf(constr.nodeStart);
            
            for(let i=0; i<len; i++) {
                const char = sequence[constr.seqStart + i];
                const targetNodeId = consensus.mapping[startIdxInCons + i];
                
                // Sanity check char matches (seeds should match)
                // If mismatch (due to seed hash collision or whatever), we force it or bubble.
                // Assuming seed is perfect match:
                const node = graph.nodes[targetNodeId];
                node.seqs.add(seqIdx);
                
                if (prevNodeId !== -1 && prevNodeId !== targetNodeId) {
                    graph.addEdge(prevNodeId, targetNodeId);
                }
                prevNodeId = targetNodeId;
            }
            currSeqPos = constr.seqEnd;
        }
    }
}

/**
 * Standard K-mer finder reused for Seq-to-Consensus
 */
function findSeqToGraphAnchors(s1, s2) {
    const K = 15; // Seed size
    // Use the logic from `kmer_seeded_nwalign` roughly
    if (s1.length < K || s2.length < K) return [];

    const getUniqueKmers = (str) => {
        const counts = new Map();
        const indices = new Map();
        for (let i = 0; i <= str.length - K; i++) {
            const sub = str.substring(i, i + K);
            counts.set(sub, (counts.get(sub) || 0) + 1);
            indices.set(sub, i);
        }
        return { counts, indices };
    };

    const map1 = getUniqueKmers(s1);
    const map2 = getUniqueKmers(s2);

    let matches = [];
    for (const [kmer, count] of map1.counts) {
        if (count === 1 && map2.counts.get(kmer) === 1) { // Unique seeds only
            matches.push({
                i: map1.indices.get(kmer),
                j: map2.indices.get(kmer),
                len: K
            });
        }
    }
    
    // Sort and LIS
    matches.sort((a, b) => a.i - b.i);
    // Reuse LIS from context if available, otherwise simple implementation
    const lis = simpleLIS(matches);
    
    // Check consistency (must be increasing in both i and j)
    return mergeAndErodeAnchors(lis, 2); // Erode slightly to allow local flexibility
}

function simpleLIS(matches) {
    if(matches.length === 0) return [];
    // Standard O(N log N) LIS on 'j' coordinate
    const tails = []; 
    const parent = new Int32Array(matches.length).fill(-1);
    for (let i = 0; i < matches.length; i++) {
        const val = matches[i].j;
        let left = 0, right = tails.length;
        while (left < right) {
            const mid = (left + right) >>> 1;
            if (matches[tails[mid]].j < val) left = mid + 1;
            else right = mid;
        }
        if (left < tails.length) {
            tails[left] = i;
            parent[i] = (left > 0) ? tails[left - 1] : -1;
        } else {
            tails.push(i);
            parent[i] = (tails.length > 1) ? tails[tails.length - 2] : -1;
        }
    }
    const result = [];
    let curr = tails[tails.length - 1];
    while (curr !== -1) {
        result.push(matches[curr]);
        curr = parent[curr];
    }
    return result.reverse();
}

/**
 * Extract a subgraph of nodes topologically between Start and End.
 * This is an optimization. If start is -1, include from beginning.
 */
function getSubgraph(graph, startId, endId, topoOrder, topoIndexMap) {
    const nodes = [];
    const startIndex = (startId === -1) ? -1 : topoIndexMap[startId];
    const endIndex = (endId === null) ? Infinity : topoIndexMap[endId];

    for (const u of topoOrder) {
        const idx = topoIndexMap[u];
        if (idx > startIndex && idx < endIndex) {
            nodes.push(u);
        }
    }
    return nodes;
}

/**
 * Core Logic: Sequence to Graph DP
 * Returns a list of operations: ['MATCH', nodeId], ['INSERT', char], ['DELETE', nodeId]
 */
function runSeqToGraphDP(seq, graphNodes, fullGraph, [gapOpen, gapExtend, matchCost, mismatchCost]) {
    const seqLen = seq.length;
    const subNodeCount = graphNodes.length;
    
    // Map subgraph node index back to real ID
    const nodeMap = graphNodes; 
    // Inverse map for quick lookup of "index in subgraph" from "real ID"
    // Since subNodeCount might be small relative to full graph, use a Map or sparse array.
    const realIdToSubIdx = new Map();
    graphNodes.forEach((id, idx) => realIdToSubIdx.set(id, idx));

    // DP Matrices
    // M[i][j]: Best score aligning seq[0..i] to subgraph node j
    // We strictly use topological order for 'j'.
    const M = Array(seqLen + 1).fill().map(() => Array(subNodeCount).fill(-Infinity));
    
    // Initialization
    // We assume we are entering the subgraph from 'prevNodeId' (handled outside)
    // So row 0 (empty seq) aligns to... 
    // This part is tricky in sub-graph DP. We assume "Free Entry" into the first layer 
    // of the subgraph relative to the implicit previous node, but penalties apply.
    // For simplicity: We treat this as a standard global alignment within the box defined by anchors.
    
    // Initialize Row 0 (Deletions in Seq / Moves in Graph)
    // Score depends on incoming edges from outside the subgraph or within.
    // Simplifying assumption: Standard Affine Gap initialization
    // M[0][j] = gapOpen + gapExtend * (distance from start of subgraph) -> Heuristic
    
    // Since we are constrained, we can set M[0][j] = -Inf mostly, except strict first nodes?
    // Let's use a standard trick: Add a virtual Start Node and End Node?
    // Or just fill Row 0 based on "Delete" penalty accumulation.
    
    // Let's iterate.
    // Init local scores
    const SCORES = Array(subNodeCount).fill(-Infinity); // Best score ending at node j with full alignment so far

    // Need "Deletion" state (Match in graph, gap in seq) and "Insertion" state (Match in seq, gap in graph)
    // To keep it single-function simple, we stick to a simpler gap-affine model with 3 states per cell is best,
    // but 1 state + checking neighbors is easier to code for graph.
    
    // Let's use the layout: M[seqIndex][nodeIndex].
    // Value is tuple: [score, traceType, tracePrevNodeIndex]
    
    // Initialize
    // i=0 corresponds to "before first char of seq"
    
    // Need to handle predecessors.
    // Preds of node u in subgraph might include nodes OUTSIDE subgraph (the anchor).
    // We treat 'anchor' as having score 0 at seq index -1.
    
    // Let's go simple:
    // M[i][u] = max(
    //   M[i-1][u] + gap_cost (Insertion: consume seq, stay at u? No, create new node),
    //   max_v(M[i-1][v]) + match_cost(seq[i], u) (Match/Mismatch),
    //   max_v(M[i][v]) + gap_cost (Deletion: skip u? No, move to u without consuming seq)
    // )

    // Correct POA Recurrence:
    // D[i][u] = Score aligning seq[1..i] ending at node u.
    // D[i][u] = max(
    //    max_{v in pred(u)} D[i-1][v] + score(seq[i], u), // Match/Mismatch
    //    D[i-1][u] + gap_extend (Insertion - wait, this assumes u is linear. In POA insertion is a new branch),
    //    max_{v in pred(u)} D[i][v] + gap_extend // Deletion
    // )
    
    // Actually, Insertion in POA means aligning seq[i] to a GAP in the graph.
    // This is usually modeled by D[i][u] (Match), I[i][u] (Ins), E[i][u] (Del).
    // But since we are building the graph, "Insertion" means we will create a NEW node.
    // During alignment, we model Insertion as aligning to a virtual state.
    
    // Let's use a simplified heuristic for "approximate" MSA as requested:
    // We treat the subgraph nodes as candidates. 
    // If we insert, we come from M[i-1][v] but don't advance graph index? 
    // No, standard NW:
    // Match: (i-1, v) -> (i, u)
    // Del: (i, v) -> (i, u)  (consuming graph node u without seq char)
    // Ins: (i-1, u) -> (i, u) (consuming seq char at node u? No, "floating" after u)

    const scoreMat = new Float32Array((seqLen + 1) * subNodeCount).fill(-1e9);
    const traceMat = new Int8Array((seqLen + 1) * subNodeCount).fill(0); // 1: Diag, 2: Ins, 3: Del
    const tracePred = new Int32Array((seqLen + 1) * subNodeCount).fill(-1); // Stores index of v for Diag/Del

    // Helper to access matrix
    const idx = (i, j) => i * subNodeCount + j;

    // Initialization: i=0 (Seq empty)
    // We can travel through the graph (Deletions)
    // If a node u is connected to start anchor, score is gapOpen + gapExtend.
    // If connected to another node v in subgraph, score is score[v] + gapExtend.
    
    for (let u = 0; u < subNodeCount; u++) {
        const uId = graphNodes[u];
        const preds = fullGraph.nodes[uId].prev;
        
        let maxPrev = -1e9;
        let bestP = -1;

        // Check connection to Start Anchor (outside subgraph)
        // If predecessor is NOT in graphNodes, it's the anchor (or effectively so).
        const isStartConnected = preds.some(p => !realIdToSubIdx.has(p));
        
        if (isStartConnected) {
            // Connected to anchor
            maxPrev = gapOpen; // First gap
        }

        // Check internal predecessors (Topo order ensures they are computed)
        for (const pId of preds) {
            if (realIdToSubIdx.has(pId)) {
                const pIdx = realIdToSubIdx.get(pId);
                const s = scoreMat[idx(0, pIdx)];
                if (s > -1e9) {
                    const val = s + gapExtend; // Extend gap
                    if (val > maxPrev) {
                        maxPrev = val;
                        bestP = pIdx;
                    }
                }
            }
        }

        scoreMat[idx(0, u)] = maxPrev;
        traceMat[idx(0, u)] = 3; // Delete
        tracePred[idx(0, u)] = bestP; 
    }

    // Main Loop
    for (let i = 1; i <= seqLen; i++) {
        const char = seq[i-1];
        
        for (let u = 0; u < subNodeCount; u++) {
            const uId = graphNodes[u];
            const nodeChar = fullGraph.nodes[uId].char;
            const match = (char === nodeChar) ? matchCost : mismatchCost;
            const preds = fullGraph.nodes[uId].prev;
            
            // 1. Match/Mismatch (Diag) from predecessors
            let maxDiag = -1e9;
            let bestDiagP = -1;

            const isStartConnected = preds.some(p => !realIdToSubIdx.has(p));
            if (isStartConnected) {
                // If aligned to anchor, prev score is 0 (relative to block start)
                // But if i > 1, it implies insertion before this match?
                // Logic: M[i-1][anchor] is effectively handled.
                // If i=1, score is match.
                if (i === 1) maxDiag = match; 
            }

            for (const pId of preds) {
                if (realIdToSubIdx.has(pId)) {
                    const pIdx = realIdToSubIdx.get(pId);
                    const s = scoreMat[idx(i-1, pIdx)];
                    if (s > -1e9) {
                        if (s + match > maxDiag) {
                            maxDiag = s + match;
                            bestDiagP = pIdx;
                        }
                    }
                }
            }

            // 2. Insertion (Up) - Seq moves, Graph stays (Conceptually: new node after predecessors)
            // But in this DP grid (Seq x GraphNodes), "Up" means M[i-1][u].
            // This represents aligning seq[i] to node u AGAIN (unlikely) 
            // or aligning seq[i] AFTER node u but before next node.
            // In Sequence-to-Graph, Insertion is usually modeled as a state separate from the graph nodes.
            // Simplified: M[i-1][u] + gapExtend.
            // This effectively maps multiple seq chars to one graph node.
            let maxIns = scoreMat[idx(i-1, u)] + gapExtend;
            
            // 3. Deletion (Left) - Seq stays, Graph moves
            // From predecessors M[i][v] + gapExtend
            let maxDel = -1e9;
            let bestDelP = -1;
            
            // Check anchor for deletion?
            if (isStartConnected && i === 0) { /* Handled in init */ } 
            
            for (const pId of preds) {
                if (realIdToSubIdx.has(pId)) {
                    const pIdx = realIdToSubIdx.get(pId);
                    const s = scoreMat[idx(i, pIdx)];
                    if (s > -1e9) {
                        if (s + gapExtend > maxDel) {
                            maxDel = s + gapExtend;
                            bestDelP = pIdx;
                        }
                    }
                }
            }

            // Pick Max
            // Bias order: Match > Del > Ins
            let bestScore = maxDiag;
            let type = 1;
            let pred = bestDiagP;

            if (maxDel > bestScore) {
                bestScore = maxDel;
                type = 3;
                pred = bestDelP;
            }
            if (maxIns > bestScore) { // Open/Extend logic simplified here
                bestScore = maxIns;
                type = 2;
                pred = u; // Points to self in previous row
            }

            scoreMat[idx(i, u)] = bestScore;
            traceMat[idx(i, u)] = type;
            tracePred[idx(i, u)] = pred;
        }
    }

    // Traceback
    // Find best end score in last row
    let maxEndScore = -1e9;
    let currU = -1;
    for(let u=0; u<subNodeCount; u++) {
        if (scoreMat[idx(seqLen, u)] > maxEndScore) {
            maxEndScore = scoreMat[idx(seqLen, u)];
            currU = u;
        }
    }

    // Reconstruction
    // Path: Array of { type: 'MATCH'|'INS'|'DEL', char: string, nodeId: number }
    const path = [];
    let currI = seqLen;
    
    while (currI > 0 || currU !== -1) {
        // Stop if we hit the "virtual anchor" boundary
        if (currI === 0 && currU === -1) break; // Should be handled by loop conditions
        
        // If currU is -1, we have remaining sequence to insert?
        if (currU === -1) {
            path.push({ type: 'INS', char: seq[currI-1], nodeId: null });
            currI--;
            continue;
        }

        const type = traceMat[idx(currI, currU)];
        const pred = tracePred[idx(currI, currU)];

        if (currI === 0) {
            // Must be Deletion chain or start
            // If type is 0 or uninit, break
            if (scoreMat[idx(currI, currU)] <= -1e8) break; 
        }

        if (type === 1) { // Match
            path.push({ type: 'MATCH', char: seq[currI-1], nodeId: graphNodes[currU] });
            currI--;
            currU = pred; // Go to prev node
        } else if (type === 2) { // Ins (Seq moves, graph stays same node context)
            path.push({ type: 'INS', char: seq[currI-1], nodeId: graphNodes[currU] }); // Ins *after* node U? Or *at* node U?
            // In this simplified DP, Ins at [i][u] comes from [i-1][u]. 
            // It effectively means an insertion appearing "on top of" or "after" u.
            // We'll treat it as insertion aligned to u (bubble).
            currI--;
            // currU stays same
        } else if (type === 3) { // Del (Graph moves, seq stays)
            path.push({ type: 'DEL', char: null, nodeId: graphNodes[currU] });
            currU = pred;
        } else {
            // If we get stuck (boundary), consume remaining seq as insertions
            if (currI > 0) {
                path.push({ type: 'INS', char: seq[currI-1], nodeId: null });
                currI--;
            } else {
                break;
            }
        }
    }

    return path.reverse();
}

/**
 * Applies the alignment path to the graph.
 * Returns the ID of the last visited/created node.
 */
function applyAlignmentToGraph(graph, alignPath, subSeq, startNodeId, seqIdx) {
    let currNodeId = startNodeId;

    for (const op of alignPath) {
        if (op.type === 'MATCH') {
            const node = graph.nodes[op.nodeId];
            if (node.char === op.char) {
                // Perfect match, fuse
                node.seqs.add(seqIdx);
                if (currNodeId !== -1 && currNodeId !== op.nodeId) {
                    graph.addEdge(currNodeId, op.nodeId);
                }
                currNodeId = op.nodeId;
            } else {
                // Mismatch -> Create Bubble Node
                // Check if a node with this char already exists from currNodeId? (Merging bubbles)
                let existingNext = -1;
                // Simple heuristic: don't merge parallel bubbles yet, just create new.
                // Merging requires checking all next neighbors of currNodeId.
                const newNode = graph.createNode(op.char);
                newNode.seqs.add(seqIdx);
                if (currNodeId !== -1) graph.addEdge(currNodeId, newNode.id);
                // Also edge to next of original?
                // The path will handle the 'next' connection in subsequent steps.
                currNodeId = newNode.id;
            }
        } else if (op.type === 'INS') {
            // Insertion: Sequence has char, Graph has nothing (or we are extending bubble)
            const newNode = graph.createNode(op.char);
            newNode.seqs.add(seqIdx);
            if (currNodeId !== -1) graph.addEdge(currNodeId, newNode.id);
            currNodeId = newNode.id;
        } else if (op.type === 'DEL') {
            // Deletion: Sequence skips node.
            // We don't visit the node, but we need to ensure connectivity.
            // Effectively we just skip updating currNodeId to this node?
            // No, we need to bypass it. 
            // Wait, if the path says DELETE node X, it means we consumed X in graph but not in seq.
            // We just don't add seqIdx to node X.
            // But we do need to bridge the gap.
            // The next MATCH will add an edge from currNodeId to the MatchNode.
            // So we don't change currNodeId! 
            // EXCEPT: If we have consecutive deletions, we are traversing the graph.
            // We shouldn't add edges for deletions.
            // Just advance the graph pointer logic implicitly?
            // Actually, for graph connectivity tracking in `currNodeId`, we should NOT update it
            // if we are skipping the node.
            // HOWEVER, the DP path follows graph edges.
            // op.nodeId is the node we skipped.
            // We technically "passed" it.
        }
    }
    return currNodeId;
}

// --- MSA Generation ---

function generateMSAFromGraph(graph) {
    const topo = graph.getTopologicalSort();
    const seqCount = graph.sequences.length;
    
    // Heuristic: Assign columns.
    // Ideally, perform graph coloring or critical path analysis.
    // Simple greedy approach:
    // Iterate topo order. Assign column C. 
    // If node u has edge to v, col(v) > col(u).
    const colMap = new Int32Array(graph.nodeCount).fill(0);
    
    for (const u of topo) {
        let minCol = colMap[u]; // initialized 0 or from predecessors
        for (const vId of graph.nodes[u].next) {
            if (colMap[vId] <= minCol) {
                colMap[vId] = minCol + 1;
            }
        }
    }

    const maxCol = Math.max(...colMap);
    const msa = Array(seqCount).fill().map(() => Array(maxCol + 1).fill('-'));

    // Fill MSA
    // We need to trace each sequence through the graph.
    // Since we didn't store full paths, we assume the graph structure + seqs set is enough?
    // No, branching makes it ambiguous which path a seq took if multiple paths have same char.
    // BUT: We stored `seqs` set on nodes.
    // If the graph is a proper POA, a sequence shouldn't split?
    // Actually, if we have identical parallel branches, it's ambiguous.
    // Solution: Re-align sequences to the final graph? (Standard trick).
    // OR: Track path during construction. (Hard with progressive).
    
    // Let's use the node occupancy.
    // Problem: If bubble has Top path A-B-C and Bot path A-X-C.
    // Seq 1 takes top, Seq 2 takes bot.
    // Topo sort: A, B, X, C.
    // Cols: A=0, B=1, X=1 (if parallel), C=2.
    // Correct.
    
    // Issue: colMap above makes X=1 and B=1 if they both come from A.
    // But how do we know row placement?
    // We just fill msa[seqIdx][colMap[u]] = char.
    
    for (const u of topo) {
        const node = graph.nodes[u];
        const col = colMap[u];
        for (const seqIdx of node.seqs) {
            // Collision check: if msa[seqIdx][col] is occupied?
            // This happens if bad topological linearization (e.g. B and X share column).
            // If collision, shift X to col+1? 
            // This requires dynamic column expansion.
            
            // Lazy fix: Just overwrite (bad) or shift global columns (complex).
            // Better greedy column assignment:
            // When visiting u, set col(u) = max(prev_cols) + 1.
            // Ensure no two nodes visited by SAME sequence share a column.
        }
    }
    
    // Improved Column Assignment
    // 1. Calculate length of path for each sequence.
    // 2. This is hard.
    
    // Let's stick to the "Re-trace" method. 
    // It guarantees correctness.
    // Since we have the graph and the original strings, we can just run 
    // a simplified "Trace" (exact match finder) for each sequence through the graph.
    // Because the sequence *constructed* the graph, a perfect path exists.
    
    const finalMSA = graph.sequences.map(() => []);
    
    // Columns logic needs to be global.
    // Let's build a "Global Column" list.
    // Each entry in list is a Set of Node IDs.
    // Nodes in same column cannot share a sequence.
    
    const columns = []; // Array<Set<NodeId>>
    const nodeToCol = new Int32Array(graph.nodeCount).fill(-1);

    for (const u of topo) {
        const node = graph.nodes[u];
        
        // Find earliest valid column
        // Must be > column of all predecessors
        // Must not contain any node v where intersect(node.seqs, v.seqs) is not empty
        
        let validCol = 0;
        
        // Pred constraint
        for(const p of node.prev) {
            if (nodeToCol[p] >= validCol) {
                validCol = nodeToCol[p] + 1;
            }
        }
        
        // Conflict constraint
        while (true) {
            if (validCol >= columns.length) {
                columns.push(new Set()); // Create new column
                break; // Valid
            }
            
            const colNodes = columns[validCol];
            let conflict = false;
            for(const existingId of colNodes) {
                const existing = graph.nodes[existingId];
                // Check intersection
                for(const s of node.seqs) {
                    if(existing.seqs.has(s)) {
                        conflict = true; 
                        break;
                    }
                }
                if(conflict) break;
            }
            
            if (!conflict) break;
            validCol++;
        }
        
        columns[validCol].add(u);
        nodeToCol[u] = validCol;
    }

    // Output
    const numCols = columns.length;
    const resultStrs = graph.sequences.map(() => new Array(numCols).fill('-'));

    for(let u=0; u<graph.nodeCount; u++) {
        const col = nodeToCol[u];
        const char = graph.nodes[u].char;
        for(const s of graph.nodes[u].seqs) {
            resultStrs[s][col] = char;
        }
    }

    return resultStrs.map(arr => arr.join(''));
}


// --- API ---

function multiSequenceAlign(sequences) {
    if (!sequences || sequences.length === 0) return [];
    if (sequences.length === 1) return [sequences[0]];

    const graph = new POAGraph();
    
    // Params: gapOpen, gapExtend, match, mismatch, boundary
    const nwParams = [-10.0, -1.0, 1.0, -1.0, 10]; 

    // Init with first
    graph.initFirstSequence(sequences[0]);

    // Progressive Align
    for (let i = 1; i < sequences.length; i++) {
        addSequenceToGraph(graph, sequences[i], nwParams);
    }

    // Generate MSA
    return generateMSAFromGraph(graph);
}


/**
 * Refined MSA Wrapper
 * Wraps multiSequenceAlign with an iterative refinement process.
 * 
 * 1. Generates initial MSA.
 * 2. Iteratively identifies "gappy" blocks with stable flanks.
 * 3. Re-aligns sequences within those blocks to a local reference using Double DP.
 * 4. Stitches blocks back and cleans up.
 */
function refinedMSA(sequences, iterations = 3) {
    // 1. Initial Alignment
    let currentMSA = multiSequenceAlign(sequences);

    if (!currentMSA || currentMSA.length === 0) return [];
    if (currentMSA.length === 1) return currentMSA;

    for (let iter = 0; iter < iterations; iter++) {
        // 2. Remove gap-only columns first to sanitize input for this iteration
        currentMSA = removeGapOnlyColumns(currentMSA);

        // 3. Identify Blocks (Gappy regions + Flanks)
        // Returns list of ranges: [{start, end, type: 'refine'|'stable'}]
        const blocks = identifyRefinementBlocks(currentMSA, 20); // 20bp flank

        const nextMSA = currentMSA.map(() => "");

        for (const block of blocks) {
            // Extract the slice for this block
            const blockSlice = currentMSA.map(seq => seq.substring(block.start, block.end));

            if (block.type === 'stable') {
                // Just copy stable blocks as-is
                for (let i = 0; i < nextMSA.length; i++) {
                    nextMSA[i] += blockSlice[i];
                }
            } else {
                // 4. Refine the Block
                const refinedSlice = refineBlockSlice(blockSlice);
                for (let i = 0; i < nextMSA.length; i++) {
                    nextMSA[i] += refinedSlice[i];
                }
            }
        }

        currentMSA = nextMSA;
    }

    // Final Cleanup
    return removeGapOnlyColumns(currentMSA);
}

/**
 * Identifies regions needing refinement.
 * A column is "unstable" if it contains gaps.
 * Contiguous unstable columns are grouped and padded by `flankSize`.
 * Overlapping regions are merged.
 */
function identifyRefinementBlocks(msa, flankSize) {
    const len = msa[0].length;
    const isUnstable = new Uint8Array(len).fill(0);

    // 1. Mark unstable columns (contain any gap)
    for (let j = 0; j < len; j++) {
        for (let i = 0; i < msa.length; i++) {
            if (msa[i][j] === '-') {
                isUnstable[j] = 1;
                break;
            }
        }
    }

    // 2. Dilation (Add flanks)
    // We want to group gaps that are close, and ensure we grab stable context.
    // Simple approach: Create a mask dilated by flankSize.
    const mask = new Uint8Array(len).fill(0);
    for (let j = 0; j < len; j++) {
        if (isUnstable[j]) {
            const start = Math.max(0, j - flankSize);
            const end = Math.min(len, j + flankSize + 1);
            for (let k = start; k < end; k++) {
                mask[k] = 1;
            }
        }
    }

    // 3. Convert mask to block ranges
    const blocks = [];
    let currentStart = 0;
    let inRefineBlock = (mask[0] === 1);

    for (let j = 1; j < len; j++) {
        const isRefine = (mask[j] === 1);
        if (isRefine !== inRefineBlock) {
            blocks.push({
                start: currentStart,
                end: j,
                type: inRefineBlock ? 'refine' : 'stable'
            });
            currentStart = j;
            inRefineBlock = isRefine;
        }
    }
    // Final block
    blocks.push({
        start: currentStart,
        end: len,
        type: inRefineBlock ? 'refine' : 'stable'
    });

    return blocks;
}

/**
 * Refines a specific vertical slice of the MSA.
 * 1. Ungaps sequences to get raw content.
 * 2. Picks longest sequence as reference.
 * 3. Aligns others to reference using Pairwise DoubleDP.
 * 4. Merges results (Left-piling insertions).
 */
function refineBlockSlice(msaSlice) {
    const numSeqs = msaSlice.length;
    
    // 1. Ungap strings
    const rawSeqs = msaSlice.map(s => s.split('-').join(''));
    
    // Optimization: If block was purely gaps or empty for some reason, return original
    if (rawSeqs.every(s => s.length === 0)) return msaSlice;

    // 2. Find Reference (Sequence with max non-gap characters)
    let bestIdx = 0;
    let maxLen = -1;
    for (let i = 0; i < numSeqs; i++) {
        if (rawSeqs[i].length > maxLen) {
            maxLen = rawSeqs[i].length;
            bestIdx = i;
        }
    }
    const refSeq = rawSeqs[bestIdx];

    // Data structures for Star Alignment
    // Normalized to Ref Coordinates: 0 to refSeq.length
    // refMatches[pos][seqIdx] = char (or '-' if deleted)
    // insertions[pos][seqIdx] = string (inserted BEFORE pos)
    
    const refLen = refSeq.length;
    const refMatches = Array(refLen).fill().map(() => Array(numSeqs).fill('-'));
    const insertions = Array(refLen + 1).fill().map(() => Array(numSeqs).fill(''));

    // 3. Pairwise Align everyone to Reference
    // NW Params reused from library defaults
    const nwParams = [-10.0, -0.2, 1.0, -0.7, 10]; 

    for (let i = 0; i < numSeqs; i++) {
        if (i === bestIdx) {
            // Align Ref to Ref (Trivial)
            for (let k = 0; k < refLen; k++) refMatches[k][i] = refSeq[k];
            continue;
        }

        const query = rawSeqs[i];
        if (query.length === 0) continue; // Leaves row as all gaps

        // Align (Ref, Query)
        // doubleDP_nwalign returns [alignedRef, alignedQuery]
        const [alnRef, alnQuery] = doubleDP_nwalign(refSeq, query, ...nwParams);

        // 4. Map Alignment to Star Structure
        let refPos = 0;
        let activeIns = ""; // Buffer for insertion string being built

        for (let k = 0; k < alnRef.length; k++) {
            const rChar = alnRef[k];
            const qChar = alnQuery[k];

            if (rChar !== '-') {
                // Flush any active insertion to the bucket BEFORE this ref char
                if (activeIns.length > 0) {
                    insertions[refPos][i] = activeIns;
                    activeIns = "";
                }
                
                // Record Match/Mismatch/Delete state
                // If qChar is '-', it's a deletion relative to ref
                refMatches[refPos][i] = qChar;
                refPos++;
            } else {
                // Gap in Ref -> Insertion in Query
                // Accumulate insertion string
                if (qChar !== '-') {
                    activeIns += qChar;
                }
            }
        }
        // Flush trailing insertion (after last ref char)
        if (activeIns.length > 0) {
            insertions[refPos][i] = activeIns;
        }
    }

    // 5. Reconstruct Matrix
    const finalCols = []; // array of strings (columns)

    for (let k = 0; k <= refLen; k++) {
        // A. Handle Insertions at this boundary
        // Find max insertion length at this position
        let maxInsLen = 0;
        for (let i = 0; i < numSeqs; i++) {
            if (insertions[k][i].length > maxInsLen) maxInsLen = insertions[k][i].length;
        }

        if (maxInsLen > 0) {
            for (let pos = 0; pos < maxInsLen; pos++) {
                let colStr = "";
                for (let i = 0; i < numSeqs; i++) {
                    const insStr = insertions[k][i];
                    // Left pile: char at pos, or gap
                    const char = (pos < insStr.length) ? insStr[pos] : '-';
                    colStr += char;
                }
                finalCols.push(colStr);
            }
        }

        // B. Handle Reference Position Match (if not at end)
        if (k < refLen) {
            let colStr = "";
            for (let i = 0; i < numSeqs; i++) {
                colStr += refMatches[k][i];
            }
            finalCols.push(colStr);
        }
    }

    // Transpose back to row-strings
    const newSlice = Array(numSeqs).fill("");
    for (const col of finalCols) {
        for (let i = 0; i < numSeqs; i++) {
            newSlice[i] += col[i];
        }
    }

    return newSlice;
}

/**
 * Removes columns that contain only gaps.
 */
function removeGapOnlyColumns(msa) {
    if (msa.length === 0) return [];
    const len = msa[0].length;
    const keepCol = new Uint8Array(len).fill(0);

    for (let j = 0; j < len; j++) {
        for (let i = 0; i < msa.length; i++) {
            if (msa[i][j] !== '-') {
                keepCol[j] = 1;
                break;
            }
        }
    }

    const newMSA = msa.map(() => "");
    for (let j = 0; j < len; j++) {
        if (keepCol[j]) {
            for (let i = 0; i < msa.length; i++) {
                newMSA[i] += msa[i][j];
            }
        }
    }
    return newMSA;
}