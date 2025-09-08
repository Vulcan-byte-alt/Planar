class PlanarGraph {
            constructor() {
                this.vertices = [];
                this.edges = [];
                this.periphery = [];
                this.nextId = 1;
                this.maxVertices = 10000;
                this.displayShowIndex = false;
                this.selectedColor = 1;
                this.colors = {
                    1: '#FF6B6B',
                    2: '#4ECDC4',
                    3: '#45B7D1',
                    4: '#96CEB4'
                };
            }

            addVertex(x, y, color = 1, touchingVertices = []) {
                if (this.vertices.length >= this.maxVertices) {
                    throw new Error('Maximum vertices exceeded');
                }


                const vertex = {
                    id: this.nextId++,
                    x: x,
                    y: y,
                    color: color,
                    radius: this.calculateRadius(),
                    touchingVertices: touchingVertices
                };

                this.vertices.push(vertex);

                this.updatePeriphery(vertex, touchingVertices);
                return vertex;
            }

            calculateRadius() {
                const baseRadius = 15;
                const digitCount = this.nextId.toString().length;
                return baseRadius + (digitCount - 1) * 3;
            }

            addEdge(v1Id, v2Id) {
                // Planarity check: do not add if this edge would cross any existing edge
                const v1 = this.getVertex(v1Id);
                const v2 = this.getVertex(v2Id);
                for (const edge of this.edges) {
                    const a = this.getVertex(edge.v1);
                    const b = this.getVertex(edge.v2);
                    if (
                        (v1Id !== edge.v1 && v1Id !== edge.v2 && v2Id !== edge.v1 && v2Id !== edge.v2) &&
                        PlanarGraph.segmentsIntersect(v1, v2, a, b)
                    ) {
                        // Edge would cross another edge - show warning
                        console.warn(`Edge (${v1Id}, ${v2Id}) would cross existing edge (${edge.v1}, ${edge.v2}) and cannot be added.`);
                        updateStatus(`Edge (${v1Id}, ${v2Id}) would cross existing edge (${edge.v1}, ${edge.v2}) and cannot be added.`);
                        // --- IGNORE
                        return null;
                    }
                }
                const edge = { v1: v1Id, v2: v2Id };
                this.edges.push(edge);
                return edge;
            }

            static segmentsIntersect(p1, p2, q1, q2) {
                // Helper for planarity check
                function ccw(a, b, c) {
                    return (c.y - a.y) * (b.x - a.x) > (b.y - a.y) * (c.x - a.x);
                }
                return (
                    ccw(p1, q1, q2) !== ccw(p2, q1, q2) &&
                    ccw(p1, p2, q1) !== ccw(p1, p2, q2)
                );
            }

            updatePeriphery(newVertex, touchingVertices) {
                if (touchingVertices.length < 2) return;

                // Only update if at least two touching vertices (Vp, Vq) are given
                const vpId = touchingVertices[0];
                const vqId = touchingVertices[touchingVertices.length - 1];
                const vpIndex = this.periphery.indexOf(vpId);
                const vqIndex = this.periphery.indexOf(vqId);
                if (vpIndex === -1 || vqIndex === -1) return;

                // Build new periphery: keep Vp, insert newVertex after Vp, keep Vq, remove all between Vp and Vq (clockwise)
                let newPeriphery = [];
                let i = vpIndex;
                newPeriphery.push(this.periphery[i]); // Vp
                newPeriphery.push(newVertex.id); // Insert new vertex after Vp
                i = (i + 1) % this.periphery.length;
                // Skip all vertices between Vp and Vq (exclusive)
                while (i !== vqIndex) {
                    i = (i + 1) % this.periphery.length;
                    if (i === vqIndex) break;
                }
                newPeriphery.push(this.periphery[vqIndex]); // Vq
                // Add the rest of the periphery (after Vq, if any)
                i = (vqIndex + 1) % this.periphery.length;
                while (i !== vpIndex) {
                    newPeriphery.push(this.periphery[i]);
                    i = (i + 1) % this.periphery.length;
                }
                this.periphery = newPeriphery;
            }

            getVertex(id) {
                return this.vertices.find(v => v.id === id);
            }

            isValidVertexSelection(vpId, vqId) {
                const vpIndex = this.periphery.indexOf(vpId);
                const vqIndex = this.periphery.indexOf(vqId);

                if (vpIndex === -1 || vqIndex === -1) return false;
                if (vpIndex === vqIndex) return false;

                return true;
            }

            getVerticesBetween(vpId, vqId) {
                const vpIndex = this.periphery.indexOf(vpId);
                const vqIndex = this.periphery.indexOf(vqId);

                if (vpIndex === -1 || vqIndex === -1) return [];

                const vertices = [];
                let current = vpIndex;

                while (current !== vqIndex) {
                    vertices.push(this.periphery[current]);
                    current = (current + 1) % this.periphery.length;
                }
                vertices.push(this.periphery[vqIndex]);

                return vertices;
            }

            /**
             * Returns the list of vertex IDs between Vp and Vq (inclusive) in both directions.
             * @param {number} vpId - Start vertex ID
             * @param {number} vqId - End vertex ID
             * @returns {{cw: number[], ccw: number[]}} - Object with clockwise and counterclockwise arrays
             */
            getVerticesBetweenBothDirections(vpId, vqId) {
                const periphery = [...this.periphery];
                console.log('getVerticesBetweenBothDirections', vpId, vqId, periphery);
                const n = periphery.length;
                const vpIndex = periphery.indexOf(vpId);
                const vqIndex = periphery.indexOf(vqId);
                if (vpIndex === -1 || vqIndex === -1) return { cw: [], ccw: [], shortest: [], longest: [] };
                // Clockwise
                let cw = [];
                let i = vpIndex;
                while (i !== vqIndex) {
                    cw.push(periphery[i]);
                    i = (i + 1) % n;
                }
                cw.push(periphery[vqIndex]);
                // Counterclockwise
                let ccw = [];
                i = vpIndex;
                while (i !== vqIndex) {
                    ccw.push(periphery[i]);
                    i = (i - 1 + n) % n;
                }
                ccw.push(periphery[vqIndex]);
                // Shortest
                let shortest = cw.length <= ccw.length ? cw : ccw;
                let longest = cw.length > ccw.length ? cw : ccw;
                console.log('Vertices between both directions:', {
                    cw: cw,
                    ccw: ccw,
                    shortest: shortest,
                    longest: longest
                });
                return { cw, ccw, shortest, longest };
            }
        }

        class GraphRenderer {
            constructor(canvas, graph) {
                this.canvas = canvas;
                this.ctx = canvas.getContext('2d');
                this.graph = graph;
                this.zoom = 1;
                this.offsetX = 0;
                this.offsetY = 0;
                this.isDragging = false;
                this.lastMouseX = 0;
                this.lastMouseY = 0;
                this.maxVisibleVertex = Infinity;
                this.isForceDirected = true; // Enable force-directed layout by default
                this.simulationSteps = 0;
                this.simulationMaxSteps = 1000;
                this.simulationRunning = false;
                this.simulationInterval = null;
                this.numVertexSolutions = 3; // Configurable number of candidate solutions

                this.resizeCanvas();
                this.setupEventListeners();
            }

            addVertexBetween(vpId, vqId) {
                console.log(`Adding vertex between ${vpId} and ${vqId}`);
                // Find indices in periphery
                const periphery = this.graph.periphery;
                const vpIndex = periphery.indexOf(vpId);
                const vqIndex = periphery.indexOf(vqId);
                if (vpIndex === -1 || vqIndex === -1 || vpIndex === vqIndex) return;

                // Get vertices
                const vp = this.graph.getVertex(vpId);
                const vq = this.graph.getVertex(vqId);

                // Use shortest path between vp and vq (inclusive) from getVerticesBetweenBothDirections
                const { shortest, longest } = this.graph.getVerticesBetweenBothDirections(vpId, vqId);

                // Use the helper to find the new vertex position
                const pos = this.findNewVertexPositionBasedOffCentroid(vpId, vqId, this.numVertexSolutions);
                // const pos = this.findNewVertexPosition(vpId, vqId);
                if (!pos) return;
                const newX = pos.x;
                const newY = pos.y;
                console.log(`New vertex position: (${newX}, ${newY})`, vp, vq);


                // Use shortest path between vp and vq (inclusive) from getVerticesBetweenBothDirections
                // const { shortest: verticesBetween } = this.graph.getVerticesBetweenBothDirections(vpId, vqId);

                // // Remove edges between consecutive vertices in verticesBetween (periphery segment)
                // if (verticesBetween.length > 2) { // inclusive of Vp and Vq
                //     for (let i = 0; i < verticesBetween.length - 1; i++) {
                //         const vA = verticesBetween[i];
                //         const vB = verticesBetween[i + 1];
                //         // Remove edge vA-vB or vB-vA
                //         this.graph.edges = this.graph.edges.filter(e => {
                //             return !((e.v1 === vA && e.v2 === vB) || (e.v1 === vB && e.v2 === vA));
                //         });
                //     }
                //     // Also remove edge from last to first if periphery is closed (optional, only if needed)
                //     const vA = verticesBetween[verticesBetween.length - 1];
                //     const vB = verticesBetween[0];
                //     this.graph.edges = this.graph.edges.filter(e => {
                //         return !((e.v1 === vA && e.v2 === vB) || (e.v1 === vB && e.v2 === vA));
                //     });
                // }

                // Add new vertex
                const newVertex = this.graph.addVertex(newX, newY, this.graph.selectedColor, shortest);

                // Connect new vertex to all vertices in between (including vp and vq)
                for (const vId of shortest) {
                    this.graph.addEdge(newVertex.id, vId);
                }

                // Update periphery: remove verticesBetween, insert newVertex after vp
                // Update periphery: use all vertices from the longest path, then add the new vertex at the end
                // const { longest, shortest } = this.graph.getVerticesBetweenBothDirections(vpId, vqId);
                let newPeriphery = [...longest, newVertex.id];
                console.log('longest path:', longest, 'shortest path:', shortest, 'new vertex:', newVertex.id, 'new periphery:', newPeriphery);
                this.graph.periphery = newPeriphery;

                // Start force-directed rendering
                this.startForceDirectedSimulation();

                this.updateStats();
                this.draw();
            }

            resizeCanvas() {
                const container = this.canvas.parentElement;
                this.canvas.width = container.clientWidth;
                this.canvas.height = container.clientHeight;
                this.draw();
            }

            setupEventListeners() {
                window.addEventListener('resize', () => this.resizeCanvas());

                this.canvas.addEventListener('mousedown', (e) => {
                    this.isDragging = true;
                    this.lastMouseX = e.clientX;
                    this.lastMouseY = e.clientY;
                });

                this.canvas.addEventListener('mousemove', (e) => {
                    if (this.isDragging) {
                        const deltaX = e.clientX - this.lastMouseX;
                        const deltaY = e.clientY - this.lastMouseY;
                        this.offsetX += deltaX;
                        this.offsetY += deltaY;
                        this.lastMouseX = e.clientX;
                        this.lastMouseY = e.clientY;
                        this.draw();
                    }
                });

                this.canvas.addEventListener('mouseup', () => {
                    this.isDragging = false;
                });

                this.canvas.addEventListener('click', (e) => {
                    if (currentMode === 'selectVp') {
                        this.selectVertex(e, 'vp');
                    } else if (currentMode === 'selectVq') {
                        this.selectVertex(e, 'vq');
                    }
                });

                this.canvas.addEventListener('wheel', (e) => {
                    e.preventDefault();
                    const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
                    this.zoom *= zoomFactor;
                    this.draw();
                });
            }

            findNewVertexPositionBasedOffCentroid(vpId, vqId, numSolutions = 3) {
                const vp = this.graph.getVertex(vpId);
                const vq = this.graph.getVertex(vqId);
                if (!vp || !vq) {
                    console.error(`Invalid vertex IDs: ${vpId}, ${vqId}`);
                    return null;
                }

                // Get the shortest path between vp and vq (this is the chord)
                const { shortest, longest } = this.graph.getVerticesBetweenBothDirections(vpId, vqId);
                if (shortest.length < 2) {
                    console.error(`Not enough vertices between ${vpId} and ${vqId}`);
                    return null;
                }

                // Calculate the midpoint of the chord (vp to vq)
                const midX = (vp.x + vq.x) / 2;
                const midY = (vp.y + vq.y) / 2;

                // Calculate the perpendicular direction to the chord
                const chordDx = vq.x - vp.x;
                const chordDy = vq.y - vp.y;
                const chordLength = Math.sqrt(chordDx * chordDx + chordDy * chordDy);

                // Generate multiple candidate solutions
                const candidates = this.generateCandidatePositions(vp, vq, midX, midY, chordLength, numSolutions);

                // Filter valid candidates (no edge crossings, outside polygon, no peripheral edge crossings)
                const validCandidates = candidates.filter(candidate => {
                    return !this.wouldCreateCrossingEdges(candidate.x, candidate.y, shortest) &&
                        !this.isPointInsidePolygon(candidate.x, candidate.y) &&
                        !this.wouldCrossPeripheralEdges(candidate.x, candidate.y, shortest);
                });

                if (validCandidates.length === 0) {
                    console.warn('No valid candidates found, using fallback position');
                    return this.getFallbackPosition(vp, vq, midX, midY, chordLength);
                }

                // Choose the best candidate based on distance to chord midpoint
                const bestCandidate = this.selectBestCandidate(validCandidates, midX, midY);

                console.log(`Generated ${candidates.length} candidates, ${validCandidates.length} valid, selected:`, bestCandidate);
                return bestCandidate;
            }

            generateCandidatePositions(vp, vq, midX, midY, chordLength, numSolutions) {
                const candidates = [];

                // Calculate perpendicular directions to the chord
                const chordDx = vq.x - vp.x;
                const chordDy = vq.y - vp.y;
                const perpX = -chordDy;
                const perpY = chordDx;
                const perpLength = Math.sqrt(perpX * perpX + perpY * perpY);
                const normPerpX = perpX / perpLength;
                const normPerpY = perpY / perpLength;

                // Determine primary outside direction
                const primaryDirection = this.determineOutsideDirection(vp.id, vq.id, normPerpX, normPerpY);

                // Generate candidates with varying strategies
                for (let i = 0; i < numSolutions; i++) {
                    let candidate;

                    if (i === 0) {
                        // Primary solution: perpendicular to chord in outside direction
                        const baseDistance = chordLength * 0.5;
                        candidate = {
                            x: midX + primaryDirection.x * baseDistance,
                            y: midY + primaryDirection.y * baseDistance,
                            strategy: 'perpendicular-primary',
                            distance: baseDistance
                        };
                    } else if (i === 1) {
                        // Secondary solution: farther out in same direction
                        const baseDistance = chordLength * 0.8;
                        candidate = {
                            x: midX + primaryDirection.x * baseDistance,
                            y: midY + primaryDirection.y * baseDistance,
                            strategy: 'perpendicular-far',
                            distance: baseDistance
                        };
                    } else {
                        // Additional solutions: slight angular variations
                        const angle = (i - 2) * (Math.PI / 12); // 15-degree increments
                        const rotatedDirX = primaryDirection.x * Math.cos(angle) - primaryDirection.y * Math.sin(angle);
                        const rotatedDirY = primaryDirection.x * Math.sin(angle) + primaryDirection.y * Math.cos(angle);
                        const baseDistance = chordLength * (0.5 + (i - 2) * 0.1);

                        candidate = {
                            x: midX + rotatedDirX * baseDistance,
                            y: midY + rotatedDirY * baseDistance,
                            strategy: `angled-${i}`,
                            distance: baseDistance
                        };
                    }

                    candidates.push(candidate);
                }

                return candidates;
            }

            selectBestCandidate(validCandidates, midX, midY) {
                // Score candidates based on multiple criteria
                const scoredCandidates = validCandidates.map(candidate => {
                    const distanceFromMidpoint = Math.sqrt(
                        (candidate.x - midX) ** 2 + (candidate.y - midY) ** 2
                    );

                    // Calculate how far the candidate is from existing vertices (prefer more spacing)
                    let minDistanceToVertex = Infinity;
                    for (const vertex of this.graph.vertices) {
                        const dist = Math.sqrt((candidate.x - vertex.x) ** 2 + (candidate.y - vertex.y) ** 2);
                        minDistanceToVertex = Math.min(minDistanceToVertex, dist);
                    }

                    // Scoring: prefer moderate distance from midpoint, good spacing from other vertices
                    const distanceScore = 1 / (1 + Math.abs(distanceFromMidpoint - candidate.distance));
                    const spacingScore = Math.min(minDistanceToVertex / 50, 1); // Normalize to [0,1]
                    const strategyBonus = candidate.strategy === 'perpendicular-primary' ? 0.1 : 0;

                    const totalScore = distanceScore * 0.4 + spacingScore * 0.5 + strategyBonus;

                    return {
                        ...candidate,
                        score: totalScore,
                        distanceFromMidpoint,
                        minDistanceToVertex
                    };
                });

                // Sort by score (highest first) and return the best
                scoredCandidates.sort((a, b) => b.score - a.score);
                return scoredCandidates[0];
            }

            getFallbackPosition(vp, vq, midX, midY, chordLength) {
                // Calculate perpendicular direction and place vertex far outside
                const chordDx = vq.x - vp.x;
                const chordDy = vq.y - vp.y;
                const perpX = -chordDy;
                const perpY = chordDx;
                const perpLength = Math.sqrt(perpX * perpX + perpY * perpY);
                const normPerpX = perpX / perpLength;
                const normPerpY = perpY / perpLength;

                const outsideDirection = this.determineOutsideDirection(vp.id, vq.id, normPerpX, normPerpY);
                const fallbackDistance = chordLength * 2; // Large distance to ensure it's outside

                return {
                    x: midX + outsideDirection.x * fallbackDistance,
                    y: midY + outsideDirection.y * fallbackDistance,
                    strategy: 'fallback',
                    distance: fallbackDistance
                };
            }

            wouldCreateCrossingEdges(newX, newY, touchingVertices) {
                // Check if edges from new vertex to touching vertices would cross existing edges
                for (const vId of touchingVertices) {
                    const touchingVertex = this.graph.getVertex(vId);
                    if (!touchingVertex) continue;

                    // Check this potential edge against all existing edges
                    for (const edge of this.graph.edges) {
                        const edgeV1 = this.graph.getVertex(edge.v1);
                        const edgeV2 = this.graph.getVertex(edge.v2);
                        if (!edgeV1 || !edgeV2) continue;

                        // Skip if the existing edge shares a vertex with our potential edge
                        if (edge.v1 === vId || edge.v2 === vId) continue;

                        // Check for intersection
                        if (PlanarGraph.segmentsIntersect(
                            { x: newX, y: newY },
                            touchingVertex,
                            edgeV1,
                            edgeV2
                        )) {
                            return true; // Would create crossing
                        }
                    }
                }
                return false; // No crossings detected
            }

            wouldCrossPeripheralEdges(newX, newY, touchingVertices) {
                // Check if edges from new vertex to touching vertices would cross peripheral edges
                for (const vId of touchingVertices) {
                    const touchingVertex = this.graph.getVertex(vId);
                    if (!touchingVertex) continue;

                    // Check this potential edge against all peripheral edges
                    for (let i = 0; i < this.graph.periphery.length; i++) {
                        const currentId = this.graph.periphery[i];
                        const nextId = this.graph.periphery[(i + 1) % this.graph.periphery.length];

                        const currentVertex = this.graph.getVertex(currentId);
                        const nextVertex = this.graph.getVertex(nextId);
                        if (!currentVertex || !nextVertex) continue;

                        // Skip if the peripheral edge shares a vertex with our potential edge
                        if (currentId === vId || nextId === vId) continue;

                        // Check for intersection between potential edge and peripheral edge
                        if (PlanarGraph.segmentsIntersect(
                            { x: newX, y: newY },
                            touchingVertex,
                            currentVertex,
                            nextVertex
                        )) {
                            return true; // Would cross peripheral edge
                        }
                    }
                }
                return false; // No peripheral edge crossings detected
            }

            /**
             * Checks if adding a new vertex at (newX, newY) connected to touchingVertices would cause any touching vertex to have more than 4 edges.
             * @param {number} newX
             * @param {number} newY
             * @param {number[]} touchingVertices
             * @returns {boolean} True if any touching vertex would exceed 4 edges
             */
            wouldExceedMaxEdges(newX, newY, touchingVertices, maxEdges = 4) {
                for (const vId of touchingVertices) {
                    let edgeCount = 0;
                    for (const edge of this.graph.edges) {
                        if (edge.v1 === vId || edge.v2 === vId) edgeCount++;
                    }
                    // Adding the new edge from new vertex to vId
                    if (edgeCount + 1 > maxEdges) {
                        return true;
                    }
                }
                return false;
            }

            determineOutsideDirection(vpId, vqId, normPerpX, normPerpY) {
                // Get the vertices that form the chord
                const vp = this.graph.getVertex(vpId);
                const vq = this.graph.getVertex(vqId);

                // Get the shortest path (the chord we're replacing)
                const { shortest } = this.graph.getVerticesBetweenBothDirections(vpId, vqId);

                // If the shortest path has more than 2 vertices, we need to place the new vertex
                // on the side opposite to the intermediate vertices
                if (shortest.length > 2) {
                    // Take a vertex in the middle of the shortest path
                    const midIndex = Math.floor(shortest.length / 2);
                    const intermediateVertex = this.graph.getVertex(shortest[midIndex]);

                    // Calculate midpoint of chord
                    const chordMidX = (vp.x + vq.x) / 2;
                    const chordMidY = (vp.y + vq.y) / 2;

                    // Vector from chord midpoint to intermediate vertex
                    const toIntermediateX = intermediateVertex.x - chordMidX;
                    const toIntermediateY = intermediateVertex.y - chordMidY;

                    // Determine which perpendicular direction points away from intermediate vertex
                    const dot1 = normPerpX * toIntermediateX + normPerpY * toIntermediateY;
                    const dot2 = -normPerpX * toIntermediateX + -normPerpY * toIntermediateY;

                    // Choose the direction that points away from the intermediate vertex
                    if (dot1 > 0) {
                        return { x: -normPerpX, y: -normPerpY }; // Opposite direction
                    } else {
                        return { x: normPerpX, y: normPerpY };
                    }
                } else {
                    // For a simple chord (only 2 vertices), use centroid method as fallback
                    const allVertices = this.graph.vertices;
                    let centroidX = 0, centroidY = 0;
                    for (const v of allVertices) {
                        centroidX += v.x;
                        centroidY += v.y;
                    }
                    centroidX /= allVertices.length;
                    centroidY /= allVertices.length;

                    const chordMidX = (vp.x + vq.x) / 2;
                    const chordMidY = (vp.y + vq.y) / 2;

                    const toCentroidX = centroidX - chordMidX;
                    const toCentroidY = centroidY - chordMidY;

                    const dotProduct = normPerpX * toCentroidX + normPerpY * toCentroidY;

                    return dotProduct > 0 ?
                        { x: -normPerpX, y: -normPerpY } :
                        { x: normPerpX, y: normPerpY };
                }
            }

            isPointInsidePolygon(x, y) {
                console.log(`Checking if point (${x}, ${y}) is inside polygon...`);
                // Use ray casting algorithm for any polygon (convex or non-convex)
                const polygonVertices = this.graph.periphery.map(id => this.graph.getVertex(id));
                if (polygonVertices.length < 3) return false;

                let inside = false;
                const n = polygonVertices.length;

                for (let i = 0, j = n - 1; i < n; j = i++) {
                    const vi = polygonVertices[i];
                    const vj = polygonVertices[j];

                    if (((vi.y > y) !== (vj.y > y)) &&
                        (x < (vj.x - vi.x) * (y - vi.y) / (vj.y - vi.y) + vi.x)) {
                        inside = !inside;
                    }
                }
                console.log(`Point (${x}, ${y}) is ${inside ? 'inside' : 'outside'} the polygon.`);
                return inside;
            }

            findNewVertexPosition(vpId, vqId) {
                // find shorter and longer paths between vp and vq
                // find the vectors,
                // 1. Vp and Vq  => Vpq
                // 2. Vp-1 and Vp => Vp-1p in the direction of the shortest periphery
                // 3. Vq-1 and Vq => Vq-1q in the direction of the shortest periphery
                // Find the vector sum of Vp-1p and Vq-1q, then normalize it. This is the direction of the new vertex.
                // Find the middle point between Vp and Vq, then offset it in the direction of the vector sum.

                const vp = this.graph.getVertex(vpId);
                const vq = this.graph.getVertex(vqId);
                if (!vp || !vq) {
                    console.error(`Invalid vertex IDs: ${vpId}, ${vqId}`);
                    return null;
                }
                // Get the shortest path between vp and vq
                const { shortest, longest } = this.graph.getVerticesBetweenBothDirections(vpId, vqId);
                if (shortest.length < 2) {
                    console.error(`Not enough vertices between ${vpId} and ${vqId}`);
                    return null;
                }

                // Calculate direction vector along the Vp-1 and Vp in longest path
                const Vp = this.graph.getVertex(longest[0]);
                const Vq = this.graph.getVertex(longest[longest.length - 1]);
                const VpMinus1 = this.graph.getVertex(longest[longest.length - 2]);
                const VqMinus1 = this.graph.getVertex(longest[1]);
                if (!Vp || !Vq || !VpMinus1 || !VqMinus1) {
                    console.error(`Invalid vertices in longest path: ${longest}`);
                    return null;
                }

                // Calculate direction vectors
                const dirVp = { x: Vp.x - VpMinus1.x, y: Vp.y - VpMinus1.y };
                const dirVq = { x: Vq.x - VqMinus1.x, y: Vq.y - VqMinus1.y };
                // Normalize direction vectors
                const lenVp = Math.sqrt(dirVp.x * dirVp.x + dirVp.y * dirVp.y) || 1;
                const lenVq = Math.sqrt(dirVq.x * dirVq.x + dirVq.y * dirVq.y) || 1;

                dirVp.x /= lenVp;
                dirVp.y /= lenVp;

                dirVq.x /= lenVq;
                dirVq.y /= lenVq;

                // Calculate the midpoint between Vp and Vq
                const midX = (vp.x + vq.x) / 2;
                const midY = (vp.y + vq.y) / 2;

                // Offset by the distance between vp and vq in the direction of the shortest periphery
                const dist = Math.sqrt((vp.x - vq.x) ** 2 + (vp.y - vq.y) ** 2);
                const offset = dist * 0.5; // You can adjust this factor as needed

                // Calculate the new vertex position
                const newX = midX + (dirVp.x + dirVq.x) * offset;
                const newY = midY + (dirVp.y + dirVq.y) * offset;
                console.log(`New vertex position: (${newX}, ${newY})`, vp, vq);

                return { x: newX, y: newY };
            }

            selectVertex(e, type) {
                // Get mouse position relative to canvas
                const rect = this.canvas.getBoundingClientRect();
                const x = (e.clientX - rect.left - this.offsetX) / this.zoom;
                const y = (e.clientY - rect.top - this.offsetY) / this.zoom;

                // Find closest periphery vertex within a reasonable radius
                let minDist = Infinity;
                let selected = null;
                for (const vId of this.graph.periphery) {
                    const v = this.graph.getVertex(vId);
                    const dist = Math.sqrt((v.x - x) ** 2 + (v.y - y) ** 2);
                    if (dist < v.radius + 10 && dist < minDist) {
                        minDist = dist;
                        selected = v;
                    }
                }
                if (!selected) {
                    updateStatus('No periphery vertex selected. Try again.');
                    return;
                }

                if (type === 'vp') {
                    selectedVp = selected.id;
                    currentMode = 'selectVq';
                    updateStatus('Vp selected. Now select Vq (second vertex) from the periphery.');
                    updateModeIndicator('Select Vq');
                } else if (type === 'vq') {
                    selectedVq = selected.id;
                    if (!this.graph.isValidVertexSelection(selectedVp, selectedVq)) {
                        updateStatus('Invalid selection. Vq must be different and on the periphery.');
                        return;
                    }
                    // Instead of addVertexBetween, use findRandomVertexPositionOutsidePolygon
                    this.findRandomVertexPositionOutsidePolygon(selectedVp, selectedVq);
                    selectedVp = null;
                    selectedVq = null;
                    currentMode = 'ready';
                    updateStatus('Random vertex added outside periphery between selected vertices.');
                    updateModeIndicator('Ready');
                }
                this.draw();
            }

            startForceDirectedSimulation() {
                // return;
                // Reset all vertex velocities to zero before starting simulation
                for (const v of this.graph.vertices) {
                    v.vx = 0;
                    v.vy = 0;
                }
                if (this.simulationRunning) return;
                this.simulationRunning = true;
                this.simulationSteps = 0;
                this.simulationInterval = setInterval(() => {
                    this.forceDirectedStep();
                    this.draw();
                    this.simulationSteps++;
                    if (this.simulationSteps > this.simulationMaxSteps) {
                        clearInterval(this.simulationInterval);
                        this.simulationRunning = false;
                    }
                }, 1000 / 60);
            }

            forceDirectedStep() {
                const vertices = this.graph.vertices;
                const edges = this.graph.edges;
                const k = 100; // Ideal edge length for internal edges
                const kPeriphery = 150; // Higher ideal length for periphery edges
                const repulsion = 25000; // Repulsion constant
                const damping = 0.55;
                const containmentStrength = 100;
                const angleResistanceStrength = 5;
                const idealAngle = Math.PI / 3;
                // Initialize velocity
                if (!vertices[0].vx) {
                    for (const v of vertices) {
                        v.vx = 0;
                        v.vy = 0;
                    }
                }
                // Repulsion between all pairs
                for (let i = 0; i < vertices.length; i++) {
                    for (let j = i + 1; j < vertices.length; j++) {
                        const v1 = vertices[i], v2 = vertices[j];
                        let dx = v2.x - v1.x;
                        let dy = v2.y - v1.y;
                        let dist = Math.sqrt(dx * dx + dy * dy) + 0.1;
                        let force = repulsion / (dist * dist);
                        let fx = force * dx / dist;
                        let fy = force * dy / dist;
                        v1.vx -= fx;
                        v1.vy -= fy;
                        v2.vx += fx;
                        v2.vy += fy;
                    }
                }
                // Spring force for edges (classic spring model)
                for (const edge of edges) {
                    const v1 = this.graph.getVertex(edge.v1);
                    const v2 = this.graph.getVertex(edge.v2);
                    let dx = v2.x - v1.x;
                    let dy = v2.y - v1.y;
                    let dist = Math.sqrt(dx * dx + dy * dy) + 0.1;
                    // Check if edge is a periphery edge
                    let isPeripheryEdge = false;
                    for (let i = 0; i < this.graph.periphery.length; i++) {
                        const idA = this.graph.periphery[i];
                        const idB = this.graph.periphery[(i + 1) % this.graph.periphery.length];
                        if ((edge.v1 === idA && edge.v2 === idB) || (edge.v1 === idB && edge.v2 === idA)) {
                            isPeripheryEdge = true;
                            break;
                        }
                    }
                    let idealLength = isPeripheryEdge ? kPeriphery : k;
                    // Spring force: proportional to (current length - ideal length)
                    let springConstant = 2;
                    // let springConstant = 5;
                    let springContantForPeriphery = 1;
                    // let springContantForPeriphery = 1;
                    if (isPeripheryEdge) {
                        springConstant = springContantForPeriphery;
                    }
                    let force = springConstant * (dist - idealLength);
                    let fx = force * dx / dist;
                    let fy = force * dy / dist;
                    v1.vx += fx;
                    v1.vy += fy;
                    v2.vx -= fx;
                    v2.vy -= fy;
                }
                // --- Angle resistance for periphery wedges ---
                const periphery = this.graph.periphery;
                for (let i = 0; i < periphery.length; i++) {
                    const prevId = periphery[(i - 1 + periphery.length) % periphery.length];
                    const currId = periphery[i];
                    const nextId = periphery[(i + 1) % periphery.length];
                    const prev = this.graph.getVertex(prevId);
                    const curr = this.graph.getVertex(currId);
                    const next = this.graph.getVertex(nextId);
                    if (!prev || !curr || !next) continue;
                    // Vectors from curr to prev and next
                    const v1x = prev.x - curr.x;
                    const v1y = prev.y - curr.y;
                    const v2x = next.x - curr.x;
                    const v2y = next.y - curr.y;
                    // Normalize
                    const len1 = Math.sqrt(v1x * v1x + v1y * v1y) || 1;
                    const len2 = Math.sqrt(v2x * v2x + v2y * v2y) || 1;
                    const u1x = v1x / len1;
                    const u1y = v1y / len1;
                    const u2x = v2x / len2;
                    const u2y = v2y / len2;
                    // Angle between vectors
                    let dot = u1x * u2x + u1y * u2y;
                    dot = Math.max(-1, Math.min(1, dot)); // Clamp for safety
                    const angle = Math.acos(dot);
                    // Resistance force proportional to deviation from ideal angle
                    const angleError = angle - idealAngle;
                    const forceMag = angleResistanceStrength * angleError;
                    // Perpendicular directions for prev and next
                    const perp1x = -u1y, perp1y = u1x;
                    const perp2x = u2y, perp2y = -u2x;
                    prev.vx += forceMag * perp1x;
                    prev.vy += forceMag * perp1y;
                    next.vx += forceMag * perp2x;
                    next.vy += forceMag * perp2y;
                    curr.vx -= forceMag * (perp1x + perp2x) / 2;
                    curr.vy -= forceMag * (perp1y + perp2y) / 2;
                }
                // --- Containment force: repulsive force from periphery vertices ---
                const peripheryVertices = this.graph.periphery.map(id => this.graph.getVertex(id));
                for (const v of vertices) {
                    let totalFx = 0, totalFy = 0;
                    for (const pv of peripheryVertices) {
                        if (v === pv) continue; // skip self
                        let dx = v.x - pv.x;
                        let dy = v.y - pv.y;
                        let dist = Math.sqrt(dx * dx + dy * dy) + 0.1;
                        // Repulsion strength increases as vertex gets closer to periphery
                        let force = containmentStrength * (1 / dist) * 2; // 2x multiplier for extra effect
                        totalFx += force * dx / dist;
                        totalFy += force * dy / dist;
                    }
                    v.vx += totalFx;
                    v.vy += totalFy;
                }
                // Update positions
                for (const v of vertices) {
                    // Freeze inner vertices: only update periphery
                    if (this.graph.periphery.includes(v.id)) {
                        v.x += v.vx * 0.02;
                        v.y += v.vy * 0.02;
                        v.vx *= damping;
                        v.vy *= damping;
                    } else {
                        // Check if vertex is inside the polygon
                        if (this.isPointInsidePolygon(v.x, v.y)) {
                            // Freeze vertices that are inside the polygon
                            v.vx = 0;
                            v.vy = 0;
                        } else {
                            // For vertices outside the polygon, apply normal physics
                            v.x += v.vx * 0.02;
                            v.y += v.vy * 0.02;
                            v.vx *= damping;
                            v.vy *= damping;
                        }
                    }
                }
            }

            redrawPeriphery() {
                // Redistribute periphery vertices evenly around a circle, keeping center and scale
                if (this.graph.periphery.length < 3) return;
                const visibleVertices = this.graph.periphery.map(id => this.graph.getVertex(id));
                // Find centroid
                let cx = 0, cy = 0;
                for (const v of visibleVertices) { cx += v.x; cy += v.y; }
                cx /= visibleVertices.length; cy /= visibleVertices.length;
                // Find average radius
                let avgR = 0;
                for (const v of visibleVertices) {
                    avgR += Math.sqrt((v.x - cx) ** 2 + (v.y - cy) ** 2);
                }
                avgR = avgR / visibleVertices.length * 1.1;
                // Place periphery vertices evenly
                for (let i = 0; i < visibleVertices.length; i++) {
                    const angle = (2 * Math.PI * i) / visibleVertices.length;
                    visibleVertices[i].x = cx + avgR * Math.cos(angle);
                    visibleVertices[i].y = cy + avgR * Math.sin(angle);
                }
                this.draw();
                updateStatus('Periphery redrawn and homogenized');
            }

            draw() {
                if (this.isForceDirected && this.graph.vertices.length > 0) {
                    if (!this.simulationRunning) this.startForceDirectedSimulation();
                }
                this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
                this.ctx.save();
                this.ctx.translate(this.offsetX, this.offsetY);
                this.ctx.scale(this.zoom, this.zoom);

                // Draw edges
                this.drawEdges();


                // Highlight periphery
                this.highlightPeriphery();

                // Draw vertices
                this.drawVertices();

                this.ctx.restore();
            }

            drawEdges() {
                // Draw edges as quadratic Bézier curves for compaction
                this.ctx.lineWidth = 6;
                for (const edge of this.graph.edges) {
                    const v1 = this.graph.getVertex(edge.v1);
                    const v2 = this.graph.getVertex(edge.v2);
                    if (!v1 || !v2 || v1.id > this.maxVisibleVertex || v2.id > this.maxVisibleVertex) continue;
                    // Control point: midpoint offset orthogonal to v1-v2
                    const mx = (v1.x + v2.x) / 2;
                    const my = (v1.y + v2.y) / 2;
                    const dx = v2.y - v1.y;
                    const dy = v1.x - v2.x;
                    const norm = Math.sqrt(dx * dx + dy * dy) || 1;
                    const offset = 0.18 * Math.sqrt((v1.x - v2.x) ** 2 + (v1.y - v2.y) ** 2);
                    const cx = mx + (dx / norm) * offset;
                    const cy = my + (dy / norm) * offset;
                    this.ctx.strokeStyle = '#666';
                    this.ctx.beginPath();
                    this.ctx.moveTo(v1.x, v1.y);
                    this.ctx.lineTo(v2.x, v2.y);
                    // this.ctx.quadraticCurveTo(cx, cy, v2.x, v2.y);
                    this.ctx.stroke();
                }
            }

            drawVertices() {
                for (const vertex of this.graph.vertices) {
                    if (vertex.id > this.maxVisibleVertex) continue;
                    // Draw vertex circle
                    this.ctx.fillStyle = this.graph.colors[vertex.color];
                    this.ctx.strokeStyle = '#333';
                    this.ctx.lineWidth = 2;
                    this.ctx.beginPath();
                    this.ctx.arc(vertex.x, vertex.y, vertex.radius, 0, 2 * Math.PI);
                    this.ctx.fill();
                    this.ctx.stroke();
                    // Draw vertex label (always show vertex number)
                    this.ctx.fillStyle = '#333';
                    this.ctx.font = `${Math.max(12, vertex.radius * 0.8)}px Arial`;
                    this.ctx.textAlign = 'center';
                    this.ctx.textBaseline = 'middle';
                    this.ctx.fillText(vertex.id.toString(), vertex.x, vertex.y);


                    // this.ctx.fillText(vertex.id.toString(), vertex.x, vertex.y - 8);

                    // // Draw position coordinates
                    // const posText = `(${Math.round(vertex.x)},${Math.round(vertex.y)})`;
                    // this.ctx.font = `${Math.max(8, vertex.radius * 0.5)}px Arial`;
                    // this.ctx.fillText(posText, vertex.x, vertex.y + 8);
                }
            }

            highlightPeriphery() {
                this.ctx.strokeStyle = '#FFD700';
                this.ctx.lineWidth = 4;

                for (let i = 0; i < this.graph.periphery.length; i++) {
                    const currentId = this.graph.periphery[i];
                    const nextId = this.graph.periphery[(i + 1) % this.graph.periphery.length];

                    const current = this.graph.getVertex(currentId);
                    const next = this.graph.getVertex(nextId);

                    if (!current || !next || current.id > this.maxVisibleVertex || next.id > this.maxVisibleVertex) continue;

                    this.ctx.beginPath();
                    this.ctx.moveTo(current.x, current.y);
                    this.ctx.lineTo(next.x, next.y);
                    this.ctx.stroke();
                }
            }

            centerAndResize() {
                if (this.graph.vertices.length === 0) return;

                const visibleVertices = this.graph.vertices.filter(v => v.id <= this.maxVisibleVertex);
                if (visibleVertices.length === 0) return;

                const minX = Math.min(...visibleVertices.map(v => v.x));
                const maxX = Math.max(...visibleVertices.map(v => v.x));
                const minY = Math.min(...visibleVertices.map(v => v.y));
                const maxY = Math.max(...visibleVertices.map(v => v.y));

                let graphWidth = maxX - minX;
                let graphHeight = maxY - minY;
                // Prevent zero width/height (which would cause zoom to be Infinity or NaN)
                if (graphWidth < 1) graphWidth = 1;
                if (graphHeight < 1) graphHeight = 1;

                const padding = 100;
                const availableWidth = this.canvas.width - padding * 2;
                const availableHeight = this.canvas.height - padding * 2;

                const scaleX = availableWidth / graphWidth;
                const scaleY = availableHeight / graphHeight;
                this.zoom = Math.min(scaleX, scaleY, 2);
                if (!isFinite(this.zoom) || this.zoom <= 0) this.zoom = 1;

                const centerX = (minX + maxX) / 2;
                const centerY = (minY + maxY) / 2;

                this.offsetX = this.canvas.width / 2 - centerX * this.zoom;
                this.offsetY = this.canvas.height / 2 - centerY * this.zoom;

                this.draw();
            }

            updateStats() {
                const visibleVertices = this.graph.vertices.filter(v => v.id <= this.maxVisibleVertex);
                const visibleEdges = this.graph.edges.filter(e => {
                    const v1 = this.graph.getVertex(e.v1);
                    const v2 = this.graph.getVertex(e.v2);
                    return v1 && v2 && v1.id <= this.maxVisibleVertex && v2.id <= this.maxVisibleVertex;
                });
                // Output all state variables for debugging
                const stateVars = {
                    zoom: this.zoom,
                    offsetX: this.offsetX,
                    offsetY: this.offsetY,
                    isDragging: this.isDragging,
                    lastMouseX: this.lastMouseX,
                    lastMouseY: this.lastMouseY,
                    maxVisibleVertex: this.maxVisibleVertex,
                    isForceDirected: this.isForceDirected,
                    simulationSteps: this.simulationSteps,
                    simulationMaxSteps: this.simulationMaxSteps,
                    simulationRunning: this.simulationRunning,
                    simulationInterval: !!this.simulationInterval,
                    periphery: this.graph.periphery.slice(),
                    selectedColor: this.graph.selectedColor,
                    displayShowIndex: this.graph.displayShowIndex,
                    colors: this.graph.colors,
                    verticesCount: this.graph.vertices.length,
                    edgesCount: this.graph.edges.length,
                    edges: this.graph.edges.map(e => ({ v1: e.v1, v2: e.v2 })) // Add edges info
                };
                let stateHtml = '<hr><b>Debug State:</b><br><pre style="font-size:11px;white-space:pre-wrap;">' + JSON.stringify(stateVars, null, 2) + '</pre>';
                document.getElementById('stats').innerHTML = `
            Vertices: ${visibleVertices.length}<br>
            Edges: ${visibleEdges.length}<br>
            Periphery: ${this.graph.periphery.length}
        ` + stateHtml;
            }

            toggleForceDirected() {
                this.isForceDirected = !this.isForceDirected;
                if (this.isForceDirected) {
                    this.startForceDirectedSimulation();
                    updateStatus('Force-directed layout enabled');
                } else {
                    updateStatus('Force-directed layout disabled');
                }
            }

            /**
             * Finds a random position outside the periphery polygon using bounding box, centroid, and random polar coordinates.
             * @param {number} [k1=0.75] - Minimum radius factor.
             * @param {number} [k2=1.5] - Maximum radius factor.
             * @returns {{x: number, y: number}} - Valid random position outside the polygon.
             */
            findRandomVertexPositionOutsidePolygon(vpId, vqId, k1 = 0.5, k2 = 2.5) {
                console.log('Finding random vertex position outside polygon...');
                const peripheryVertices = this.graph.periphery.map(id => this.graph.getVertex(id));
                if (peripheryVertices.length < 3) return null;

                // Find bounding box
                let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
                for (const v of peripheryVertices) {
                    if (v.x < xmin) xmin = v.x;
                    if (v.x > xmax) xmax = v.x;
                    if (v.y < ymin) ymin = v.y;
                    if (v.y > ymax) ymax = v.y;
                }

                // Calculate centroid
                let Cx = 0, Cy = 0;
                for (const v of peripheryVertices) {
                    Cx += v.x;
                    Cy += v.y;
                }
                Cx /= peripheryVertices.length;
                Cy /= peripheryVertices.length;

                console.log('bounding box:', { xmin, xmax, ymin, ymax }, 'centroid:', { Cx, Cy });

                // Max dimension
                const D = Math.max(xmax - xmin, ymax - ymin);
                const rmin = k1 * D;
                const rmax = k2 * D;

                // restrict random point between anglse made by vp and vq
                const vp = this.graph.getVertex(vpId);
                const vq = this.graph.getVertex(vqId);
                if (!vp || !vq) {
                    console.error(`Invalid vertex IDs: ${vpId}, ${vqId}`);
                    return null;
                }
                // Calculate angle range between vp and vq
                const angleVp = Math.atan2(vp.y - Cy, vp.x - Cx);
                const angleVq = Math.atan2(vq.y - Cy, vq.x - Cx);
                let angleMin = Math.min(angleVp, angleVq);
                let angleMax = Math.max(angleVp, angleVq);

                console.log('angle range:', { angleMin, angleMax });
                console.log('radius range:', { rmin, rmax });
                let pos = null;

                // vP, vQ are the two vertices between which we want to add a new vertex
                // find the shortest periphery.
                // vP, vP+1 = vector1
                // vQ, vQ-1 = vector2
                // find the angle by,
                // perpendicularVector1 = perpendicular vector1 towards the inside of the polygon
                // perpendicularVector2 = perpendicular vector2 towards the inside of the polygon
                // normPerp1 = normalize perpendicularVector1
                // normPerp2 = normalize perpendicularVector2
                // Cx, Cy = centroid of the polygon
                // rmin, rmax = radius range
                // find the angle between the two perpendicular vectors

                const vector1 = {
                    x: vp.x - Cx,
                    y: vp.y - Cy
                };
                const vector2 = {
                    x: vq.x - Cx,
                    y: vq.y - Cy
                };

                // Normalize the vectors
                const len1 = Math.sqrt(vector1.x * vector1.x + vector1.y * vector1.y) || 1;
                const len2 = Math.sqrt(vector2.x * vector2.x + vector2.y * vector2.y) || 1;
                const normPerp1 = {
                    x: -vector1.y / len1,
                    y: vector1.x / len1
                };
                const normPerp2 = {
                    x: -vector2.y / len2,
                    y: vector2.x / len2
                };

                // Calculate the angle between the two normalized perpendicular vectors
                const dotProduct = normPerp1.x * normPerp2.x + normPerp1.y * normPerp2.y;
                const angleBetween = Math.acos(Math.max(-1, Math.min(1, dotProduct))); // Clamp to [-1, 1]
                console.log('angle between perpendicular vectors:', angleBetween, 'dot product:', dotProduct);


                // Use shortest path between vp and vq (inclusive) from getVerticesBetweenBothDirections
                const { shortest, longest } = this.graph.getVerticesBetweenBothDirections(vpId, vqId);
                for (let attempt = 0; attempt < 100; attempt++) {
                    // const theta = angleMin + Math.random() * (angleMax - angleMin);
                    // Calculate the angle in the middle of the two vertices
                    // Calculate the angle in the middle of the two vertices
                    let theta = angleBetween;
                    // let theta = (angleMin + angleMax) / 2;
                    // Add a small random variation to avoid always placing at the exact middle
                    theta += (Math.random() - 0.5) * (angleMax - angleMin) * 0.3;

                    // Use larger radius values to ensure points are further from centroid
                    const r = rmin + Math.random() * (rmax - rmin); // Random radius between rmin and rmax

                    console.log('random angle:', theta, Math.cos(theta), Math.sin(theta), 'radius:', r);
                    const Px = Cx + r * Math.cos(theta) * Math.pow(-1, attempt); // Add some random scaling
                    const Py = Cy + r * Math.sin(theta) * Math.pow(-1, attempt);
                    if (!this.isPointInsidePolygon(Px, Py)) {
                        if (!this.wouldCreateCrossingEdges(Px, Py, shortest)) {
                            // Check if any touching vertex would exceed 4 edges
                            if (!this.wouldExceedMaxEdges(Px, Py, shortest, 5)) {
                                console.log('valid position found:', Px, Py);
                                pos = { x: Px, y: Py };
                                break; // Found a valid position
                            } else {
                                console.log('position would exceed max edges for a touching vertex, trying again', Px, Py);
                            }
                        } else {
                            console.log(attempt, 'position would create crossing edges, trying again', Px, Py);
                        }
                    }
                }
                // return { x: Cx + rmax, y: Cy };
                if (!pos) {
                    updateStatus('Could not find a valid position for random vertex');
                    return;
                }

                const newVertex = graph.addVertex(pos.x, pos.y, graph.selectedColor, shortest);
                // Add new vertex
                // const newVertex = this.graph.addVertex(newX, newY, this.graph.selectedColor, shortest);

                // Connect new vertex to all vertices in between (including vp and vq)
                for (const vId of shortest) {
                    this.graph.addEdge(newVertex.id, vId);
                }

                // Update periphery: remove verticesBetween, insert newVertex after vp
                // Update periphery: use all vertices from the longest path, then add the new vertex at the end
                // const { longest, shortest } = this.graph.getVerticesBetweenBothDirections(vpId, vqId);
                let newPeriphery = [...longest, newVertex.id];
                console.log('longest path:', longest, 'shortest path:', shortest, 'new vertex:', newVertex.id, 'new periphery:', newPeriphery);
                this.graph.periphery = newPeriphery;

                // Start force-directed rendering
                this.startForceDirectedSimulation();

                this.updateStats();
                this.draw();
            }
        }

        // Global variables
        let graph = new PlanarGraph();
        let renderer;
        let currentMode = 'ready';
        let selectedVp = null;
        let selectedVq = null;
        let autoAddInterval = null;

        // Initialize
        window.onload = function () {
            const canvas = document.getElementById('canvas');
            renderer = new GraphRenderer(canvas, graph);
            updateStatus('Application ready. Click "S - Start Triangle" to begin.');
        };

        // Command functions
        function startBasicGraph() {
            graph = new PlanarGraph();
            renderer.graph = graph;
            renderer.maxVisibleVertex = Infinity;

            // Create triangle vertices
            const centerX = renderer.canvas.width / 2;
            const centerY = renderer.canvas.height / 2;
            const radius = 100;

            const v1 = graph.addVertex(centerX, centerY - radius, 1);
            const v2 = graph.addVertex(centerX - radius * 0.866, centerY + radius * 0.5, 1);
            const v3 = graph.addVertex(centerX + radius * 0.866, centerY + radius * 0.5, 1);

            // Add edges
            graph.addEdge(v1.id, v2.id);
            graph.addEdge(v2.id, v3.id);
            graph.addEdge(v3.id, v1.id);

            // Set periphery
            graph.periphery = [v1.id, v2.id, v3.id];

            renderer.draw();
            renderer.updateStats();
            updateStatus('Basic triangle created');
        }

        function addRandomVertex() {
            if (graph.periphery.length < 3) {
                updateStatus('Need at least 3 vertices in periphery');
                return;
            }
            // Select two random adjacent vertices in periphery
            const randomIndex = Math.floor(Math.random() * graph.periphery.length);
            const randomIndex2 = Math.floor(Math.random() * graph.periphery.length);
            const vpId = graph.periphery[randomIndex];
            console.log(`Selected random index ${randomIndex}, ${randomIndex2}}`);
            // const vqId = graph.periphery[(randomIndex + 1) % graph.periphery.length];
            const vqId = graph.periphery[randomIndex2];
            if (vpId === vqId) {
                updateStatus('Selected same vertex for Vp and Vq, trying again');
                return;
            } else {
                console.log(`Selected Vp: ${vpId}, Vq: ${vqId}`);
            }

            // Use the new method to find a random position outside the polygon
            renderer.findRandomVertexPositionOutsidePolygon(vpId, vqId);
            renderer.draw();
            updateStatus('Random vertex added outside periphery');
        }

        function enterAddVertexMode() {
            if (graph.periphery.length < 2) {
                updateStatus('Need at least 2 vertices in periphery');
                return;
            }

            currentMode = 'selectVp';
            selectedVp = null;
            selectedVq = null;
            updateStatus('Select Vp (first vertex) from the periphery');
            updateModeIndicator('Select Vp');
        }

        function goToVertex() {
            const input = document.getElementById('goToVertex');
            const vertexNum = parseInt(input.value);

            if (isNaN(vertexNum) || vertexNum < 1) {
                updateStatus('Invalid vertex number');
                return;
            }

            renderer.maxVisibleVertex = vertexNum;
            renderer.draw();
            renderer.updateStats();
            updateStatus(`Showing vertices up to ${vertexNum}`);
        }

        function zoomIn() {
            renderer.zoom *= 1.2;
            renderer.draw();
        }

        function zoomOut() {
            renderer.zoom *= 0.8;
            renderer.draw();
        }

        function centerAndResize() {
            renderer.centerAndResize();
            updateStatus('Graph centered and resized');
        }

        function toggleDisplay() {
            graph.displayShowIndex = !graph.displayShowIndex;
            renderer.draw();
            updateStatus(graph.displayShowIndex ? 'Showing vertex indices' : 'Showing colors only');
        }

        function selectColor(colorNum) {
            graph.selectedColor = colorNum;

            // Update UI
            document.querySelectorAll('.color-box').forEach(box => {
                box.classList.remove('selected');
            });
            document.querySelector(`[data-color="${colorNum}"]`).classList.add('selected');

            updateStatus(`Selected color ${colorNum}`);
        }

        function updateStatus(message) {
            document.getElementById('status').textContent = message;
        }

        function updateModeIndicator(mode) {
            document.getElementById('modeIndicator').textContent = mode;
        }

        // Add Redraw command to UI
        function redrawPeriphery() {
            renderer.redrawPeriphery();
        }

        // Add color palette setup
        function setColorPalette(newPalette) {
            // newPalette: {1: '#hex', 2: '#hex', ...}
            Object.assign(graph.colors, newPalette);
            renderer.draw();
            updateStatus('Color palette updated');
        }

        function startAutoAdd() {
            if (!autoAddInterval) {
                autoAddInterval = setInterval(addRandomVertex, 3000);
                document.getElementById('startAutoAdd').disabled = true;
                document.getElementById('stopAutoAdd').disabled = false;
            }
        }

        function stopAutoAdd() {
            if (autoAddInterval) {
                clearInterval(autoAddInterval);
                autoAddInterval = null;
                document.getElementById('startAutoAdd').disabled = false;
                document.getElementById('stopAutoAdd').disabled = true;
            }
        }

        // document.addEventListener('DOMContentLoaded', function() {
        //     setTimeout(() => {
        //         startBasicGraph();
        //         // setInterval(addRandomVertex, 1000); // Commented out to use manual/auto controls
        //     }, 5000);
        // });