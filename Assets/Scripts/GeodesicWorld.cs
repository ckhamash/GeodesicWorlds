using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class GeodesicWorld : MonoBehaviour {

	public Vertex vertPrefab;
	//public Surface surface; // use for viewing in editor
	public List<Vertex> worldTiles;
	public int hexA, hexB; // defines the geodesic type
	public float radius;
	private float _edgeErrorValue = .001f;
	// Use this for initialization
	void Start() {
		GeodesicGenerator.GenerateGeodesic(hexA, hexB, vertPrefab);
		return;

		transform.position = Vector3.zero;
		transform.eulerAngles = Vector3.zero;

		// Initialize corner vertices (positions & neighbors)
		List<Vertex> worldCorners = InitializeIcosahedron(); // These define the corners of the icosahedron
		radius = worldCorners[0].transform.position.magnitude;

		// Create vertices for a face (bounding rotated positions in a triangle)
		List<Vertex> faceV, edgeV, cornerV; // temp face
		InitializeFaceVertices(out faceV, out edgeV, out cornerV);


		// get orientation for each face: 20 sets of 3 corners
		// for each worldCorner get each neighbor and its clockwise neighbor
		List<Vector3> cornerVec = new List<Vector3>();
		cornerVec.Add(cornerV[0].transform.position);
		cornerVec.Add(cornerV[1].transform.position);
		cornerVec.Add(cornerV[2].transform.position);

		List<Vertex> uniqueEdges = new List<Vertex>();
		List<Vertex> uniqueCorners = new List<Vertex>();
		List<Vertex> connectedEdges = new List<Vertex>(); // list of most recent completely connected edges (used for mending the seams)
		// Copy and place each face and remove all resulting duplicates
		for (int i = 0; i < worldCorners.Count; i++) {
			for (int j = 0; j < worldCorners[i].neighbors.Count; j++) {
				// if index < i then face is already complete
				List<int> worldCornerIndex = new List<int>();
				int neighbor = worldCorners.IndexOf(worldCorners[i].neighbors[j]);
				int adjacentNeighbor = worldCorners.IndexOf(worldCorners[i].neighbors[(j + 1) % worldCorners[i].neighbors.Count]);
				worldCornerIndex.Add(i);
				worldCornerIndex.Add(neighbor);
				worldCornerIndex.Add(adjacentNeighbor);
				if (neighbor > i && adjacentNeighbor > i) {
					// Place and orient the faces
					List<Vector3> worldVec = new List<Vector3>();
					worldVec.Add(worldCorners[i].transform.position);
					worldVec.Add(worldCorners[neighbor].transform.position);
					worldVec.Add(worldCorners[adjacentNeighbor].transform.position);

					List<Vertex> newFace = faceV;
					List<Vertex> newEdges = edgeV;
					List<Vertex> newCorners = cornerV;
					CopyFace(ref newFace, ref newEdges, ref newCorners);
					PlaceFace(newFace, cornerVec, worldVec);
					print("Place Face: " + worldCornerIndex[0] + ", " + worldCornerIndex[1] + ", " + worldCornerIndex[2]);

					// order newCorners to match worldCornerIndex order
					List<Vertex> tempNewCorners = new List<Vertex>();
					for (int m = 0; m < worldCornerIndex.Count; m++) {
						for (int n = 0; n < newCorners.Count; n++) {
							if (Vector3.Distance(newCorners[n].transform.position, worldCorners[worldCornerIndex[m]].transform.position) < _edgeErrorValue) {
								tempNewCorners.Add(newCorners[n]);
								break;
							}
						}
					}
					newCorners = tempNewCorners;

					List<Vertex> toDestroy = new List<Vertex>(); // these will be removed from newFace
					// check duplicate edges

					List<Vertex> toAdd = new List<Vertex>();
					foreach (Vertex n in newEdges) {
						bool combined = false;
						foreach (Vertex c in uniqueEdges) {
							if (Vector3.Distance(n.transform.position, c.transform.position) < _edgeErrorValue) {
								// combine neighbors to already existing duplicate
								// check # of real neighbors (not null)
								if (hexA == 0 || hexB == 0) { // edge has four neighbors
									// get index first of new and index last of pre-existing (these are the same vertex)
									int indexFirstNew = 0, indexLastOld = 0;
									for (int k = 0; k < n.neighbors.Count; k++) { // get index of first neighbor (clockwise)
										if (n.neighbors[k] == null && n.neighbors[(k + 1) % n.neighbors.Count] != null) {
											indexFirstNew = (k + 1) % n.neighbors.Count;
											break;
										}
									}
									for (int k = 0; k < c.neighbors.Count; k++) { // get index of last neighbor (clockwise)
										if (c.neighbors[k] != null && c.neighbors[(k + 1) % c.neighbors.Count] == null) {
											indexLastOld = k;
											break;
										}
									}
									// add next three new neighbors if old are null
									for (int k = 1; k <= 3; k++) {
										if (c.neighbors[(indexLastOld + k) % c.neighbors.Count] == null) {
											Vertex newNeighbor = n.neighbors[(indexFirstNew + k) % n.neighbors.Count];
											c.neighbors[(indexLastOld + k) % c.neighbors.Count] = newNeighbor;
											newNeighbor.neighbors[newNeighbor.neighbors.IndexOf(n)] = c;
										}
									}
								}
								else { // edge has three neighbors
									// get index first of new and index first empty of pre-existing (these are the same vertex)
									int indexFirstNew = 0, indexFirstNullOld = 0;
									for (int k = 0; k < n.neighbors.Count; k++) { // get index of first neighbor (clockwise)
										if (n.neighbors[k] == null && n.neighbors[(k + 1) % n.neighbors.Count] != null) {
											indexFirstNew = (k + 1) % n.neighbors.Count;
											break;
										}
									}
									for (int k = 0; k < c.neighbors.Count; k++) { // get index of last neighbor (clockwise)
										if (c.neighbors[k] != null && c.neighbors[(k + 1) % c.neighbors.Count] == null) {
											indexFirstNullOld = (k + 1) % c.neighbors.Count;
											break;
										}
									}
									// add three clockwise neighbors
									for (int k = 0; k < 3; k++) {
										Vertex newNeighbor = n.neighbors[(indexFirstNew + k) % n.neighbors.Count];
										c.neighbors[(indexFirstNullOld + k) % c.neighbors.Count] = newNeighbor;
										newNeighbor.neighbors[newNeighbor.neighbors.IndexOf(n)] = c;
									}
								}

								// destroy the new duplicate and remove the old from the duplicate list
								combined = true;
								toDestroy.Add(n);
								uniqueEdges.Remove(c);
								connectedEdges.Add(c);
								break;
							}
						}
						if (!combined) // to be added to duplicates
							toAdd.Add(n);
					}
					foreach (Vertex v in toAdd)
						uniqueEdges.Add(v);

					// add new corners' neighbors to existing corners, then destroy
					// use respective worldCorners neighbors to determine placement of neighbor in list 
					// note: only works because neighbors are defined clockwise
					for (int n = 0; n < newCorners.Count; n++) {
						bool match = false;
						//  neighbor index of new corner should match the neighbor index of its respective world corner
						int indexNewNeighbor = worldCorners[worldCornerIndex[n]].neighbors.IndexOf(worldCorners[worldCornerIndex[(n + 1) % worldCornerIndex.Count]]);
						foreach (Vertex c in uniqueCorners) {
							if (Vector3.Distance(newCorners[n].transform.position, c.transform.position) < _edgeErrorValue) {
								// num neighbors dont matter, only add first neighbor
								// only if neighbor is not already defined
								if (c.neighbors[indexNewNeighbor] == null) {
									int indexFirst = 0;
									for (int k = 0; k < newCorners[n].neighbors.Count; k++) { // get index of first neighbor (clockwise)
										if (newCorners[n].neighbors[k] == null && newCorners[n].neighbors[(k + 1) % newCorners[n].neighbors.Count] != null) {
											indexFirst = (k + 1) % newCorners[n].neighbors.Count;
											break;
										}
									}

									// need to pair new neighbor
									newCorners[n].neighbors[indexFirst].neighbors[newCorners[n].neighbors[indexFirst].neighbors.IndexOf(newCorners[n])] = c;
									// put new corner's neighbor in c.neighbors[newNeighborIndex]
									c.neighbors[indexNewNeighbor] = newCorners[n].neighbors[indexFirst];
								}
								match = true;
								// destroy duplicate corners
								toDestroy.Add(newCorners[n]);
								break;
							}
						}
						// new corner does not match any existing corners
						if (!match) { 
							// need to adjust neighbors to match pattern above (make sure neighbor exists)
							int indexFirst = 0;
							for (int k = 0; k < newCorners[n].neighbors.Count; k++) { // get index of first neighbor (clockwise)
								if (newCorners[n].neighbors[k] == null && newCorners[n].neighbors[(k + 1) % newCorners[n].neighbors.Count] != null) {
									indexFirst = (k + 1) % newCorners[n].neighbors.Count;
									break;
								}
							}
							int indexDiff =  indexNewNeighbor - indexFirst;
							if (indexDiff < 0)
								indexDiff += newCorners[n].neighbors.Count;
							List<Vertex> offsetNeighbors = new List<Vertex>(newCorners[n].neighbors);
							for (int k = 0; k < newCorners[n].neighbors.Count; k++) {
								offsetNeighbors[(k + indexDiff) % newCorners[n].neighbors.Count] = newCorners[n].neighbors[k];
							}
							newCorners[n].neighbors = offsetNeighbors;

							uniqueCorners.Add(newCorners[n]);
						}
					}
					
					for (int k = toDestroy.Count - 1; k >= 0; k--) {
						newFace.Remove(toDestroy[k]);
						DestroyImmediate(toDestroy[k].gameObject);
					}
					worldTiles.AddRange(newFace);
				}
			}
		}
		// destroy template face
		for (int i = faceV.Count - 1; i >= 0; i--) {
			DestroyImmediate(faceV[i].gameObject);
		}

		// Connect neighbors at seams
		connectedEdges.AddRange(uniqueCorners);
		while (connectedEdges.Count > 0) {
			List<Vertex> toConnect = new List<Vertex>();
			List<Vertex> toRemove = new List<Vertex>();
			print(connectedEdges.Count);
			foreach (Vertex connected in connectedEdges) {
				for (int i = 0; i < connected.neighbors.Count; i++) {
					bool completed = true;
					if (connected.neighbors[i] != null && connected.neighbors[(i + 1) % connected.neighbors.Count] != null && !connected.neighbors[i].neighbors.Contains(connected.neighbors[(i + 1) % connected.neighbors.Count])) {
						completed = false;
						Vertex leftNeighbor = connected.neighbors[i];
						Vertex rightNeighbor = connected.neighbors[(i + 1) % connected.neighbors.Count];
						// index of left's missing neighbor
						int leftIndex = (leftNeighbor.neighbors.IndexOf(connected) + leftNeighbor.neighbors.Count - 1) % leftNeighbor.neighbors.Count;
						// index of right's missing neighbor
						int rightIndex = (rightNeighbor.neighbors.IndexOf(connected) + 1) % rightNeighbor.neighbors.Count;

						leftNeighbor.neighbors[leftIndex] = rightNeighbor;
						rightNeighbor.neighbors[rightIndex] = leftNeighbor;

						// add completed neighbors to list
						int numNullLeft = 0;
						for (int j = 0; j < leftNeighbor.neighbors.Count; j++) {
							if (leftNeighbor.neighbors[j] == null)
								numNullLeft++;
						}
						int numNullRight = 0;
						for (int j = 0; j < rightNeighbor.neighbors.Count; j++) {
							if (rightNeighbor.neighbors[j] == null)
								numNullRight++;
						}
						if (numNullLeft < 2 && numNullLeft <= numNullRight) {
							toConnect.Add(leftNeighbor);
						}
						else if (numNullRight < 2) {
							toConnect.Add(rightNeighbor);
						}
					}
					if (completed) // this vertex and its neighbors are properly connected
						toRemove.Add(connected);
				}
			}

			foreach (Vertex v in toRemove)
				connectedEdges.Remove(v);
			connectedEdges.AddRange(toConnect);
		}

		// Adjust positions to create sphere shape
		// create a weight index based on distance from radius
		List<float> weightIndex = new List<float>();
		foreach (Vertex v in worldTiles) {
			weightIndex.Add((radius / v.transform.position.magnitude) - 0.5f);
			v.transform.position *= radius / v.transform.position.magnitude;
		}

		// Centroid Smoothing
		// move vertices to centroid of 3d shell formed with neighbors' opposite neighbors
		print(Vector3.Distance(worldTiles[0].transform.position, worldTiles[0].neighbors[0].transform.position));
		for (int n = 0; n < 0; n++) {
			// adjust all simultaneously
			List<Vector3> newPositions = new List<Vector3>();
			foreach (Vertex v in worldTiles) {
				// ignore corners
				if (v.neighbors.Count < 6) {
					newPositions.Add(v.transform.position);
					continue;
				}

				List<Vector3> hexCorners = new List<Vector3>();
				List<float> hexWeights = new List<float>();
				for (int i = 0; i < v.neighbors.Count; i++) {
					/*if (v.neighbors[i].neighbors.Count < 6) {
						int index = (v.neighbors[i].neighbors.IndexOf(v) + 2);
						hexCorners.Add(
							(v.neighbors[i].neighbors[index % v.neighbors[i].neighbors.Count].transform.position
							+ v.neighbors[i].neighbors[(index + 1) % v.neighbors[i].neighbors.Count].transform.position)
							/ 2);
						//hexWeights.Add(weightIndex[worldTiles.IndexOf(v.neighbors[i])]);
					}
					else {
						int index = (v.neighbors[i].neighbors.IndexOf(v) + 3) % v.neighbors[i].neighbors.Count;
						hexCorners.Add(v.neighbors[i].neighbors[index].transform.position);
						//hexWeights.Add(weightIndex[worldTiles.IndexOf(v.neighbors[i].neighbors[index])]);
					}*/
					hexCorners.Add(v.neighbors[i].transform.position);
					foreach (Vertex neighbor in v.neighbors) {
						float weight = 0;
						foreach (Vertex neighborSecond in neighbor.neighbors) {
							float dist = Vector3.Distance(neighbor.transform.position, neighborSecond.transform.position);
							weight += dist;
						}
						weight /= neighbor.neighbors.Count; // avg dist of neighbor's neighbors
						weight *= weight * weight * weight;
						hexWeights.Add(weight);
					}
				}

				Vector3 numeratorSum = Vector3.zero;
				float denominatorSum = 0;
				for (int i = 0; i < hexCorners.Count; i++) {
					Vector3 a = v.transform.position;
					Vector3 b = hexCorners[i % hexCorners.Count];
					Vector3 c = hexCorners[(i + 1) % hexCorners.Count];
					Vector3 avg = (a + b + c) / 3;
					float doubleArea = Vector3.Magnitude(Vector3.Cross(b - a, c - a));
					doubleArea *= (hexWeights[i] + hexWeights[(i + 1) % hexWeights.Count]) / 2; // apply avg weight of vectors
					numeratorSum += avg * doubleArea;
					denominatorSum += doubleArea;
				}

				Vector3 centroid = numeratorSum / denominatorSum;
				centroid *= radius / centroid.magnitude; // position on sphere surface
				newPositions.Add(centroid);
			}
			for (int i = 0; i < worldTiles.Count; i++) {
				worldTiles[i].transform.position = newPositions[i];
			}
			print(Vector3.Distance(worldTiles[0].transform.position, worldTiles[0].neighbors[0].transform.position));
		}

		// Radius Weighted Smoothing
		// smooth positions for a more centered position between neighbors
		for (int n = 0; n < 0; n++) {
			// adjust all simultaneously
			List<Vector3> newPositions = new List<Vector3>();
			newPositions.Clear();
			foreach (Vertex v in worldTiles) {
				// ignore corners
				if (v.neighbors.Count < 6) {
					newPositions.Add(v.transform.position);
					continue;
				}

				Vector3 position = Vector3.zero;
				foreach (Vertex neighbor in v.neighbors) {
					position += neighbor.transform.position * weightIndex[worldTiles.IndexOf(neighbor)]; // weight index defined before creating sphere
				}

				position /= v.neighbors.Count;
				position *= radius / position.magnitude;
				newPositions.Add(position);
			}

			for (int i = 0; i < worldTiles.Count; i++) {
				worldTiles[i].transform.position = newPositions[i];
			}

			print(Vector3.Distance(worldTiles[0].transform.position, worldTiles[0].neighbors[0].transform.position));
		}

		// Neighbor Distance Smoothing
		// smooth positions for a more centered position between neighbors
		for (int n = 0; n < 100; n++) { // n < (A + B) * 2 recommended number of passes
			// adjust all simultaneously
			List<Vector3> newPositions = new List<Vector3>();
			newPositions.Clear();
			foreach (Vertex v in worldTiles) {
				// ignore corners
				if (v.neighbors.Count < 6) {
					newPositions.Add(v.transform.position);
					continue;
				}

				Vector3 position = Vector3.zero;
				foreach (Vertex neighbor in v.neighbors) {
					float weight = 0;
					foreach (Vertex neighborSecond in neighbor.neighbors) {
						float dist = Vector3.Distance(neighbor.transform.position, neighborSecond.transform.position);
						weight += dist;
					}
					weight /= neighbor.neighbors.Count; // avg dist of neighbor's neighbors
					weight *= weight * weight;
					position += neighbor.transform.position * weight;
				}

				position /= v.neighbors.Count;
				position *= radius / position.magnitude;
				newPositions.Add(position);
			}
			for (int i = 0; i < worldTiles.Count; i++) {
				worldTiles[i].transform.position = newPositions[i];
			}
			print(Vector3.Distance(worldTiles[0].transform.position, worldTiles[0].neighbors[0].transform.position));
		}

		for (int i = worldCorners.Count - 1; i >= 0; i--) {
			DestroyImmediate(worldCorners[i].gameObject);
		}
	}
	
	// Update is called once per frame
	void Update() {
	}

	// sets up cornerVerts with their initial positions and neighbors
	List<Vertex> InitializeIcosahedron() {
		List<Vertex> icoVerts = new List<Vertex>();
		for (int i = 0; i < 12; i++) {
			Vertex vert = Instantiate(vertPrefab);
            vert.transform.parent = transform;
			//vert.maxNeighbors = 5;
			icoVerts.Add(vert);
		}
		float goldenRatio = (1 + Mathf.Sqrt(5))/2;
		// set positions
		icoVerts[0].transform.position = new Vector3(0, 1, goldenRatio); // 1, 8, 4, 5, 11 
		icoVerts[1].transform.position = new Vector3(0, -1, goldenRatio); // 0, 11, 6, 7, 8
		icoVerts[2].transform.position = new Vector3(0, -1, -goldenRatio); // 3, 9, 7, 6, 10
		icoVerts[3].transform.position = new Vector3(0, 1, -goldenRatio); // 2, 10, 5, 4, 9

		icoVerts[4].transform.position = new Vector3(1, goldenRatio, 0); // 5, 0, 8, 9, 3
		icoVerts[5].transform.position = new Vector3(-1, goldenRatio, 0); // 4, 3, 10, 11, 0
		icoVerts[6].transform.position = new Vector3(-1, -goldenRatio, 0); // 7, 1, 11, 10, 2
		icoVerts[7].transform.position = new Vector3(1, -goldenRatio, 0); // 6, 2, 9, 8, 1

		icoVerts[8].transform.position = new Vector3(goldenRatio, 0, 1); // 9, 4, 0, 1, 7
		icoVerts[9].transform.position = new Vector3(goldenRatio, 0, -1); // 8, 7, 2, 3, 4
		icoVerts[10].transform.position = new Vector3(-goldenRatio, 0, -1); // 11, 5, 3, 2, 6
		icoVerts[11].transform.position = new Vector3(-goldenRatio, 0, 1); // 10, 6, 1, 0, 5
		
		// define neighbors clockwise
		icoVerts[0].neighbors.Add(icoVerts[1]); icoVerts[0].neighbors.Add(icoVerts[8]);	icoVerts[0].neighbors.Add(icoVerts[4]); icoVerts[0].neighbors.Add(icoVerts[5]); icoVerts[0].neighbors.Add(icoVerts[11]);
		icoVerts[1].neighbors.Add(icoVerts[0]); icoVerts[1].neighbors.Add(icoVerts[11]); icoVerts[1].neighbors.Add(icoVerts[6]); icoVerts[1].neighbors.Add(icoVerts[7]); icoVerts[1].neighbors.Add(icoVerts[8]);
		icoVerts[2].neighbors.Add(icoVerts[3]); icoVerts[2].neighbors.Add(icoVerts[9]);	icoVerts[2].neighbors.Add(icoVerts[7]); icoVerts[2].neighbors.Add(icoVerts[6]); icoVerts[2].neighbors.Add(icoVerts[10]);
		icoVerts[3].neighbors.Add(icoVerts[2]);	icoVerts[3].neighbors.Add(icoVerts[10]); icoVerts[3].neighbors.Add(icoVerts[5]); icoVerts[3].neighbors.Add(icoVerts[4]); icoVerts[3].neighbors.Add(icoVerts[9]);

		icoVerts[4].neighbors.Add(icoVerts[5]); icoVerts[4].neighbors.Add(icoVerts[0]); icoVerts[4].neighbors.Add(icoVerts[8]); icoVerts[4].neighbors.Add(icoVerts[9]); icoVerts[4].neighbors.Add(icoVerts[3]);
		icoVerts[5].neighbors.Add(icoVerts[4]); icoVerts[5].neighbors.Add(icoVerts[3]); icoVerts[5].neighbors.Add(icoVerts[10]); icoVerts[5].neighbors.Add(icoVerts[11]); icoVerts[5].neighbors.Add(icoVerts[0]);
		icoVerts[6].neighbors.Add(icoVerts[7]); icoVerts[6].neighbors.Add(icoVerts[1]); icoVerts[6].neighbors.Add(icoVerts[11]); icoVerts[6].neighbors.Add(icoVerts[10]); icoVerts[6].neighbors.Add(icoVerts[2]);
		icoVerts[7].neighbors.Add(icoVerts[6]); icoVerts[7].neighbors.Add(icoVerts[2]); icoVerts[7].neighbors.Add(icoVerts[9]); icoVerts[7].neighbors.Add(icoVerts[8]); icoVerts[7].neighbors.Add(icoVerts[1]);

		icoVerts[8].neighbors.Add(icoVerts[9]); icoVerts[8].neighbors.Add(icoVerts[4]); icoVerts[8].neighbors.Add(icoVerts[0]); icoVerts[8].neighbors.Add(icoVerts[1]); icoVerts[8].neighbors.Add(icoVerts[7]);
		icoVerts[9].neighbors.Add(icoVerts[8]); icoVerts[9].neighbors.Add(icoVerts[7]); icoVerts[9].neighbors.Add(icoVerts[2]); icoVerts[9].neighbors.Add(icoVerts[3]); icoVerts[9].neighbors.Add(icoVerts[4]);
		icoVerts[10].neighbors.Add(icoVerts[11]); icoVerts[10].neighbors.Add(icoVerts[5]); icoVerts[10].neighbors.Add(icoVerts[3]); icoVerts[10].neighbors.Add(icoVerts[2]); icoVerts[10].neighbors.Add(icoVerts[6]);
		icoVerts[11].neighbors.Add(icoVerts[10]); icoVerts[11].neighbors.Add(icoVerts[6]); icoVerts[11].neighbors.Add(icoVerts[1]); icoVerts[11].neighbors.Add(icoVerts[0]); icoVerts[11].neighbors.Add(icoVerts[5]);

		// set length of edge
		float length = new Vector2(hexB + hexA * 0.5f, hexA * Mathf.Sqrt(3) / 2).magnitude;
		foreach (Vertex vert in icoVerts) {
			vert.transform.position *= length / 2; // 2 is the initial side length
		}
		return icoVerts;
	}

	// creates, places, and sets neighbors for a single face, (then orients it)
	// returns list of hexes with edges and corners defined
	void InitializeFaceVertices(out List<Vertex> verts, out List<Vertex> edges, out List<Vertex> corners) {
		verts = new List<Vertex>();
		edges = new List<Vertex>();
		corners = new List<Vertex>();
		int l = hexA + hexB + 1;
		// Create hexagonal grid of l by l vertices (bottom left is 0,0; next row is 0.5,sqrt(3)/2)
		for (int i = 0; i < l; i++) {
			for (int j = 0; j < l; j++) {
				Vertex vert = Instantiate(vertPrefab);
				vert.transform.parent = transform;
				vert.transform.position = new Vector3(j + i * 0.5f, i * Mathf.Sqrt(3) / 2);
				//vert.maxNeighbors = 6;
				//for (int k = 0; k < vert.maxNeighbors; k++)
				for (int k = 0; k < 6; k++)
					vert.neighbors.Add(null); // initialize neighbor list
				verts.Add(vert);
			}
		}
		// define neighbors clockwise: 0 = i+1, 1 = i-length+1, 2 = i-length, 3 = i-1, 4 = i+length-1, 5 = i+length
		for (int i = 0; i < verts.Count; i++) {
			if (i % l != l - 1) // check right
				verts[i].neighbors[0] = verts[i + 1];
			if (i % l != l - 1 && i >= l) // check bottom right
				verts[i].neighbors[1] = verts[i - l + 1];
			if (i >= l) // check bottom left
				verts[i].neighbors[2] = verts[i - l];
			if (i % l != 0) // check left
				verts[i].neighbors[3] = verts[i - 1];
			if (i % l != 0 && i < verts.Count - l) // check top left
				verts[i].neighbors[4] = verts[i + l - 1];
			if (i < verts.Count - l) // check top right
				verts[i].neighbors[5] = verts[i + l];
		}

		// define the triangle corners: 0 = [b], 1 = [(a+b)length], 2 = [b+a + b*length]
		corners.Add(verts[hexB]);
		corners.Add(verts[l * (hexA + hexB)]);
		corners.Add(verts[(hexA + hexB) + (l * hexB)]);
		corners[0].GetComponent<Renderer>().material.color = Color.magenta;
		corners[1].GetComponent<Renderer>().material.color = Color.magenta;
		corners[2].GetComponent<Renderer>().material.color = Color.magenta;

		foreach (Vertex c in corners) { // reduce neighbors to 5
			//c.maxNeighbors = 5;
			for (int i = 0; i < c.neighbors.Count; i++) {
				if (c.neighbors[i] == null) {
					c.neighbors.RemoveAt(i);
					print("Remove");
					break;
				}
			}
		}

		// define points on edge  (these will be duplicates)
		int gcd = GreatestCommonDivisor(hexA, hexB);
		int a = hexA / gcd;
		int b = hexB / gcd;
		for (int i = 1; i < gcd; i++) {
			edges.Add(verts[(hexB) + ((b * i * (l - 1)) + (a * i * l))]); // corner[0] --> corner[1]
			edges.Add(verts[(l * (hexA + hexB)) + ((b * i) + (a * i * (-l + 1)))]); // corner[1] --> corner[2]
			edges.Add(verts[(hexB) + ((b * i * l) + (a * i))]); // corner[0] --> corner[2]
		}
		// debug
		for (int i = 0; i < edges.Count; i++)
			edges[i].GetComponent<Renderer>().material.color = Color.cyan;

		// define triangle bounds
		PolygonCollider2D triangleBounds = gameObject.AddComponent<PolygonCollider2D>();
		Vector2[] points = new Vector2[] { corners[0].transform.position, corners[1].transform.position, corners[2].transform.position };
		triangleBounds.SetPath(0, points);


		float leftX = corners[0].transform.position.x; // leftmost X val
		float rightX = corners[0].transform.position.x; // rightmost X val
		for (int i = 0; i < corners.Count; i++) {
			if (leftX > corners[i].transform.position.x)
				leftX = corners[i].transform.position.x;
			if (rightX < corners[i].transform.position.x)
				rightX = corners[i].transform.position.x;
		}

		// remove external neighbors (keep border positions), specify vertices that fall on triangle edge(these will be duplicates on the globe)
		List<Vertex> toRemove = new List<Vertex>();
		for (int i = 0; i < verts.Count; i++) {
			// trim obvious external (left and right)
			if (verts[i].transform.position.x < leftX || verts[i].transform.position.x > rightX) {
				toRemove.Add(verts[i]);
				continue;
			}

			if (edges.Contains(verts[i]) || corners.Contains(verts[i]))
				continue;

			if (!triangleBounds.OverlapPoint(verts[i].transform.position)) // outside of triangle
				toRemove.Add(verts[i]);
		}
		// delete triangle bounds
		Destroy(gameObject.GetComponent<PolygonCollider2D>());

		for (int i = toRemove.Count - 1; i >= 0; i--) {
			verts.Remove(toRemove[i]);
			DestroyImmediate(toRemove[i].gameObject);
		}
		//Debug
		corners[0].GetComponent<Renderer>().material.color = Color.magenta;
		corners[1].GetComponent<Renderer>().material.color = Color.magenta;
		corners[2].GetComponent<Renderer>().material.color = Color.magenta;
	}

	// returns true if point lies along a line segment
	bool PointOnSegment(Vector3 point, Vector3 start, Vector3 end) {
		//print(distanceToSegment(point, start, end));
		float d = DistanceToSegmentSquared(point, start, end);
		//print(d);
		return d <= _edgeErrorValue * _edgeErrorValue;
		/*if (point == start || point == end)
			return true;
		Vector3 lineA = start - point;
		Vector3 lineB = start - end;
		// error range of arctan(0.5/l)
		float l = lineB.magnitude;
		float angleError = Mathf.Rad2Deg * Mathf.Atan2(edgeErrorValue, l);
		if (lineA.magnitude <= lineB.magnitude && Mathf.Abs(Vector3.Angle(lineA, lineB)) < angleError) // on the line
			return true;
		else
			return false;*/
	}

	// returns the shortest distance from a point to a line segment
	float DistanceToSegmentSquared(Vector2 point, Vector2 start, Vector2 end) {
		float dot = (point.x - start.x)*(end.x - start.x) + (point.y - start.y)*(end.y - start.y);
		//float dot = Vector2.Dot(point - start, end - start);
		dot /= (end - start).sqrMagnitude; // normalized
		if (dot < 0) // point exists some distance before the start of the line
			return (point - start).sqrMagnitude;
		else if (dot <= 1) {
			return (point - start).sqrMagnitude - dot * dot * (end - start).sqrMagnitude;
		}
		else // point exists some distance after the end of the line
			return (point - end).sqrMagnitude;
	}

	// return the greatest common divisor of two ints
	int GreatestCommonDivisor(int a, int b) {
		if (b == 0)
			return a;
		else
			return GreatestCommonDivisor(b, a % b);
	}

	// given an assembled face of the geodesic, place it
	void PlaceFace(List<Vertex> faceVerts, List<Vector3> faceCornerLocations, List<Vector3> worldCornerLocations) {
		// normVec = crossproduct(a-b,a-c) where abc are defined clockwise(following LHR such that normal points outward)
		// upVec = b - mid(a, c)
		// use Quaternion.LookRotation

		// calculate copy orientation
		Vector3 norm = Vector3.Cross(faceCornerLocations[0] - faceCornerLocations[1], faceCornerLocations[0] - faceCornerLocations[2]);
		Vector3 up = faceCornerLocations[1] - (faceCornerLocations[0] + faceCornerLocations[2]) / 2;

		// create temp gameobj, give it the copy orientation and position of face corner a
		GameObject temp = new GameObject();
		temp.transform.rotation = Quaternion.LookRotation(norm, up);
		temp.transform.position = faceCornerLocations[0];

		// parent verts to temp object
		foreach (Vertex vert in faceVerts) {
			vert.transform.position = temp.transform.InverseTransformPoint(vert.transform.position);
			//vert.transform.SetParent(temp.transform);
		}

		// calculate correct orientation
		norm = Vector3.Cross(worldCornerLocations[0] - worldCornerLocations[1], worldCornerLocations[0] - worldCornerLocations[2]);
		up = worldCornerLocations[1] - (worldCornerLocations[0] + worldCornerLocations[2]) / 2;

		// manipulate temp obj to correct position and orientation
		temp.transform.rotation = Quaternion.LookRotation(norm, up);
		temp.transform.position = worldCornerLocations[0];

		// parent verts to world
		foreach (Vertex vert in faceVerts) {
			vert.transform.position = temp.transform.TransformPoint(vert.transform.position);
			vert.transform.SetParent(this.transform);
		}
		
		Destroy(temp);
	}

	// Duplicates a face including defined edges and corners with appropriate neighbors from the new list
	void CopyFace(ref List<Vertex> verts, ref List<Vertex> edges, ref List<Vertex> corners) {
		List<Vertex> copyVerts = new List<Vertex>();
		List<Vertex> copyEdges = new List<Vertex>();
		List<Vertex> copyCorners = new List<Vertex>();

		// copy position and neighbor count
		for (int i = 0; i < verts.Count; i++) {
			Vertex copy = Instantiate(vertPrefab);
			copy.transform.position = verts[i].transform.position;
			//copy.maxNeighbors = verts[i].maxNeighbors;
			copy.GetComponent<Renderer>().material.color = verts[i].GetComponent<Renderer>().material.color; // debug
			copyVerts.Add(copy);
			if (edges.Contains(verts[i]))
				copyEdges.Add(copy);
			if (corners.Contains(verts[i]))
				copyCorners.Add(copy);
		}

		// get neighbors based on index in verts
		for (int i = 0; i < verts.Count; i++) {
			for (int j = 0; j < verts[i].neighbors.Count; j++) {
				int index = verts.IndexOf(verts[i].neighbors[j]);
				if (index < 0) // null value
					copyVerts[i].neighbors.Add(null);
				else
					copyVerts[i].neighbors.Add(copyVerts[index]);
			}
		}

		verts = copyVerts;
		edges = copyEdges;
		corners = copyCorners;
	}
}
