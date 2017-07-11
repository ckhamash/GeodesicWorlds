﻿using UnityEngine;
using System.Collections;
using System.Collections.Generic;

// class for generating Class I, II, & III geodesic spheres as a collection of Vertex objects

public static class GeodesicGenerator {

	private static int _faceLengthX = 1, _faceLengthY = 1; // number of vertices to traverse between "corners" of the sphere
	private static Vertex _vertexPrefab; // instantiate new Vertices with this prefab
	private static List<Vertex> _geodesicCorners; // the vertices that define the underlying geodesic shape.
	private static List<Vertex> _templateVertices, _templateEdges, _templateCorners; // temp face
	private static float _edgeErrorValue = .001f; // used for finding duplicates

	public static void SetDimensions(int x, int y) {
		if (x < 1 && y < 1) // div by 0 if both are 0
			y = 1;
		_faceLengthX = Mathf.Max(x, 0);
		_faceLengthY = Mathf.Max(y, 0);
	}

	public static int GetDimensionX() {
		return _faceLengthX;
	}

	public static int GetDimensionY() {
		return _faceLengthY;
	}

	// Generates a connected vertex mesh of a geodesic sphere based on a regular icosahedron
	public static void GenerateGeodesic(int x, int y, Vertex prefab) {
		SetDimensions(x, y);
		InitializePrefab(prefab);
		InitializeIcosahedron();
		InitializeFaceVertices();
		
		// destroy all leftover objects when finished (those in the static lists, they will not be kept)
	}

	private static void InitializePrefab(Vertex prefab) {
		_vertexPrefab = prefab;
	}

	// Create the initial vertices of the icosahedron
	private static void InitializeIcosahedron() {
		_geodesicCorners = new List<Vertex>();
		for (int i = 0; i < 12; i++) {
			Vertex vert = Object.Instantiate(_vertexPrefab);
			//vert.transform.parent = transform;
			//vert.maxNeighbors = 5;
			_geodesicCorners.Add(vert);
		}

		float goldenRatio = (1 + Mathf.Sqrt(5)) / 2;
		// set positions
		_geodesicCorners[0].transform.position = new Vector3(0, 1, goldenRatio); // 1, 8, 4, 5, 11 
		_geodesicCorners[1].transform.position = new Vector3(0, -1, goldenRatio); // 0, 11, 6, 7, 8
		_geodesicCorners[2].transform.position = new Vector3(0, -1, -goldenRatio); // 3, 9, 7, 6, 10
		_geodesicCorners[3].transform.position = new Vector3(0, 1, -goldenRatio); // 2, 10, 5, 4, 9

		_geodesicCorners[4].transform.position = new Vector3(1, goldenRatio, 0); // 5, 0, 8, 9, 3
		_geodesicCorners[5].transform.position = new Vector3(-1, goldenRatio, 0); // 4, 3, 10, 11, 0
		_geodesicCorners[6].transform.position = new Vector3(-1, -goldenRatio, 0); // 7, 1, 11, 10, 2
		_geodesicCorners[7].transform.position = new Vector3(1, -goldenRatio, 0); // 6, 2, 9, 8, 1

		_geodesicCorners[8].transform.position = new Vector3(goldenRatio, 0, 1); // 9, 4, 0, 1, 7
		_geodesicCorners[9].transform.position = new Vector3(goldenRatio, 0, -1); // 8, 7, 2, 3, 4
		_geodesicCorners[10].transform.position = new Vector3(-goldenRatio, 0, -1); // 11, 5, 3, 2, 6
		_geodesicCorners[11].transform.position = new Vector3(-goldenRatio, 0, 1); // 10, 6, 1, 0, 5

		// define neighbors clockwise
		_geodesicCorners[0].neighbors.Add(_geodesicCorners[1]); _geodesicCorners[0].neighbors.Add(_geodesicCorners[8]); _geodesicCorners[0].neighbors.Add(_geodesicCorners[4]); _geodesicCorners[0].neighbors.Add(_geodesicCorners[5]); _geodesicCorners[0].neighbors.Add(_geodesicCorners[11]);
		_geodesicCorners[1].neighbors.Add(_geodesicCorners[0]); _geodesicCorners[1].neighbors.Add(_geodesicCorners[11]); _geodesicCorners[1].neighbors.Add(_geodesicCorners[6]); _geodesicCorners[1].neighbors.Add(_geodesicCorners[7]); _geodesicCorners[1].neighbors.Add(_geodesicCorners[8]);
		_geodesicCorners[2].neighbors.Add(_geodesicCorners[3]); _geodesicCorners[2].neighbors.Add(_geodesicCorners[9]); _geodesicCorners[2].neighbors.Add(_geodesicCorners[7]); _geodesicCorners[2].neighbors.Add(_geodesicCorners[6]); _geodesicCorners[2].neighbors.Add(_geodesicCorners[10]);
		_geodesicCorners[3].neighbors.Add(_geodesicCorners[2]); _geodesicCorners[3].neighbors.Add(_geodesicCorners[10]); _geodesicCorners[3].neighbors.Add(_geodesicCorners[5]); _geodesicCorners[3].neighbors.Add(_geodesicCorners[4]); _geodesicCorners[3].neighbors.Add(_geodesicCorners[9]);

		_geodesicCorners[4].neighbors.Add(_geodesicCorners[5]); _geodesicCorners[4].neighbors.Add(_geodesicCorners[0]); _geodesicCorners[4].neighbors.Add(_geodesicCorners[8]); _geodesicCorners[4].neighbors.Add(_geodesicCorners[9]); _geodesicCorners[4].neighbors.Add(_geodesicCorners[3]);
		_geodesicCorners[5].neighbors.Add(_geodesicCorners[4]); _geodesicCorners[5].neighbors.Add(_geodesicCorners[3]); _geodesicCorners[5].neighbors.Add(_geodesicCorners[10]); _geodesicCorners[5].neighbors.Add(_geodesicCorners[11]); _geodesicCorners[5].neighbors.Add(_geodesicCorners[0]);
		_geodesicCorners[6].neighbors.Add(_geodesicCorners[7]); _geodesicCorners[6].neighbors.Add(_geodesicCorners[1]); _geodesicCorners[6].neighbors.Add(_geodesicCorners[11]); _geodesicCorners[6].neighbors.Add(_geodesicCorners[10]); _geodesicCorners[6].neighbors.Add(_geodesicCorners[2]);
		_geodesicCorners[7].neighbors.Add(_geodesicCorners[6]); _geodesicCorners[7].neighbors.Add(_geodesicCorners[2]); _geodesicCorners[7].neighbors.Add(_geodesicCorners[9]); _geodesicCorners[7].neighbors.Add(_geodesicCorners[8]); _geodesicCorners[7].neighbors.Add(_geodesicCorners[1]);

		_geodesicCorners[8].neighbors.Add(_geodesicCorners[9]); _geodesicCorners[8].neighbors.Add(_geodesicCorners[4]); _geodesicCorners[8].neighbors.Add(_geodesicCorners[0]); _geodesicCorners[8].neighbors.Add(_geodesicCorners[1]); _geodesicCorners[8].neighbors.Add(_geodesicCorners[7]);
		_geodesicCorners[9].neighbors.Add(_geodesicCorners[8]); _geodesicCorners[9].neighbors.Add(_geodesicCorners[7]); _geodesicCorners[9].neighbors.Add(_geodesicCorners[2]); _geodesicCorners[9].neighbors.Add(_geodesicCorners[3]); _geodesicCorners[9].neighbors.Add(_geodesicCorners[4]);
		_geodesicCorners[10].neighbors.Add(_geodesicCorners[11]); _geodesicCorners[10].neighbors.Add(_geodesicCorners[5]); _geodesicCorners[10].neighbors.Add(_geodesicCorners[3]); _geodesicCorners[10].neighbors.Add(_geodesicCorners[2]); _geodesicCorners[10].neighbors.Add(_geodesicCorners[6]);
		_geodesicCorners[11].neighbors.Add(_geodesicCorners[10]); _geodesicCorners[11].neighbors.Add(_geodesicCorners[6]); _geodesicCorners[11].neighbors.Add(_geodesicCorners[1]); _geodesicCorners[11].neighbors.Add(_geodesicCorners[0]); _geodesicCorners[11].neighbors.Add(_geodesicCorners[5]);

		// set length of edge
		float length = new Vector2(_faceLengthY + _faceLengthX * 0.5f, _faceLengthX * Mathf.Sqrt(3) / 2).magnitude;
		foreach (Vertex vert in _geodesicCorners) {
			vert.transform.position *= length / 2; // 2 is the initial side length
		}
	}

	// Create the template the faces of the geodesic will be created from
	private static void InitializeFaceVertices() {
		_templateVertices = new List<Vertex>();
		_templateEdges = new List<Vertex>();
		_templateCorners = new List<Vertex>();
		int l = _faceLengthX + _faceLengthY + 1;
		// Create hexagonal grid of l by l vertices (bottom left is 0,0; next row is 0.5,sqrt(3)/2)
		for (int i = 0; i < l; i++) {
			for (int j = 0; j < l; j++) {
				Vertex vert = Object.Instantiate(_vertexPrefab);
				//vert.transform.parent = transform;
				vert.transform.position = new Vector3(j + i * 0.5f, i * Mathf.Sqrt(3) / 2);
				//vert.maxNeighbors = 6;
				//for (int k = 0; k < vert.maxNeighbors; k++)
				for (int k = 0; k < 6; k++)
					vert.neighbors.Add(null); // initialize neighbor list
				_templateVertices.Add(vert);
			}
		}
		// define neighbors clockwise: 0 = i+1, 1 = i-length+1, 2 = i-length, 3 = i-1, 4 = i+length-1, 5 = i+length
		for (int i = 0; i < _templateVertices.Count; i++) {
			if (i % l != l - 1) // check right
				_templateVertices[i].neighbors[0] = _templateVertices[i + 1];
			if (i % l != l - 1 && i >= l) // check bottom right
				_templateVertices[i].neighbors[1] = _templateVertices[i - l + 1];
			if (i >= l) // check bottom left
				_templateVertices[i].neighbors[2] = _templateVertices[i - l];
			if (i % l != 0) // check left
				_templateVertices[i].neighbors[3] = _templateVertices[i - 1];
			if (i % l != 0 && i < _templateVertices.Count - l) // check top left
				_templateVertices[i].neighbors[4] = _templateVertices[i + l - 1];
			if (i < _templateVertices.Count - l) // check top right
				_templateVertices[i].neighbors[5] = _templateVertices[i + l];
		}

		// define the triangle corners: 0 = [b], 1 = [(a+b)length], 2 = [b+a + b*length]
		_templateCorners.Add(_templateVertices[_faceLengthY]);
		_templateCorners.Add(_templateVertices[l * (_faceLengthX + _faceLengthY)]);
		_templateCorners.Add(_templateVertices[(_faceLengthX + _faceLengthY) + (l * _faceLengthY)]);
		_templateCorners[0].GetComponent<Renderer>().material.color = Color.magenta;
		_templateCorners[1].GetComponent<Renderer>().material.color = Color.magenta;
		_templateCorners[2].GetComponent<Renderer>().material.color = Color.magenta;

		foreach (Vertex c in _templateCorners) { // reduce neighbors to 5
			//c.maxNeighbors = 5;
			for (int i = 0; i < c.neighbors.Count; i++) {
				if (c.neighbors[i] == null) {
					c.neighbors.RemoveAt(i);
					break;
				}
			}
		}

		// define points on edge  (these will be duplicates)
		int gcd = GreatestCommonDivisor(_faceLengthX, _faceLengthY);
		int a = _faceLengthX / gcd;
		int b = _faceLengthY / gcd;
		for (int i = 1; i < gcd; i++) {
			_templateEdges.Add(_templateVertices[(_faceLengthY) + ((b * i * (l - 1)) + (a * i * l))]); // corner[0] --> corner[1]
			_templateEdges.Add(_templateVertices[(l * (_faceLengthX + _faceLengthY)) + ((b * i) + (a * i * (-l + 1)))]); // corner[1] --> corner[2]
			_templateEdges.Add(_templateVertices[(_faceLengthY) + ((b * i * l) + (a * i))]); // corner[0] --> corner[2]
		}
		// debug
		for (int i = 0; i < _templateEdges.Count; i++)
			_templateEdges[i].GetComponent<Renderer>().material.color = Color.cyan;

		// define triangle bounds
		GameObject temp = new GameObject();
		PolygonCollider2D triangleBounds = temp.AddComponent<PolygonCollider2D>();
		Vector2[] points = new Vector2[] { _templateCorners[0].transform.position, _templateCorners[1].transform.position, _templateCorners[2].transform.position };
		triangleBounds.SetPath(0, points);


		float leftX = _templateCorners[0].transform.position.x; // leftmost X val
		float rightX = _templateCorners[0].transform.position.x; // rightmost X val
		for (int i = 0; i < _templateCorners.Count; i++) {
			if (leftX > _templateCorners[i].transform.position.x)
				leftX = _templateCorners[i].transform.position.x;
			if (rightX < _templateCorners[i].transform.position.x)
				rightX = _templateCorners[i].transform.position.x;
		}

		// remove external neighbors (keep border positions), specify vertices that fall on triangle edge(these will be duplicates on the globe)
		List<Vertex> toRemove = new List<Vertex>();
		for (int i = 0; i < _templateVertices.Count; i++) {
			// trim obvious external (left and right)
			if (_templateVertices[i].transform.position.x < leftX || _templateVertices[i].transform.position.x > rightX) {
				toRemove.Add(_templateVertices[i]);
				continue;
			}

			if (_templateEdges.Contains(_templateVertices[i]) || _templateCorners.Contains(_templateVertices[i]))
				continue;

			if (!triangleBounds.OverlapPoint(_templateVertices[i].transform.position)) // outside of triangle
				toRemove.Add(_templateVertices[i]);
		}
		// delete temp object
		Object.Destroy(temp);

		for (int i = 0; i <toRemove.Count; i++) {
			_templateVertices.Remove(toRemove[i]);
			Object.DestroyImmediate(toRemove[i].gameObject);
		}
		//Debug
		_templateCorners[0].GetComponent<Renderer>().material.color = Color.magenta;
		_templateCorners[1].GetComponent<Renderer>().material.color = Color.magenta;
		_templateCorners[2].GetComponent<Renderer>().material.color = Color.magenta;
	}

	// Return the greatest common divisor of two ints
	private static int GreatestCommonDivisor(int a, int b) {
		if (b == 0)
			return a;
		else
			return GreatestCommonDivisor(b, a % b);
	}

	// Duplicates a face including defined edges and corners with appropriate neighbors from the new list
	private static void CopyFace(ref List<Vertex> verts, ref List<Vertex> edges, ref List<Vertex> corners) {
		List<Vertex> copyVerts = new List<Vertex>();
		List<Vertex> copyEdges = new List<Vertex>();
		List<Vertex> copyCorners = new List<Vertex>();

		// copy position and neighbor count
		for (int i = 0; i < verts.Count; i++) {
			Vertex copy = Object.Instantiate(_vertexPrefab);
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

	// Given an assembled face of the geodesic, place it
	private static void MoveFace(List<Vertex> faceVerts, List<Vector3> faceCornerLocations, List<Vector3> worldCornerLocations) {
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

		// remove verts' parenting
		foreach (Vertex vert in faceVerts) {
			vert.transform.position = temp.transform.TransformPoint(vert.transform.position);
			vert.transform.SetParent(null);
		}

		Object.Destroy(temp);
	}

	private static void PlaceFaces() {
		List<Vector3> cornerVectors = new List<Vector3>();
		// first 3 vertices are the corners
		cornerVectors.Add(_templateCorners[0].transform.position);
		cornerVectors.Add(_templateCorners[1].transform.position);
		cornerVectors.Add(_templateCorners[2].transform.position);

		List<Vertex> uniqueEdges = new List<Vertex>();
		List<Vertex> uniqueCorners = new List<Vertex>();
		List<Vertex> connectedEdges = new List<Vertex>(); // list of most recent completely connected edges (used for mending the seams)
		
		// Copy and place each face and remove all resulting duplicates
		for (int i = 0; i < _geodesicCorners.Count; i++) {
			for (int j = 0; j < _geodesicCorners[i].neighbors.Count; j++) {
				// if index < i then face is already complete
				List<int> worldCornerIndex = new List<int>();
				int neighbor = _geodesicCorners.IndexOf(_geodesicCorners[i].neighbors[j]);
				int adjacentNeighbor = _geodesicCorners.IndexOf(_geodesicCorners[i].neighbors[(j + 1) % _geodesicCorners[i].neighbors.Count]);
				worldCornerIndex.Add(i);
				worldCornerIndex.Add(neighbor);
				worldCornerIndex.Add(adjacentNeighbor);
				if (neighbor > i && adjacentNeighbor > i) {
					// Place and orient the faces
					List<Vector3> worldCornerVectors = new List<Vector3>();
					worldCornerVectors.Add(_geodesicCorners[i].transform.position);
					worldCornerVectors.Add(_geodesicCorners[neighbor].transform.position);
					worldCornerVectors.Add(_geodesicCorners[adjacentNeighbor].transform.position);

					List<Vertex> newFace = _templateVertices;
					List<Vertex> newEdges = _templateEdges;
					List<Vertex> newCorners = _templateCorners;
					CopyFace(ref newFace, ref newEdges, ref newCorners);
					MoveFace(newFace, cornerVectors, worldCornerVectors);
					//print("Place Face: " + worldCornerIndex[0] + ", " + worldCornerIndex[1] + ", " + worldCornerIndex[2]);

					// rearrange newCorners to match worldCornerIndex order
					List<Vertex> tempNewCorners = new List<Vertex>();
					for (int m = 0; m < worldCornerIndex.Count; m++) {
						for (int n = 0; n < newCorners.Count; n++) {
							if (Vector3.Distance(newCorners[n].transform.position, _geodesicCorners[worldCornerIndex[m]].transform.position) < _edgeErrorValue) {
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
								// get # of real neighbors (not null)
								if (_faceLengthX == 0 || _faceLengthY == 0) { // edge vertices have four neighbors
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
						int indexNewNeighbor = _geodesicCorners[worldCornerIndex[n]].neighbors.IndexOf(_geodesicCorners[worldCornerIndex[(n + 1) % worldCornerIndex.Count]]);
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
							int indexDiff = indexNewNeighbor - indexFirst;
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
						Object.DestroyImmediate(toDestroy[k].gameObject);
					}
					worldTiles.AddRange(newFace);
				}
			}
		}
	}

	private static void AddFaceToWorld() {

	}
}