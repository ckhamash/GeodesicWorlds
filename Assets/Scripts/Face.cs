using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Face {
	private List<Vertex> _vertices, _edges, _corners;

	// Warning: very slow. Only use for initial face, create new ones using the Copy function.
	public Face(int dimensionA, int dimensionB, Vertex vertexPrefab) {
		_vertices = new List<Vertex>();
		_edges = new List<Vertex>();
		_corners = new List<Vertex>();
		int l = dimensionA + dimensionB + 1;
		// Create hexagonal grid of l by l vertices (bottom left is 0,0; next row is 0.5,sqrt(3)/2)
		for (int i = 0; i < l; i++) {
			for (int j = 0; j < l; j++) {
				Vertex vertex = Object.Instantiate(vertexPrefab);
				//vert.transform.parent = transform;
				vertex.transform.position = new Vector3(j + i * 0.5f, i * Mathf.Sqrt(3) / 2);
				//vert.maxNeighbors = 6;
				//for (int k = 0; k < vert.maxNeighbors; k++)
				for (int k = 0; k < 6; k++)
					vertex.neighbors.Add(null); // initialize neighbor list
				_vertices.Add(vertex);
			}
		}
		// define neighbors clockwise: 0 = i+1, 1 = i-length+1, 2 = i-length, 3 = i-1, 4 = i+length-1, 5 = i+length
		for (int i = 0; i < _vertices.Count; i++) {
			if (i % l != l - 1) // check right
				_vertices[i].neighbors[0] = _vertices[i + 1];
			if (i % l != l - 1 && i >= l) // check bottom right
				_vertices[i].neighbors[1] = _vertices[i - l + 1];
			if (i >= l) // check bottom left
				_vertices[i].neighbors[2] = _vertices[i - l];
			if (i % l != 0) // check left
				_vertices[i].neighbors[3] = _vertices[i - 1];
			if (i % l != 0 && i < _vertices.Count - l) // check top left
				_vertices[i].neighbors[4] = _vertices[i + l - 1];
			if (i < _vertices.Count - l) // check top right
				_vertices[i].neighbors[5] = _vertices[i + l];
		}

		// define the triangle corners: 0 = [b], 1 = [(a+b)length], 2 = [b+a + b*length]
		_corners.Add(_vertices[dimensionB]);
		_corners.Add(_vertices[l * (dimensionA + dimensionB)]);
		_corners.Add(_vertices[(dimensionA + dimensionB) + (l * dimensionB)]);
		_corners[0].GetComponent<Renderer>().material.color = Color.magenta;
		_corners[1].GetComponent<Renderer>().material.color = Color.magenta;
		_corners[2].GetComponent<Renderer>().material.color = Color.magenta;

		foreach (Vertex c in _corners) { // reduce neighbors to 5
												 //c.maxNeighbors = 5;
			for (int i = 0; i < c.neighbors.Count; i++) {
				if (c.neighbors[i] == null) {
					c.neighbors.RemoveAt(i);
					break;
				}
			}
		}

		// define points on edge  (these will be duplicates)
		int gcd = GreatestCommonDivisor(dimensionA, dimensionB);
		int a = dimensionA / gcd;
		int b = dimensionB / gcd;
		for (int i = 1; i < gcd; i++) {
			_edges.Add(_vertices[(dimensionB) + ((b * i * (l - 1)) + (a * i * l))]); // corner[0] --> corner[1]
			_edges.Add(_vertices[(l * (dimensionA + dimensionB)) + ((b * i) + (a * i * (-l + 1)))]); // corner[1] --> corner[2]
			_edges.Add(_vertices[(dimensionB) + ((b * i * l) + (a * i))]); // corner[0] --> corner[2]
		}
		// debug
		for (int i = 0; i < _edges.Count; i++)
			_edges[i].GetComponent<Renderer>().material.color = Color.cyan;

		// define triangle bounds
		GameObject temp = new GameObject();
		PolygonCollider2D triangleBounds = temp.AddComponent<PolygonCollider2D>();
		Vector2[] points = new Vector2[] { _corners[0].transform.position, _corners[1].transform.position, _corners[2].transform.position };
		triangleBounds.SetPath(0, points);


		float leftX = _corners[0].transform.position.x; // leftmost X val
		float rightX = _corners[0].transform.position.x; // rightmost X val
		for (int i = 0; i < _corners.Count; i++) {
			if (leftX > _corners[i].transform.position.x)
				leftX = _corners[i].transform.position.x;
			if (rightX < _corners[i].transform.position.x)
				rightX = _corners[i].transform.position.x;
		}

		// remove external neighbors (keep border positions), specify vertices that fall on triangle edge(these will be duplicates on the globe)
		List<Vertex> toRemove = new List<Vertex>();
		for (int i = 0; i < _vertices.Count; i++) {
			// trim obvious external (left and right)
			if (_vertices[i].transform.position.x < leftX || _vertices[i].transform.position.x > rightX) {
				toRemove.Add(_vertices[i]);
				continue;
			}

			if (_edges.Contains(_vertices[i]) || _corners.Contains(_vertices[i]))
				continue;

			if (!triangleBounds.OverlapPoint(_vertices[i].transform.position)) // outside of triangle
				toRemove.Add(_vertices[i]);
		}
		// delete temp object
		Object.Destroy(temp);

		for (int i = 0; i < toRemove.Count; i++) {
			_vertices.Remove(toRemove[i]);
			Object.DestroyImmediate(toRemove[i].gameObject);
		}
		//Debug
		_corners[0].GetComponent<Renderer>().material.color = Color.magenta;
		_corners[1].GetComponent<Renderer>().material.color = Color.magenta;
		_corners[2].GetComponent<Renderer>().material.color = Color.magenta;
	}

	// Create a copy of face with a new set of Vertex objects
	public Face(Face face, Vertex vertexPrefab) {
		_vertices = new List<Vertex>();
		_edges = new List<Vertex>();
		_corners = new List<Vertex>();

		// copy position and neighbor count
		for (int i = 0; i < face.Vertices.Count; i++) {
			Vertex copy = Object.Instantiate(vertexPrefab);
			copy.transform.position = face.Vertices[i].transform.position;
			//copy.maxNeighbors = verts[i].maxNeighbors;
			copy.GetComponent<Renderer>().material.color = face.Vertices[i].GetComponent<Renderer>().material.color; // debug
			_vertices.Add(copy);
			if (face.Edges.Contains(face.Vertices[i]))
				_edges.Add(copy);
			if (face.Corners.Contains(face.Vertices[i]))
				_corners.Add(copy);
		}

		// find neighbors by index
		for (int i = 0; i < face.Vertices.Count; i++) {
			for (int j = 0; j < face.Vertices[i].neighbors.Count; j++) {
				int index = face.Vertices.IndexOf(face.Vertices[i].neighbors[j]);
				if (index < 0) // null value
					_vertices[i].neighbors.Add(null);
				else
					_vertices[i].neighbors.Add(_vertices[index]);
			}
		}
	}

	public List<Vertex> Vertices {
		get {
			return _vertices;
		}
	}

	public List<Vertex> Edges {
		get {
			return _edges;
		}
	}

	public List<Vertex> Corners {
		get {
			return _corners;
		}
	}

	// Return the greatest common divisor of two integers
	private int GreatestCommonDivisor(int a, int b) {
		if (b == 0)
			return a;
		else
			return GreatestCommonDivisor(b, a % b);
	}
	
	/// <summary>
	/// Move this Face to a new position and orientation, where newPlane is 3 points defining a plane.
	/// The first corner of the Face will be moved to the first position defined in newPlane.
	/// 
	/// </summary>
	public void MoveTo(List<Vector3> newPlane) {
		// normVec = crossproduct(a-b,a-c) where abc are defined clockwise(following LHR such that normal points outward)
		// upVec = b - mid(a, c)
		// use Quaternion.LookRotation
		Vector3[] cornerPositions = { _corners[0].transform.position, _corners[1].transform.position, _corners[2].transform.position };

		// calculate copy orientation
		Vector3 norm = Vector3.Cross(cornerPositions[0] - cornerPositions[1], cornerPositions[0] - cornerPositions[2]);
		Vector3 up = cornerPositions[1] - (cornerPositions[0] + cornerPositions[2]) / 2;

		// create temp gameobj, give it the copy orientation and position of face corner a
		GameObject temp = new GameObject();
		temp.transform.rotation = Quaternion.LookRotation(norm, up);
		temp.transform.position = cornerPositions[0];

		// parent verts to temp object
		foreach (Vertex vertex in _vertices) {
			vertex.transform.position = temp.transform.InverseTransformPoint(vertex.transform.position);
			//vert.transform.SetParent(temp.transform);
		}

		// calculate new orientation
		norm = Vector3.Cross(newPlane[0] - newPlane[1], newPlane[0] - newPlane[2]);
		up = newPlane[1] - (newPlane[0] + newPlane[2]) / 2;

		// manipulate temp obj to new position and orientation
		temp.transform.rotation = Quaternion.LookRotation(norm, up);
		temp.transform.position = newPlane[0];

		// remove vertices' parenting
		foreach (Vertex vertex in _vertices) {
			vertex.transform.position = temp.transform.TransformPoint(vertex.transform.position);
			vertex.transform.SetParent(null);
		}

		Object.Destroy(temp);
	}

	// Move each vertex directly towards the surface of the sphere as defined by center and radius
	public void ConformToSphere(Vector3 center, float radius) {
		foreach(Vertex v in _vertices) {
			Vector3 relativePosition = v.transform.position - center;
			relativePosition *= radius / relativePosition.magnitude;
			v.transform.position = relativePosition + center;
		}
	}

	// compares edges and corners
	public void Join(Face face) {

	}
}
