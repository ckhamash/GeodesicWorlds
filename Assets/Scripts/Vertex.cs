using UnityEngine;
using System.Collections;
using System.Collections.Generic; // Lists

public class Vertex : MonoBehaviour {

	//public int maxNeighbors; // number of neighbors/sides this vertex should have
	public List<Vertex> neighbors; // list of neighbors, ordered clockwise

	// Use this for initialization
	void Start () {

	}
	
	// Update is called once per frame
	void Update () {

		// Debug: draw lines to neighbors
		for (int i = 0; i < neighbors.Count; i++)
			Debug.DrawLine(transform.position, getMidpoint(i), new Color(0, 0, 1));
	}

	// returns the position of the midpoint between this vertex to its neighbor
	Vector3 getMidpoint(int index) {
		if (index < neighbors.Count && neighbors[index] != null)
			return transform.position + ((neighbors[index].transform.position - transform.position) / 2);
		else
			return transform.position;
	}

	// getOppositeNeighbor
}
