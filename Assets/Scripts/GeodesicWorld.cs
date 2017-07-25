using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class GeodesicWorld : MonoBehaviour {

	public Vertex vertPrefab;
	//public Surface surface; // use for viewing in editor
	public List<Vertex> tiles;
	public int hexA, hexB; // defines the geodesic type
	public float radius;
	private float _edgeErrorValue = .001f;
	// Use this for initialization
	void Start() {
		tiles = GeodesicGenerator.GenerateGeodesic(hexA, hexB, vertPrefab);
		foreach (Vertex t in tiles) {
			t.transform.SetParent(transform);
		}
	}
	
	// Update is called once per frame
	void Update() {
	}
}
